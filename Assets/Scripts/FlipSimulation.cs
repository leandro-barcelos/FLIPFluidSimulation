using System.Threading.Tasks;
using UnityEngine;

public class FlipSimulation : MonoBehaviour
{

    const int NumThreads = 8;

    #region Auxiliary Structures

    enum CellType
    {
        Fluid,
        Solid,
        Air
    }

    private struct MeshProperties
    {
        // ReSharper disable once NotAccessedField.Local
        public Matrix4x4 Mat;
        // ReSharper disable once NotAccessedField.Local
        public Vector4 Color;

        public static int Size()
        {
            return
                sizeof(float) * 4 * 4 + // matrix;
                sizeof(float) * 4;      // color;
        }
    }

    #endregion

    #region Public Variables

    [Header("Simulation Parameters")]
    public float gravityAcceleration = -9.81f;
    [Range(0f, 1f)] public float flipRatio = 0.9f;
    public int numPressureIterations = 50;
    public int numParticleIterations = 2;
    public bool compensateDrift = true;
    public bool separateParticles = true;
    public float restDensity = 0.0f;
    public float collisionForce = 5.0f;

    [Header("Visualization")]
    public bool showParticles = true;
    public bool showGrid = false;

    [Header("Fluid")]
    public float density = 1000.0f;
    public float waterHeight = 0.8f;
    public float waterWidth = 0.6f;
    public float waterDepth = 0.6f;

    [Header("Particles")]
    public ComputeShader integrateParticlesShader;
    public ComputeShader handleParticleCollisionsShader;

    [Header("Grid")]
    public int gridResolution = 100;

    [Header("Rendering")]
    public float occlusionRange;
    public Material particleMaterial;
    public ComputeShader updateParticlePropertiesShader;

    [Header("Obstacles")]
    public ComputeShader generateObstaclesSDFShader;
    public int FastSweepingPasses = 3;

    #endregion

    #region Private Variables

    private float _deltaTime;
    private Vector3 _simOrigin;

    //  Grid
    private float _cellSpacing, _cellInvSpacing;
    private Vector3Int _gridDimensions;
    private int _numCells;
    private float[] _cellU, _cellV, _cellW, _cellRu, _cellRv, _cellRw, _cellPrevU, _cellPrevV, _cellPrevW, _cellS;
    private ComputeBuffer _cellType;
    private Color[] _cellColor;

    // Particles
    private int _maxParticles, _particleNumCells;
    private float _particleRadius, _particleInvSpacing;
    private Vector3Int _particleCellDimensions;
    private ComputeBuffer _particlePos, _particleVel;
    private float[] _particleDensity;

    // Rendering
    private Mesh _particleMesh;
    private ComputeBuffer _meshPropertiesBuffer, _argsBuffer;
    private Bounds _bounds;

    // Obstacles
    // private List<MeshFilter> _obstacles;
    private ComputeBuffer _obstaclesSDF;

    #endregion

    void Start()
    {
        _simOrigin = transform.position - transform.localScale / 2;
        _deltaTime = Time.fixedDeltaTime;
        // _obstacles = new List<MeshFilter>();

        InitializeRendering();
        InitializeGrid();
        InitializeParticles();

    }

    void Update()
    {
        IntegrateParticles();

        if (separateParticles)
        {
            // PushParticlesApart();
        }

        foreach (MeshFilter meshFilter in FindObjectsByType<MeshFilter>(FindObjectsSortMode.None))
        {
            if (meshFilter.gameObject.name == "Simulation") continue;
            GenerateObstaclesSDF(meshFilter);
        }

        HandleParticleCollisions(); // TODO: handle obstacle collisions

        // TODO: update cell types
        // TODO: transfer velocities to grid
        // TODO: update particle density
        // TODO: solve incompressibility

        // TransferVelocityToParticle();

        // TODO: update grid properties

        UpdateMeshProperties(); // TODO: Update colors as well

        if (showParticles)
            Graphics.DrawMeshInstancedIndirect(_particleMesh, 0, particleMaterial, _bounds, _argsBuffer);
    }

    private void OnDisable()
    {
        _meshPropertiesBuffer?.Release();
        _argsBuffer?.Release();
        _particlePos?.Release();
        _particleVel?.Release();
        _obstaclesSDF?.Release();
    }

    #region Initialization

    private void InitializeRendering()
    {
        _particleMesh = Resources.GetBuiltinResource<Mesh>("Cube.fbx");
        _bounds = new Bounds(transform.position, Vector3.one * (occlusionRange + 1));
    }

    private void InitializeGrid()
    {
        var scale = transform.localScale;
        _cellSpacing = scale.y / gridResolution;

        _gridDimensions = new Vector3Int(
            Mathf.FloorToInt(scale.x / _cellSpacing) + 1,
            Mathf.FloorToInt(scale.y / _cellSpacing) + 1,
            Mathf.FloorToInt(scale.z / _cellSpacing) + 1
        );

        _cellSpacing = Mathf.Max(scale.x / _gridDimensions.x, scale.y / _gridDimensions.y, scale.z / _gridDimensions.z);
        _cellInvSpacing = 1.0f / _cellSpacing;
        _numCells = _gridDimensions.x * _gridDimensions.y * _gridDimensions.z;

        _cellU = new float[_numCells];
        _cellV = new float[_numCells];
        _cellW = new float[_numCells];
        _cellRu = new float[_numCells];
        _cellRv = new float[_numCells];
        _cellRw = new float[_numCells];
        _cellS = new float[_numCells];
        _cellType = new ComputeBuffer(_numCells, sizeof(CellType));
        _cellColor = new Color[_numCells];

        for (int i = 0; i < _gridDimensions.x; i++)
        {
            for (int j = 0; j < _gridDimensions.y; j++)
            {
                for (int k = 0; k < _gridDimensions.z; k++)
                {
                    float s = 1.0f;
                    if (i == 0 || j == 0 || k == 0 || i == _gridDimensions.x - 1 || k == _gridDimensions.z - 1)
                    {
                        s = 0.0f;
                    }
                    _cellS[i * _gridDimensions.y * _gridDimensions.z + j * _gridDimensions.z + k] = s;
                }
            }
        }
    }

    private void InitializeParticles()
    {
        _particleRadius = 0.3f * _cellSpacing;
        var scale = transform.localScale;
        _particleInvSpacing = 1f / (2.2f * _particleRadius);

        _particleCellDimensions = new Vector3Int(
            Mathf.FloorToInt(waterWidth * scale.x * _particleInvSpacing) + 1,
            Mathf.FloorToInt(waterHeight * scale.y * _particleInvSpacing) + 1,
            Mathf.FloorToInt(waterDepth * scale.z * _particleInvSpacing) + 1
        );

        _particleNumCells = _particleCellDimensions.x * _particleCellDimensions.y * _particleCellDimensions.z;

        var numX = Mathf.FloorToInt((waterWidth * scale.x - 2.0f * _cellSpacing - 2.0f * _particleRadius) / _particleRadius * 2);
        var numY = Mathf.FloorToInt((waterHeight * scale.y - 2.0f * _cellSpacing - 2.0f * _particleRadius) / _particleRadius * 2);
        var numZ = Mathf.FloorToInt((waterDepth * scale.z - 2.0f * _cellSpacing - 2.0f * _particleRadius) / _particleRadius * 2);
        _maxParticles = numX * numY * numZ;

        _particlePos = new ComputeBuffer(_maxParticles, sizeof(float) * 3);
        _particleVel = new ComputeBuffer(_maxParticles, sizeof(float) * 3);
        _particleDensity = new float[_maxParticles];

        InitializeParticlePositions();
    }

    private void InitializeParticlePositions()
    {
        uint[] args = { 0, 0, 0, 0, 0 };
        args[0] = _particleMesh.GetIndexCount(0);
        args[1] = (uint)_maxParticles;
        args[2] = _particleMesh.GetIndexStart(0);
        args[3] = _particleMesh.GetBaseVertex(0);

        _argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _argsBuffer.SetData(args);

        _meshPropertiesBuffer = new ComputeBuffer(_maxParticles, MeshProperties.Size());

        MeshProperties[] properties = new MeshProperties[_maxParticles];

        Vector3 particleScale = new Vector3(_particleRadius * 2, _particleRadius * 2, _particleRadius * 2);
        Quaternion rotation = Quaternion.identity;
        Vector3 offset = transform.localScale * -0.5f;

        var particlePos = new Vector3[_maxParticles];

        int index = 0;
        for (int i = 0; i < _particleCellDimensions.x; i++)
        {
            for (int j = 0; j < _particleCellDimensions.y; j++)
            {
                for (int k = 0; k < _particleCellDimensions.z; k++)
                {
                    float x = _particleRadius + i * (2.0f * _particleRadius) + offset.x;
                    float y = _particleRadius + j * (2.0f * _particleRadius) + offset.y;
                    float z = _particleRadius + k * (2.0f * _particleRadius) + offset.z;

                    particlePos[index] = new Vector3(x, y, z);

                    MeshProperties props = new MeshProperties
                    {
                        Mat = Matrix4x4.TRS(particlePos[index], rotation, particleScale),
                        Color = Color.white
                    };

                    properties[index] = props;
                    index++;
                }
            }
        }

        _particlePos.SetData(particlePos);
        _meshPropertiesBuffer.SetData(properties);
        particleMaterial.SetBuffer(ShaderIDs.Properties, _meshPropertiesBuffer);
    }

    #endregion

    #region Simulation

    private void IntegrateParticles()
    {
        integrateParticlesShader.SetFloat("_TimeStep", Time.fixedDeltaTime);
        integrateParticlesShader.SetFloat("_Gravity", gravityAcceleration);

        integrateParticlesShader.SetVector("_ParticleSize", (Vector3)_particleCellDimensions);

        integrateParticlesShader.SetBuffer(0, "_ParticlePos", _particlePos);
        integrateParticlesShader.SetBuffer(0, "_ParticleVel", _particleVel);

        integrateParticlesShader.Dispatch(0, Mathf.CeilToInt((float)_particleCellDimensions.x / NumThreads), Mathf.CeilToInt((float)_particleCellDimensions.y / NumThreads), Mathf.CeilToInt((float)_particleCellDimensions.z / NumThreads));
    }

    private void UpdateMeshProperties()
    {

        updateParticlePropertiesShader.SetVector(ShaderIDs.ParticleSize, (Vector3)_particleCellDimensions);
        updateParticlePropertiesShader.SetBuffer(0, ShaderIDs.ParticlePos, _particlePos);
        updateParticlePropertiesShader.SetBuffer(0, ShaderIDs.Properties, _meshPropertiesBuffer);

        updateParticlePropertiesShader.Dispatch(0,
            Mathf.CeilToInt((float)_particleCellDimensions.x / NumThreads),
            Mathf.CeilToInt((float)_particleCellDimensions.y / NumThreads),
            Mathf.CeilToInt((float)_particleCellDimensions.z / NumThreads));
    }

    private void HandleParticleCollisions()
    {
        handleParticleCollisionsShader.SetFloat(ShaderIDs.ParticleRadius, _particleRadius);
        handleParticleCollisionsShader.SetFloat(ShaderIDs.CellSpacing, _cellSpacing);
        handleParticleCollisionsShader.SetFloat(ShaderIDs.CellInvSpacing, _cellInvSpacing);
        handleParticleCollisionsShader.SetFloat(ShaderIDs.TimeStep, Time.fixedDeltaTime);

        handleParticleCollisionsShader.SetVector(ShaderIDs.ParticleSize, (Vector3)_particleCellDimensions);
        handleParticleCollisionsShader.SetVector(ShaderIDs.SimOrigin, _simOrigin);
        handleParticleCollisionsShader.SetVector(ShaderIDs.GridDimensions, (Vector3)_gridDimensions);

        handleParticleCollisionsShader.SetBuffer(0, ShaderIDs.ParticlePos, _particlePos);
        handleParticleCollisionsShader.SetBuffer(0, ShaderIDs.ParticleVel, _particleVel);
        handleParticleCollisionsShader.SetBuffer(0, ShaderIDs.CellType, _cellType);
        handleParticleCollisionsShader.SetBuffer(0, ShaderIDs.ObstaclesSDF, _obstaclesSDF);

        handleParticleCollisionsShader.Dispatch(0, Mathf.CeilToInt((float)_particleCellDimensions.x / NumThreads), Mathf.CeilToInt((float)_particleCellDimensions.y / NumThreads), Mathf.CeilToInt((float)_particleCellDimensions.z / NumThreads));
    }

    private void GenerateObstaclesSDF(MeshFilter meshFilter)
    {
        _obstaclesSDF = new ComputeBuffer(_numCells, sizeof(float));
        ComputeBuffer closestTriangles = new(_numCells, sizeof(int));
        ComputeBuffer intersectionCounts = new(_numCells, sizeof(int));
        int triangleCount = meshFilter.mesh.triangles.Length;
        ComputeBuffer triangles = new(triangleCount, sizeof(int));
        triangles.SetData(meshFilter.mesh.triangles);
        ComputeBuffer vertices = new(meshFilter.mesh.vertexCount, sizeof(float) * 3);
        vertices.SetData(meshFilter.mesh.vertices);

        InitSDFNearGeometry(meshFilter.transform.localToWorldMatrix, closestTriangles, intersectionCounts, triangleCount, triangles, vertices);

        FastSweeping(closestTriangles, vertices, triangles);

        CalculateSDFSign(intersectionCounts);

        closestTriangles.Release();
        intersectionCounts.Release();
        triangles.Release();
        vertices.Release();
    }

    private void CalculateSDFSign(ComputeBuffer intersectionCounts)
    {
        int signKernel = generateObstaclesSDFShader.FindKernel("CalculateSDFSign");

        generateObstaclesSDFShader.SetVector(ShaderIDs.GridDimensions, (Vector3)_gridDimensions);

        generateObstaclesSDFShader.SetBuffer(signKernel, ShaderIDs.ObstaclesSDF, _obstaclesSDF);
        generateObstaclesSDFShader.SetBuffer(signKernel, ShaderIDs.CellType, _cellType);
        generateObstaclesSDFShader.SetBuffer(signKernel, ShaderIDs.IntersectionCounts, intersectionCounts);

        generateObstaclesSDFShader.Dispatch(signKernel, 1, Mathf.CeilToInt((float)_gridDimensions.y / NumThreads), Mathf.CeilToInt((float)_gridDimensions.z / NumThreads));
    }

    private void InitSDFNearGeometry(Matrix4x4 transform, ComputeBuffer closestTriangles, ComputeBuffer intersectionCounts, int triangleCount, ComputeBuffer triangles, ComputeBuffer vertices)
    {
        int initKernel = generateObstaclesSDFShader.FindKernel("Init");
        generateObstaclesSDFShader.SetMatrix(ShaderIDs.Transform, transform);
        generateObstaclesSDFShader.SetVector(ShaderIDs.SimOrigin, _simOrigin);
        generateObstaclesSDFShader.SetVector(ShaderIDs.GridDimensions, (Vector3)_gridDimensions);
        generateObstaclesSDFShader.SetFloat(ShaderIDs.CellSpacing, _cellSpacing);
        generateObstaclesSDFShader.SetInt(ShaderIDs.TriangleCount, triangleCount);

        generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.CellType, _cellType);
        generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.ObstaclesSDF, _obstaclesSDF);
        generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.ClosestTriangles, closestTriangles);
        generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.IntersectionCounts, intersectionCounts);
        generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.Vertices, vertices);
        generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.Triangles, triangles);

        generateObstaclesSDFShader.Dispatch(initKernel, Mathf.CeilToInt((float)_gridDimensions.x / NumThreads), Mathf.CeilToInt((float)_gridDimensions.y / NumThreads), Mathf.CeilToInt((float)_gridDimensions.z / NumThreads));
                    }

    private void FastSweeping(ComputeBuffer closestTriangles, ComputeBuffer vertices, ComputeBuffer triangles)
    {
        int fastSweepKernel = generateObstaclesSDFShader.FindKernel("FastSweeping");

        // Set parameters
        generateObstaclesSDFShader.SetVector(ShaderIDs.SimOrigin, _simOrigin);
        generateObstaclesSDFShader.SetVector(ShaderIDs.GridDimensions, (Vector3)_gridDimensions);
        generateObstaclesSDFShader.SetFloat(ShaderIDs.CellSpacing, _cellSpacing);

        // Set buffers
        generateObstaclesSDFShader.SetBuffer(fastSweepKernel, ShaderIDs.ObstaclesSDF, _obstaclesSDF);
        generateObstaclesSDFShader.SetBuffer(fastSweepKernel, ShaderIDs.ClosestTriangles, closestTriangles);
        generateObstaclesSDFShader.SetBuffer(fastSweepKernel, ShaderIDs.Vertices, vertices);
        generateObstaclesSDFShader.SetBuffer(fastSweepKernel, ShaderIDs.Triangles, triangles);

        int[] iDirections = { 1, -1, 1, -1, 1, -1, 1, -1 };
        int[] jDirections = { 1, 1, 1, 1, -1, -1, -1, -1 };
        int[] kDirections = { 1, -1, -1, 1, 1, -1, -1, 1 };

        // Dispatch in 8 directions
        for (int _ = 0; _ < FastSweepingPasses; _++)
        {
            for (int dir = 0; dir < 8; dir++)
        {
                Vector3 direction = new(iDirections[dir], jDirections[dir], kDirections[dir]);

                // Set directional parameters
                generateObstaclesSDFShader.SetVector(ShaderIDs.Direction, direction);

                // Dispatch compute shader for current direction
                generateObstaclesSDFShader.Dispatch(fastSweepKernel, Mathf.CeilToInt((float)_gridDimensions.x / NumThreads), Mathf.CeilToInt((float)_gridDimensions.y / NumThreads), Mathf.CeilToInt((float)_gridDimensions.z / NumThreads));
            }
        }
    }

    // private void PushParticlesApart()
    // {
    //     ConcurrentDictionary<int, ConcurrentBag<int>> particlesInCell = new();

    //     // Populate the dictionary with particles in cells
    //     Parallel.For(0, _maxParticles, i =>
    //     {
    //         var cellIndex = ParticleToCell1D(_particlePos[i]);
    //         lock (particlesInCell)
    //         {
    //             if (!particlesInCell.TryGetValue(cellIndex, out var particles))
    //             {
    //                 particles = new ConcurrentBag<int>();
    //                 particlesInCell[cellIndex] = particles;
    //             }
    //             particles.Add(i);
    //         }
    //     });

    //     var maxDist = 2f * _particleRadius;
    //     var maxDist2 = maxDist * maxDist;

    //     for (var iter = 0; iter < numParticleIterations; iter++)
    //     {
    //         Parallel.For(0, _maxParticles, i =>
    //         {
    //             var cellPos = ParticleToCell3D(_particlePos[i]);

    //             var x0 = Mathf.Max(cellPos[0] - 1, 0);
    //             var y0 = Mathf.Max(cellPos[1] - 1, 0);
    //             var z0 = Mathf.Max(cellPos[2] - 1, 0);
    //             var x1 = Mathf.Min(cellPos[0] + 1, _gridDimensions.x - 1);
    //             var y1 = Mathf.Min(cellPos[1] + 1, _gridDimensions.y - 1);
    //             var z1 = Mathf.Min(cellPos[2] + 1, _gridDimensions.z - 1);

    //             for (var xi = x0; xi <= x1; xi++)
    //             {
    //                 for (var yi = y0; yi <= y1; yi++)
    //                 {
    //                     for (var zi = z0; zi <= z1; zi++)
    //                     {
    //                         var cellNr = To1D(xi, yi, zi);
    //                         if (particlesInCell.TryGetValue(cellNr, out var particles))
    //                         {
    //                             foreach (int particleIndex in particles)
    //                             {
    //                                 if (particleIndex == i) continue;

    //                                 var posNeighbor = _particlePos[particleIndex];
    //                                 var particleToNeighbor = posNeighbor - _particlePos[i];
    //                                 var distanceSquared = particleToNeighbor.sqrMagnitude;
    //                                 if (distanceSquared > maxDist2 || distanceSquared == 0.0f) continue;

    //                                 var distance = Mathf.Sqrt(distanceSquared);
    //                                 var overlap = 0.5f * (maxDist - distance) / distance;
    //                                 particleToNeighbor *= overlap;

    //                                 lock (_particlePos)
    //                                 {
    //                                     _particlePos[i] -= particleToNeighbor;
    //                                     _particlePos[particleIndex] += particleToNeighbor;
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         });
    //     }
    // }


    #endregion

    #region Utils

    private int ParticleToCell1D(Vector3 position)
    {
        var xCell = Mathf.FloorToInt((position.x - _simOrigin.x) * _cellInvSpacing);
        var yCell = Mathf.FloorToInt((position.y - _simOrigin.y) * _cellInvSpacing);
        var zCell = Mathf.FloorToInt((position.z - _simOrigin.z) * _cellInvSpacing);

        xCell = Mathf.Clamp(xCell, 0, _gridDimensions.x - 1);
        yCell = Mathf.Clamp(yCell, 0, _gridDimensions.y - 1);
        zCell = Mathf.Clamp(zCell, 0, _gridDimensions.z - 1);

        return To1D(xCell, yCell, zCell);
    }

    private int[] ParticleToCell3D(Vector3 position)
    {
        var xCell = Mathf.FloorToInt((position.x - _simOrigin.x) * _cellInvSpacing);
        var yCell = Mathf.FloorToInt((position.y - _simOrigin.y) * _cellInvSpacing);
        var zCell = Mathf.FloorToInt((position.z - _simOrigin.z) * _cellInvSpacing);

        return new int[] { xCell, yCell, zCell };
    }

    private bool IsValidIndex(int i, int j, int k)
    {
        return i >= 0 && i < _gridDimensions.x &&
               j >= 0 && j < _gridDimensions.y &&
               k >= 0 && k < _gridDimensions.z;
    }

    private int To1D(int x, int y, int z)
    {
        return (z * _gridDimensions.x * _gridDimensions.y) + (y * _gridDimensions.x) + x;
    }

    private int[] To3D(int idx)
    {
        int z = idx / (_gridDimensions.x * _gridDimensions.y);
        idx -= (z * _gridDimensions.x * _gridDimensions.y);
        int y = idx / _gridDimensions.x;
        int x = idx % _gridDimensions.x;
        return new[] { x, y, z };
    }

    #endregion
}
