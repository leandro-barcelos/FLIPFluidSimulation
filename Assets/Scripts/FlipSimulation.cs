using System.Threading.Tasks;
using UnityEngine;

public class FlipSimulation : MonoBehaviour
{
    private static readonly int Properties = Shader.PropertyToID("_Properties");
    private static readonly int ParticleSize = Shader.PropertyToID("_ParticleSize");
    private static readonly int ParticlePos = Shader.PropertyToID("_ParticlePos");
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
    public int binarySearchSteps = 5;

    [Header("Visualization")]
    public bool showParticles = true;
    public bool showGrid = false;

    [Header("Fluid")]
    public float density = 1000.0f;
    public float waterHeight = 0.8f;
    public float waterWidth = 0.6f;
    public float waterDepth = 0.6f;

    [Header("Grid")]
    public int gridResolution = 100;

    [Header("Rendering")]
    public float occlusionRange;
    public Material particleMaterial;
    public ComputeShader updateParticlePropertiesShader;

    #endregion

    #region Private Variables

    private float _deltaTime;

    //  Grid
    private float _cellSpacing, _cellInvSpacing;
    private Vector3Int _gridDimensions;
    private int _numCells;
    private float[] _cellU, _cellV, _cellW, _cellRu, _cellRv, _cellRw, _cellPrevU, _cellPrevV, _cellPrevW, _cellS;
    private CellType[] _cellType;
    private Color[] _cellColor;

    // Particles
    private int _maxParticles, _particleNumCells;
    private float _particleRadius, _particleInvSpacing;
    private Vector3Int _particleCellDimensions;
    private Vector3[] _particlePos, _particleVel;
    private float[] _particleDensity;
    private int[] _numParticlesInCell, _firstParticleInCell, _particleCellIDs;

    // Rendering
    private Mesh _particleMesh;
    private ComputeBuffer _meshPropertiesBuffer, _argsBuffer;
    private Bounds _bounds;

    // Obstacles
    // private List<MeshFilter> _obstacles;
    private float[] _obstaclesSDF;

    #endregion

    void Start()
    {
        _deltaTime = Time.fixedDeltaTime;
        // _obstacles = new List<MeshFilter>();

        InitializeRendering();
        InitializeGrid();
        InitializeParticles();

        foreach (MeshFilter meshFilter in FindObjectsByType<MeshFilter>(FindObjectsSortMode.None))
        {
            if (meshFilter.gameObject.name == "Simulation") continue;
            GenerateObstaclesSDF(meshFilter);
        }
    }

    void Update()
    {
        IntegrateParticles();

        // TODO: push particles apart

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
        _cellType = new CellType[_numCells];
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
        _numParticlesInCell = new int[_particleNumCells];
        _firstParticleInCell = new int[_particleNumCells + 1];

        var numX = Mathf.FloorToInt((waterWidth * scale.x - 2.0f * _cellSpacing - 2.0f * _particleRadius) / _particleRadius * 2);
        var numY = Mathf.FloorToInt((waterHeight * scale.y - 2.0f * _cellSpacing - 2.0f * _particleRadius) / _particleRadius * 2);
        var numZ = Mathf.FloorToInt((waterDepth * scale.z - 2.0f * _cellSpacing - 2.0f * _particleRadius) / _particleRadius * 2);
        _maxParticles = numX * numY * numZ;

        _particleCellIDs = new int[_maxParticles];
        _particlePos = new Vector3[_maxParticles];
        _particleVel = new Vector3[_maxParticles];
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

                    _particlePos[index] = new Vector3(x, y, z);

                    MeshProperties props = new MeshProperties
                    {
                        Mat = Matrix4x4.TRS(_particlePos[index], rotation, particleScale),
                        Color = Color.white
                    };

                    properties[index] = props;
                    index++;
                }
            }
        }

        _meshPropertiesBuffer.SetData(properties);
        particleMaterial.SetBuffer(Properties, _meshPropertiesBuffer);
    }

    #endregion

    #region Simulation

    private void IntegrateParticles()
    {
        Parallel.For(0, _maxParticles, i =>
        {
            _particleVel[i].y += _deltaTime * gravityAcceleration;
            _particlePos[i] += _particleVel[i] * _deltaTime;
        });
    }

    private void UpdateMeshProperties()
    {
        var particlePosBuffer = new ComputeBuffer(_maxParticles, sizeof(float) * 3);
        particlePosBuffer.SetData(_particlePos);

        updateParticlePropertiesShader.SetVector(ParticleSize, (Vector3)_particleCellDimensions);
        updateParticlePropertiesShader.SetBuffer(0, ParticlePos, particlePosBuffer);
        updateParticlePropertiesShader.SetBuffer(0, Properties, _meshPropertiesBuffer);

        updateParticlePropertiesShader.Dispatch(0,
            Mathf.CeilToInt((float)_particleCellDimensions.x / NumThreads),
            Mathf.CeilToInt((float)_particleCellDimensions.y / NumThreads),
            Mathf.CeilToInt((float)_particleCellDimensions.z / NumThreads));

        particlePosBuffer.Release();
    }

    private void HandleParticleCollisions()
    {
        var simOrigin = transform.position - transform.localScale / 2;

        var minX = simOrigin.x + _cellSpacing + _particleRadius;
        var maxX = simOrigin.x + (_gridDimensions.x - 1) * _cellSpacing - _particleRadius;
        var minY = simOrigin.y + _cellSpacing + _particleRadius;
        var maxY = simOrigin.y + (_gridDimensions.y - 1) * _cellSpacing - _particleRadius;
        var minZ = simOrigin.z + _cellSpacing + _particleRadius;
        var maxZ = simOrigin.z + (_gridDimensions.z - 1) * _cellSpacing - _particleRadius;

        Parallel.For(0, _maxParticles, i =>
        {
            var xCell = Mathf.FloorToInt((_particlePos[i].x - simOrigin.x) * _cellInvSpacing);
            var yCell = Mathf.FloorToInt((_particlePos[i].y - simOrigin.y) * _cellInvSpacing);
            var zCell = Mathf.FloorToInt((_particlePos[i].z - simOrigin.z) * _cellInvSpacing);

            int cellIndex = To1D(xCell, yCell, zCell);

            if (_cellType[cellIndex] == CellType.Solid)
            {
                _particlePos[i] -= _particleVel[i] * _deltaTime;

                var start = Vector3.zero;
                var end = _particleVel[i];

                for (int j = 0; j < binarySearchSteps; j++)
                {
                    _particleVel[i] = (start + end) * 0.5f;
                    _particlePos[i] += _particleVel[i] * _deltaTime;

                    xCell = Mathf.FloorToInt((_particlePos[i].x - simOrigin.x) * _cellInvSpacing);
                    yCell = Mathf.FloorToInt((_particlePos[i].y - simOrigin.y) * _cellInvSpacing);
                    zCell = Mathf.FloorToInt((_particlePos[i].z - simOrigin.z) * _cellInvSpacing);

                    cellIndex = To1D(xCell, yCell, zCell);

                    if (_cellType[cellIndex] == CellType.Solid)
                    {
                        end = _particleVel[i];
                    }
                    else
                    {
                        start = _particleVel[i];
                    }

                    _particlePos[i] -= _particleVel[i] * _deltaTime;
                }

                if (_cellType[cellIndex] == CellType.Solid)
                {
                    _particlePos[i] += new Vector3(_obstaclesSDF[cellIndex], _obstaclesSDF[cellIndex], _obstaclesSDF[cellIndex]);
                }
            }

            if (_particlePos[i].x < minX)
            {
                _particlePos[i].x = minX;
                _particleVel[i].x = 0f;
            }
            else if (_particlePos[i].x > maxX)
            {
                _particlePos[i].x = maxX;
                _particleVel[i].x = 0f;
            }

            if (_particlePos[i].y < minY)
            {
                _particlePos[i].y = minY;
                _particleVel[i].y = 0f;
            }
            else if (_particlePos[i].y > maxY)
            {
                _particlePos[i].y = maxY;
                _particleVel[i].y = 0f;
            }

            if (_particlePos[i].z < minZ)
            {
                _particlePos[i].z = minZ;
                _particleVel[i].z = 0f;
            }
            else if (_particlePos[i].z > maxZ)
            {
                _particlePos[i].z = maxZ;
                _particleVel[i].z = 0f;
            }
        });
    }

    private void GenerateObstaclesSDF(MeshFilter meshFilter)
    {
        var simOrigin = transform.position - transform.localScale / 2;
        _obstaclesSDF = new float[_numCells];
        int[] closestTriangles = new int[_numCells];
        int[] intersectionCounts = new int[_numCells];

        Parallel.For(0, _numCells, i =>
        {
            _obstaclesSDF[i] = float.MaxValue;
            closestTriangles[i] = -1;
            intersectionCounts[i] = 0;
        });

        for (var e = 0; e < meshFilter.mesh.triangles.Length; e += 3)
        {
            var a = meshFilter.transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[e]]);
            var b = meshFilter.transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[e + 1]]);
            var c = meshFilter.transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[e + 2]]);

            for (int i = 0; i < _gridDimensions.x; i++)
            {
                for (int j = 0; j < _gridDimensions.y; j++)
                {
                    for (int k = 0; k < _gridDimensions.z; k++)
                    {
                        var index = To1D(i, j, k);
                        var p = simOrigin + new Vector3(i * _cellSpacing, j * _cellSpacing, k * _cellSpacing);
                        var q = simOrigin + new Vector3((i + 1) * _cellSpacing, j * _cellSpacing, k * _cellSpacing);

                        if (!RaycastTriangle(p, q, a, b, c, out var intersection)) continue;

                        intersectionCounts[index]++;

                        var d = Vector3.Distance(p, intersection);

                        if (d >= _obstaclesSDF[index]) continue;

                        _obstaclesSDF[index] = d;
                        closestTriangles[index] = e;
                    }
                }
            }
        }

        FastSweeping(meshFilter, _obstaclesSDF, closestTriangles);


        for (int j = 0; j < _gridDimensions.y; j++)
        {
            for (int k = 0; k < _gridDimensions.z; k++)
            {
                int c = 0;
                for (int i = 0; i < _gridDimensions.x; i++)
                {
                    int index = To1D(i, j, k);

                    if (c % 2 != 0)
                    {
                        _obstaclesSDF[index] *= -1;
                        _cellType[index] = CellType.Solid;
                    }

                    c += intersectionCounts[index];
                }
            }
        }

    }

    private void FastSweeping(MeshFilter meshFilter, float[] distances, int[] closestTriangles)
    {
        int[] iDirections = { 1, -1, 1, -1, 1, -1, 1, -1 };
        int[] jDirections = { 1, 1, 1, 1, -1, -1, -1, -1 };
        int[] kDirections = { 1, -1, -1, 1, 1, -1, -1, 1 };

        for (var dir = 0; dir < 8; dir++)
        {
            var iDir = iDirections[dir];
            var jDir = jDirections[dir];
            var kDir = kDirections[dir];

            for (int i = iDir > 0 ? 0 : _gridDimensions.x - 1;
                         iDir > 0 ? i < _gridDimensions.x : i >= 0;
                         i += iDir)
            {
                for (int j = jDir > 0 ? 0 : _gridDimensions.y - 1;
                             jDir > 0 ? j < _gridDimensions.y : j >= 0;
                             j += jDir)
                {
                    for (int k = kDir > 0 ? 0 : _gridDimensions.z - 1;
                                 kDir > 0 ? k < _gridDimensions.z : k >= 0;
                                 k += kDir)
                    {
                        UpdateDistanceAndClosestTriangle(meshFilter, distances, closestTriangles, i, j, k);
                    }
                }
            }
        }
    }

    private void UpdateDistanceAndClosestTriangle(MeshFilter meshFilter, float[] distances, int[] closestTriangles, int i, int j, int k)
    {
        var simOrigin = transform.position - transform.localScale / 2;
        float minDist = distances[To1D(i, j, k)];
        int closestTriangle = closestTriangles[To1D(i, j, k)];
        var p = simOrigin + new Vector3(i * _cellSpacing, j * _cellSpacing, k * _cellSpacing);

        // Iterate over neighbors
        int[] di = { -1, 1, 0, 0, 0, 0 };
        int[] dj = { 0, 0, -1, 1, 0, 0 };
        int[] dk = { 0, 0, 0, 0, -1, 1 };

        for (int n = 0; n < 6; n++)
        {
            int ni = i + di[n];
            int nj = j + dj[n];
            int nk = k + dk[n];

            if (IsValidIndex(ni, nj, nk))
            {
                int neighborIndex = To1D(ni, nj, nk);
                int triangle = closestTriangles[neighborIndex];
                if (triangle != -1)
                {
                    var v0 = meshFilter.transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[triangle]]);
                    var v1 = meshFilter.transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[triangle + 1]]);
                    var v2 = meshFilter.transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[triangle + 2]]);
                    float dist = CalculateDistanceToTriangle(p, v0, v1, v2);
                    if (dist < minDist)
                    {
                        minDist = dist;
                        closestTriangle = triangle;
                    }
                }
            }
        }

        distances[To1D(i, j, k)] = minDist;
        closestTriangles[To1D(i, j, k)] = closestTriangle;
    }

    #endregion

    #region Utils

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

    private float CalculateDistanceToTriangle(Vector3 p, Vector3 v0, Vector3 v1, Vector3 v2)
    {
        var v = p - v0;

        var n = CalculateSurfaceNormal(v0, v1, v2);

        return Mathf.Abs(Vector3.Dot(v, n));
    }

    private bool RaycastTriangle(Vector3 s, Vector3 e, Vector3 a, Vector3 b, Vector3 c, out Vector3 intersection)
    {
        intersection = Vector3.zero;
        var direction = e - s;
        var normalizedDirection = direction.normalized;

        Vector3 edge1 = b - a;
        Vector3 edge2 = c - a;
        Vector3 normal = Vector3.Cross(normalizedDirection, edge2);

        float dotProduct = Vector3.Dot(edge1, normal);

        if (dotProduct > -float.Epsilon && dotProduct < float.Epsilon) return false;

        float inverseDotProduct = 1.0f / dotProduct;
        Vector3 toStart = s - a;
        float trianglParam = inverseDotProduct * Vector3.Dot(toStart, normal);

        if (trianglParam is < 0f or > 1f) return false;
        Vector3 edgeCross = Vector3.Cross(toStart, edge1);
        float raycastParam = inverseDotProduct * Vector3.Dot(normalizedDirection, edgeCross);

        if (raycastParam < 0f || trianglParam + raycastParam > 1f) return false;

        float distance = inverseDotProduct * Vector3.Dot(edge2, edgeCross);
        if (distance > float.Epsilon && distance < direction.magnitude)
        {
            intersection = s + normalizedDirection * distance;
            return true;
        }

        return false;
    }

    private Vector3 CalculateSurfaceNormal(Vector3 a, Vector3 b, Vector3 c)
    {
        var u = b - a;
        var v = c - a;
        var n = Vector3.Cross(u, v);
        return n.normalized;
    }

    #endregion
}
