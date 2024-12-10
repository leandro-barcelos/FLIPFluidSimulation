using System;
using UnityEngine;
using UnityEngine.Profiling;

public class FlipSimulation : MonoBehaviour
{
    #region Auxiliary Structures
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
    public int gridResolution = 25;
    public float gravity = -9.81f;
    public float dt = 1.0f / 60.0f;
    [Range(0f, 1f)] public float flipRatio = 0.9f;
    public int numPressureIters = 50;
    public int numParticleIters = 2;
    public float overRelaxation = 1.9f;
    public bool compensateDrift = true;
    public bool separateParticles = true;
    public bool showParticles = true;
    public bool showGrid = false;
    private Fluid fluid;

    [Header("Rendering")]
    public float occlusionRange;
    public Material particleMaterial;

    [Header("Obstacles")]
    public GameObject obstacle;
    private Vector3 obstacleVel;
    // public ComputeShader generateObstaclesSDFShader;
    // public int FastSweepingPasses = 3;

    #endregion

    #region Private Variables

    // Rendering
    private Mesh _particleMesh;
    private ComputeBuffer _meshPropertiesBuffer, _argsBuffer;
    private Bounds _bounds;

    // Obstacles
    // private ComputeBuffer _obstaclesSDF;

    #endregion

    void Start()
    {
        var tankDimensions = 1.0f * transform.localScale;

        var cellSpacing = tankDimensions.y / gridResolution;
        var density = 1000.0f;

        var waterDimensions = new Vector3(0.8f, 0.6f, 0.6f);

        var particleRadius = 0.3f * cellSpacing;
        var dx = 2.0f * particleRadius;
        var dy = Mathf.Sqrt(3.0f) / 2.0f * dx;
        var dz = Mathf.Sqrt(3.0f) / 2.0f * dx;

        var numX = Mathf.FloorToInt((waterDimensions.x * tankDimensions.x - 2.0f * cellSpacing - 2.0f * particleRadius) / dx);
        var numY = Mathf.FloorToInt((waterDimensions.y * tankDimensions.y - 2.0f * cellSpacing - 2.0f * particleRadius) / dy);
        var numZ = Mathf.FloorToInt((waterDimensions.z * tankDimensions.z - 2.0f * cellSpacing - 2.0f * particleRadius) / dz);
        var maxParticles = numX * numY * numZ;

        var simOrigin = transform.position - transform.localScale / 2;

        fluid = new Fluid(simOrigin, density, tankDimensions, cellSpacing, particleRadius, maxParticles);

        // Rendering
        _particleMesh = Resources.GetBuiltinResource<Mesh>("Cube.fbx");
        _bounds = new Bounds(transform.position, Vector3.one * (occlusionRange + 1));

        uint[] args = { 0, 0, 0, 0, 0 };
        args[0] = _particleMesh.GetIndexCount(0);
        args[1] = (uint)maxParticles;
        args[2] = _particleMesh.GetIndexStart(0);
        args[3] = _particleMesh.GetBaseVertex(0);

        _argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _argsBuffer.SetData(args);

        _meshPropertiesBuffer = new ComputeBuffer(maxParticles, MeshProperties.Size());

        // Create particles
        var properties = new MeshProperties[maxParticles];
        var particlePos = new Vector3[maxParticles];

        Vector3 particleScale = new(particleRadius * 2, particleRadius * 2, particleRadius * 2);
        Quaternion rotation = Quaternion.identity;

        fluid.numParticles = numX * numY * numZ;
        var p = 0;
        for (int i = 0; i < numX; i++)
        {
            for (int j = 0; j < numY; j++)
            {
                for (int k = 0; k < numZ; k++)
                {
                    particlePos[p] = new Vector3(
                        cellSpacing + particleRadius + dx * i + (j % 2 == 0 ? 0.0f : particleRadius),
                        cellSpacing + particleRadius + dy * j,
                        cellSpacing + particleRadius + dz * k
                    ) + fluid.simOrigin;

                    MeshProperties props = new()
                    {
                        Mat = Matrix4x4.TRS(particlePos[p], rotation, particleScale),
                        Color = Color.white
                    };

                    properties[p++] = props;
                }
            }
        }

        fluid.particlePosBuffer.SetData(particlePos);
        _meshPropertiesBuffer.SetData(properties);
        particleMaterial.SetBuffer(ShaderIDs.Properties, _meshPropertiesBuffer);

        // Setup grid cells for tank
        var n = fluid.fDimensions.y;
        var m = fluid.fDimensions.z;
        for (int i = 0; i < fluid.fDimensions.x; i++)
        {
            for (int j = 0; j < fluid.fDimensions.y; j++)
            {
                for (int k = 0; k < fluid.fDimensions.z; k++)
                {
                    float s = 1.0f;
                    if (i == 0 || j == 0 || k == 0 || i == fluid.fDimensions.x - 1 || k == fluid.fDimensions.z - 1) s = 0.0f;

                    fluid.s[i * n * m + j * m + k] = s;
                }
            }
        }

        SetObstacle(new(0.25f, 0.0f, 0.25f), true);
    }

    void Update()
    {
        // _obstaclesSDF = new ComputeBuffer(fluid.fNumCells, sizeof(float));

        Vector3 obstacleScale = obstacle.transform.localScale;
        float obstacleRadius = Mathf.Max(obstacleScale.x, obstacleScale.y, obstacleScale.z) / 2;

        fluid.Simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters, overRelaxation, compensateDrift, separateParticles, obstacle.transform.position, obstacleVel, obstacleRadius, _meshPropertiesBuffer);

        // foreach (MeshFilter meshFilter in FindObjectsByType<MeshFilter>(FindObjectsSortMode.None))
        // {
        //     if (meshFilter.gameObject.name == "Simulation") continue;
        //     GenerateObstaclesSDF(meshFilter);
        // }

        if (showParticles)
            Graphics.DrawMeshInstancedIndirect(_particleMesh, 0, particleMaterial, _bounds, _argsBuffer);
    }

    private void OnDisable()
    {
        fluid.Destroy();
        _meshPropertiesBuffer?.Release();
        _argsBuffer?.Release();
        // _obstaclesSDF?.Release();
    }

    #region Simulation

    public void SetObstacle(Vector3 position, bool reset)
    {
        Vector3 vel = Vector3.zero;

        if (!reset)
            vel = (position - obstacle.transform.position) / dt;

        obstacle.transform.position = position;
        Vector3 scale = obstacle.transform.localScale;
        float r = Mathf.Max(scale.x, scale.y, scale.z) / 2;

        fluid.SetObstacle(position, vel, r);

        obstacleVel = vel;
    }

    // private void GenerateObstaclesSDF(MeshFilter meshFilter)
    // {
    //     _obstaclesSDF = new ComputeBuffer(fluid.fNumCells, sizeof(float));
    //     ComputeBuffer closestTriangles = new(fluid.fNumCells, sizeof(int));
    //     ComputeBuffer intersectionCounts = new(fluid.fNumCells, sizeof(int));
    //     int triangleCount = meshFilter.mesh.triangles.Length;
    //     ComputeBuffer triangles = new(triangleCount, sizeof(int));
    //     triangles.SetData(meshFilter.mesh.triangles);
    //     ComputeBuffer vertices = new(meshFilter.mesh.vertexCount, sizeof(float) * 3);
    //     vertices.SetData(meshFilter.mesh.vertices);

    //     InitSDFNearGeometry(meshFilter.transform.localToWorldMatrix, closestTriangles, intersectionCounts, triangleCount, triangles, vertices);

    //     FastSweeping(closestTriangles, vertices, triangles);

    //     CalculateSDFSign(intersectionCounts);

    //     closestTriangles.Release();
    //     intersectionCounts.Release();
    //     triangles.Release();
    //     vertices.Release();
    // }

    // private void CalculateSDFSign(ComputeBuffer intersectionCounts)
    // {
    //     int signKernel = generateObstaclesSDFShader.FindKernel("CalculateSDFSign");

    //     generateObstaclesSDFShader.SetVector(ShaderIDs.PDimensions, (Vector3)fluid.fDimensions);

    //     generateObstaclesSDFShader.SetBuffer(signKernel, ShaderIDs.ObstaclesSDF, _obstaclesSDF);
    //     generateObstaclesSDFShader.SetBuffer(signKernel, ShaderIDs.CellType, fluid.cellTypeBuffer);
    //     generateObstaclesSDFShader.SetBuffer(signKernel, ShaderIDs.IntersectionCounts, intersectionCounts);

    //     generateObstaclesSDFShader.Dispatch(signKernel, 1, Mathf.CeilToInt((float)fluid.fDimensions.y / NumThreads), Mathf.CeilToInt((float)fluid.fDimensions.z / NumThreads));
    // }

    // private void InitSDFNearGeometry(Matrix4x4 transform, ComputeBuffer closestTriangles, ComputeBuffer intersectionCounts, int triangleCount, ComputeBuffer triangles, ComputeBuffer vertices)
    // {
    //     int initKernel = generateObstaclesSDFShader.FindKernel("Init");
    //     generateObstaclesSDFShader.SetMatrix(ShaderIDs.Transform, transform);
    //     generateObstaclesSDFShader.SetVector(ShaderIDs.SimOrigin, fluid.simOrigin);
    //     generateObstaclesSDFShader.SetVector(ShaderIDs.PDimensions, (Vector3)fluid.fDimensions);
    //     generateObstaclesSDFShader.SetFloat(ShaderIDs.CellSpacing, fluid.h);
    //     generateObstaclesSDFShader.SetInt(ShaderIDs.TriangleCount, triangleCount);

    //     generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.CellType, fluid.cellTypeBuffer);
    //     generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.ObstaclesSDF, _obstaclesSDF);
    //     generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.ClosestTriangles, closestTriangles);
    //     generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.IntersectionCounts, intersectionCounts);
    //     generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.Vertices, vertices);
    //     generateObstaclesSDFShader.SetBuffer(initKernel, ShaderIDs.Triangles, triangles);

    //     generateObstaclesSDFShader.Dispatch(initKernel, Mathf.CeilToInt((float)fluid.fDimensions.x / NumThreads), Mathf.CeilToInt((float)fluid.fDimensions.y / NumThreads), Mathf.CeilToInt((float)fluid.fDimensions.z / NumThreads));
    // }

    // private void FastSweeping(ComputeBuffer closestTriangles, ComputeBuffer vertices, ComputeBuffer triangles)
    // {
    //     int fastSweepKernel = generateObstaclesSDFShader.FindKernel("FastSweeping");

    //     // Set parameters
    //     generateObstaclesSDFShader.SetVector(ShaderIDs.SimOrigin, fluid.simOrigin);
    //     generateObstaclesSDFShader.SetVector(ShaderIDs.PDimensions, (Vector3)fluid.fDimensions);
    //     generateObstaclesSDFShader.SetFloat(ShaderIDs.CellSpacing, fluid.h);

    //     // Set buffers
    //     generateObstaclesSDFShader.SetBuffer(fastSweepKernel, ShaderIDs.ObstaclesSDF, _obstaclesSDF);
    //     generateObstaclesSDFShader.SetBuffer(fastSweepKernel, ShaderIDs.ClosestTriangles, closestTriangles);
    //     generateObstaclesSDFShader.SetBuffer(fastSweepKernel, ShaderIDs.Vertices, vertices);
    //     generateObstaclesSDFShader.SetBuffer(fastSweepKernel, ShaderIDs.Triangles, triangles);

    //     int[] iDirections = { 1, -1, 1, -1, 1, -1, 1, -1 };
    //     int[] jDirections = { 1, 1, 1, 1, -1, -1, -1, -1 };
    //     int[] kDirections = { 1, -1, -1, 1, 1, -1, -1, 1 };

    //     // Dispatch in 8 directions
    //     for (int _ = 0; _ < FastSweepingPasses; _++)
    //     {
    //         for (int dir = 0; dir < 8; dir++)
    //         {
    //             Vector3 direction = new(iDirections[dir], jDirections[dir], kDirections[dir]);

    //             // Set directional parameters
    //             generateObstaclesSDFShader.SetVector(ShaderIDs.Direction, direction);

    //             // Dispatch compute shader for current direction
    //             generateObstaclesSDFShader.Dispatch(fastSweepKernel, Mathf.CeilToInt((float)fluid.fDimensions.x / NumThreads), Mathf.CeilToInt((float)fluid.fDimensions.y / NumThreads), Mathf.CeilToInt((float)fluid.fDimensions.z / NumThreads));
    //         }
    //     }
    // }
    #endregion
}
