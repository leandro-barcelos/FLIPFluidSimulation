using System;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Serialization;

public class FlipSimulation : MonoBehaviour
{
    private static readonly int Properties = Shader.PropertyToID("_Properties");
    const int NumThreads = 8;

    #region Auxiliary Structures

    private struct MeshProperties
    {
        // ReSharper disable once NotAccessedField.Local
        public Matrix4x4 Mat;
        // ReSharper disable once NotAccessedField.Local
        public Vector4 Color;

        public static int Size() {
            return
                sizeof(float) * 4 * 4 + // matrix;
                sizeof(float) * 4;      // color;
        }
    }

    #endregion

    #region Public Variables

    [Header("Simulation Parameters")]
    public float gravityAcceleration = -9.81f;

    [Header("Fluid")]
    public float density = 1000.0f;
    public float waterHeight = 0.8f;
    public float waterWidth = 0.6f;
    public float waterDepth = 0.6f;

    [Header("Grid")]
    public int gridResolution = 100;

    [Header("Particles")]

    [Header("Rendering")]
    public float occlusionRange;
    public Material particleMaterial;

    #endregion

    #region Private Variables

    //  Grid
    private float _cellSpacing;
    private Vector3Int _gridDimensions;
    private Vector3[][][] _gridVel;

    // Particles
    private int _maxParticles;
    private float _particleRadius;
    private Vector3Int _particleDimensions;
    private Vector3[] _particlePos, _particleVel;
    private float[] _particleDensity;
    private Color[] _particleColor;

    // Rendering
    private Mesh _particleMesh;
    private ComputeBuffer _meshPropertiesBuffer;
    private Bounds _bounds;
    private ComputeBuffer _argsBuffer;

    #endregion

    void Start()
    {
        InitializeRendering();
        InitializeGrid();
        InitializeParticles();
        InitializeBuffers();
    }

    void Update()
    {
        IntegrateParticles(Time.fixedDeltaTime);

        // TODO: push particles apart

        HandleParticleCollisions(); // TODO: handle obstacle collisions

        // TODO: update cell types
        // TODO: transfer velocities to grid
        // TODO: update particle density
        // TODO: solve incompressibility

        TransferVelocityToParticle();

        // TODO: update grid properties

        UpdateMeshProperties(); // TODO: Update colors as well

        Graphics.DrawMeshInstancedIndirect(_particleMesh, 0, particleMaterial, _bounds, _argsBuffer);
    }

    private void OnDisable() {
        _meshPropertiesBuffer.Release();
        _argsBuffer.Release();
    }

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

        _gridVel = new Vector3[_gridDimensions.x][][];

        for (int i = 0; i < _gridDimensions.x; i++)
        {
            _gridVel[i] = new Vector3[_gridDimensions.y][];
            for (int j = 0; j < _gridDimensions.y; j++)
            {
                _gridVel[i][j] = new Vector3[_gridDimensions.z];
            }
        }
    }

    private void InitializeParticles()
    {
        _particleRadius = 0.3f * _cellSpacing;

        float dx = 2.0f * _particleRadius;
        float dy = Mathf.Sqrt(3.0f) / 2.0f * dx;
        float dz = Mathf.Sqrt(2.0f / 3.0f) * dx;

        Vector3 scale = transform.localScale;
        _particleDimensions = new Vector3Int(
            Mathf.FloorToInt((waterWidth * scale.x - 2.0f * _cellSpacing - 2.0f * _particleRadius) / dx),
            Mathf.FloorToInt((waterHeight * scale.y - 2.0f * _cellSpacing - 2.0f * _particleRadius) / dy),
            Mathf.FloorToInt((waterDepth * scale.z - 2.0f * _cellSpacing - 2.0f * _particleRadius) / dz)
        );

        _maxParticles = _particleDimensions.x * _particleDimensions.y * _particleDimensions.z;
    }

    private void InitializeBuffers()
    {
        uint[] args = new uint[] { 0, 0, 0, 0, 0 };
        args[0] = _particleMesh.GetIndexCount(0);
        args[1] = (uint)_maxParticles;
        args[2] = _particleMesh.GetIndexStart(0);
        args[3] = _particleMesh.GetBaseVertex(0);

        _argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _argsBuffer.SetData(args);

        _particleVel = new Vector3[_maxParticles];
        _particlePos = new Vector3[_maxParticles];

        MeshProperties[] properties = new MeshProperties[_maxParticles];

        Vector3 particleScale = new Vector3(_particleRadius * 2, _particleRadius * 2, _particleRadius * 2);
        Quaternion rotation = Quaternion.identity;

        Vector3 offset = transform.localScale * -0.5f;

        int index = 0;
        for (int i = 0; i < _particleDimensions.x; i++)
        {
            for (int j = 0; j < _particleDimensions.y; j++)
            {
                for (int k = 0; k < _particleDimensions.z; k++)
                {
                    float x = _particleRadius + i * (2.0f * _particleRadius) + offset.x;
                    float y = _particleRadius + j * (Mathf.Sqrt(3.0f) / 2.0f * (2.0f * _particleRadius)) + offset.y;
                    float z = _particleRadius + k * (Mathf.Sqrt(2.0f / 3.0f) * (2.0f * _particleRadius)) + offset.z;

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

        _meshPropertiesBuffer = new ComputeBuffer(_maxParticles, MeshProperties.Size());
        _meshPropertiesBuffer.SetData(properties);

        particleMaterial.SetBuffer(Properties, _meshPropertiesBuffer);
    }

    private void IntegrateParticles(float dt)
    {
        throw new NotImplementedException();
    }

    private void UpdateMeshProperties()
    {
        throw new NotImplementedException();
    }

    private void HandleParticleCollisions()
    {
        throw new NotImplementedException();
    }

    private void TransferVelocityToParticle()
    {
        throw new NotImplementedException();
    }
}
