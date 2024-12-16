using System;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;

public class Parti : MonoBehaviour
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

    [Header("Parameters")]
    [Range(0f, 1f / 60f)] public float timeStep = 1f / 60f;
    public int gridWidth = 40;
    public int gridHeight = 20;
    public int gridDepth = 20;
    public int particlesPerCell = 10;
    [Range(0.2f, 3f)] public float gridCellDensity = 0.5f;
    [Range(0.5f, 0.99f)] public float flipness = 0.99f;

    [Header("Rendering")]
    public float occlusionRange;
    public Material particleMaterial, cellMaterial;
    public bool renderParticles, renderGrid;

    private Simulator simulator;

    // Mouse
    private Vector3 lastMousePlane = Vector3.zero;

    // Camera
    private CameraOrbit cameraOrbit;
    private new Camera camera;

    // Boxes
    private List<Bounds> boxes;

    // Rendering
    private Mesh _particleMesh, _cellMesh;
    private ComputeBuffer _particleMeshPropertiesBuffer, _particleArgsBuffer;
    private ComputeBuffer _cellMeshPropertiesBuffer, _cellArgsBuffer;
    private Bounds _bounds;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    private void Start()
    {
        transform.localScale = new(gridWidth, gridHeight, gridDepth);

        camera = Camera.main;
        cameraOrbit = camera.GetComponent<CameraOrbit>();
        InitializeBoxes();

        var desiredParticleCount = GetParticleCount();

        var particlesWidth = 512;
        var particlesHeight = Mathf.CeilToInt(desiredParticleCount / particlesWidth);

        var particleCount = particlesWidth * particlesHeight;

        var gridCells = gridWidth * gridHeight * gridDepth * gridCellDensity;

        var gridResolutionY = Mathf.CeilToInt(Mathf.Pow(gridCells / 2, 1.0f / 3.0f));
        var gridResolutionZ = gridResolutionY * 1;
        var gridResolutionX = gridResolutionY * 2;

        Vector3Int gridSize = new(gridWidth, gridHeight, gridDepth);
        Vector3Int gridResolution = new(gridResolutionX, gridResolutionY, gridResolutionZ);

        CreateCellPositions(gridSize, gridResolution);

        var sphereRadius = 7f / gridResolutionX;
        Vector3[] particlesPositions = CreateParticlePositions(particleCount, sphereRadius);

        simulator = new(particlesWidth, particlesHeight, particlesPositions, gridSize, gridResolution, particlesPerCell, flipness);
    }

    private Vector3[] CreateParticlePositions(int particleCount, float sphereRadius)
    {
        _particleMesh = OctahedronSphereCreator.Create(1, 1f);
        _bounds = new Bounds(transform.position, Vector3.one * (occlusionRange + 1));

        uint[] args = { 0, 0, 0, 0, 0 };
        args[0] = _particleMesh.GetIndexCount(0);
        args[1] = (uint)particleCount;
        args[2] = _particleMesh.GetIndexStart(0);
        args[3] = _particleMesh.GetBaseVertex(0);

        _particleArgsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _particleArgsBuffer.SetData(args);

        _particleMeshPropertiesBuffer = new ComputeBuffer(particleCount, MeshProperties.Size());

        Vector3[] particlesPositions = new Vector3[particleCount];
        MeshProperties[] properties = new MeshProperties[particleCount];
        var totalVolume = 0f;

        Vector3 particleScale = new(sphereRadius * 2, sphereRadius * 2, sphereRadius * 2);
        Quaternion rotation = Quaternion.identity;

        for (var i = 0; i < boxes.Count; i++)
        {
            var box = boxes[i];
            var volume = box.size.x * box.size.y * box.size.z;

            totalVolume += volume;
        }

        var particlesCreatedSoFar = 0;
        for (var i = 0; i < boxes.Count; i++)
        {
            var box = boxes[i];
            var volume = box.size.x * box.size.y * box.size.z;

            int particlesInBox;

            if (i < boxes.Count - 1)
            {
                particlesInBox = Mathf.FloorToInt(particleCount * volume / totalVolume);
            }
            else
            {
                particlesInBox = particleCount - particlesCreatedSoFar;
            }

            for (var j = 0; j < particlesInBox; j++)
            {
                var position = new Vector3(
                    Random.Range(box.min.x, box.max.x),
                    Random.Range(box.min.y, box.max.y),
                    Random.Range(box.min.z, box.max.z)
                );

                MeshProperties props = new()
                {
                    Mat = Matrix4x4.TRS(position, rotation, particleScale),
                    Color = Color.blue
                };

                particlesPositions[particlesCreatedSoFar + j] = position;
                properties[particlesCreatedSoFar + j] = props;
            }

            particlesCreatedSoFar += particlesInBox;
        }

        _particleMeshPropertiesBuffer.SetData(properties);
        particleMaterial.SetBuffer(ShaderIDs.Properties, _particleMeshPropertiesBuffer);

        return particlesPositions;
    }

    private void CreateCellPositions(Vector3Int gridSize, Vector3Int gridResolution)
    {
        int numCells = gridResolution.x * gridResolution.y * gridResolution.z;

        _cellMesh = Resources.GetBuiltinResource<Mesh>("Cube.fbx");
        _bounds = new Bounds(transform.position, Vector3.one * (occlusionRange + 1));

        uint[] args = { 0, 0, 0, 0, 0 };
        args[0] = _cellMesh.GetIndexCount(0);
        args[1] = (uint)numCells;
        args[2] = _cellMesh.GetIndexStart(0);
        args[3] = _cellMesh.GetBaseVertex(0);

        _cellArgsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _cellArgsBuffer.SetData(args);

        _cellMeshPropertiesBuffer = new ComputeBuffer(numCells, MeshProperties.Size());

        Vector3[] cellsPositions = new Vector3[numCells];
        MeshProperties[] properties = new MeshProperties[numCells];

        Vector3 cellScale = new(gridSize.x / (float)gridResolution.x, gridSize.y / (float)gridResolution.y, gridSize.z / (float)gridResolution.z);
        Quaternion rotation = Quaternion.identity;

        int index = 0;
        for (var i = 0; i < gridResolution.x; i++)
        {
            for (var j = 0; j < gridResolution.y; j++)
            {
                for (var k = 0; k < gridResolution.z; k++)
                {
                    cellsPositions[index] = new(i * cellScale.x + cellScale.x / 2, j * cellScale.y + cellScale.y / 2, k * cellScale.z + cellScale.z / 2);
                    properties[index] = new()
                    {
                        Mat = Matrix4x4.TRS(cellsPositions[index] - transform.localScale / 2, rotation, cellScale),
                        Color = new(0, 0, 0, 0)
                    };

                    index++;
                }
            }
        }

        _cellMeshPropertiesBuffer.SetData(properties);
        cellMaterial.SetBuffer(ShaderIDs.Properties, _cellMeshPropertiesBuffer);
    }

    private void InitializeBoxes()
    {
        boxes = new();
        foreach (var meshFilter in FindObjectsByType<MeshFilter>(FindObjectsSortMode.None))
        {
            if (meshFilter.gameObject.CompareTag("Spawn"))
            {
                var center = meshFilter.transform.position;
                var size = meshFilter.transform.localScale;
                var box = new Bounds(center, size);
                boxes.Add(box);
                meshFilter.gameObject.SetActive(false);
            }
        }
    }

    // Update is called once per frame
    private void Update()
    {
        UpdateMouse(out Vector3 worldSpaceMouseRay, out Vector3 worldMouseVelocity);

        simulator.Simulate(timeStep, worldMouseVelocity, camera.transform.position, worldSpaceMouseRay, _particleMeshPropertiesBuffer, _cellMeshPropertiesBuffer);

        if (renderParticles)
            Graphics.DrawMeshInstancedIndirect(_particleMesh, 0, particleMaterial, _bounds, _particleArgsBuffer);

        if (renderGrid)
            Graphics.DrawMeshInstancedIndirect(_cellMesh, 0, cellMaterial, _bounds, _cellArgsBuffer);
    }

    private void OnDestroy()
    {
        _particleMeshPropertiesBuffer?.Release();
        _particleArgsBuffer?.Release();
    }

    private void UpdateMouse(out Vector3 worldSpaceMouseRay, out Vector3 worldMouseVelocity)
    {
        worldSpaceMouseRay = Vector3.zero;
        worldMouseVelocity = Vector3.zero;

        if (Input.GetMouseButton(2))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(ray, out _))
            {
                Vector3 currentMousePlane = ray.GetPoint(cameraOrbit.distance);
                worldMouseVelocity = (currentMousePlane - lastMousePlane) / Time.deltaTime * 0.6f;
                lastMousePlane = currentMousePlane;

                worldSpaceMouseRay = ray.direction;
            }
        }
    }

    private int GetParticleCount()
    {
        var gridCells = gridWidth * gridHeight * gridDepth * gridCellDensity;

        var gridResolutionY = Mathf.CeilToInt(Mathf.Pow(gridCells / 2, 1.0f / 3.0f));
        var gridResolutionZ = gridResolutionY * 1;
        var gridResolutionX = gridResolutionY * 2;

        var totalGridCells = gridResolutionX * gridResolutionY * gridResolutionZ;

        var totalVolume = 0f;
        var cumulativeVolume = new float[boxes.Count];

        for (var i = 0; i < boxes.Count; i++)
        {
            var box = boxes[i];
            var volume = box.size.x * box.size.y * box.size.z;

            totalVolume += volume;
            cumulativeVolume[i] = totalVolume;
        }

        var fractionFilled = totalVolume / (gridWidth * gridHeight * gridDepth);

        var desiredParticleCount = fractionFilled * totalGridCells * particlesPerCell;

        return (int)desiredParticleCount;
    }


}
