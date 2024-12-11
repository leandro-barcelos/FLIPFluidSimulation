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
    public Material particleMaterial;

    private Simulator simulator;

    // Mouse
    private Vector3 lastMousePlane = Vector3.zero;

    // Camera
    private CameraOrbit cameraOrbit;
    private new Camera camera;

    // Boxes
    private List<Bounds> boxes;

    // Rendering
    private Mesh _particleMesh;
    private ComputeBuffer _meshPropertiesBuffer, _argsBuffer;
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

        var sphereRadius = 7f / gridResolutionX;

        _particleMesh = Resources.GetBuiltinResource<Mesh>("Cube.fbx");
        _bounds = new Bounds(transform.position, Vector3.one * (occlusionRange + 1));

        uint[] args = { 0, 0, 0, 0, 0 };
        args[0] = _particleMesh.GetIndexCount(0);
        args[1] = (uint)particleCount;
        args[2] = _particleMesh.GetIndexStart(0);
        args[3] = _particleMesh.GetBaseVertex(0);

        _argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _argsBuffer.SetData(args);

        _meshPropertiesBuffer = new ComputeBuffer(particleCount, MeshProperties.Size());

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
                    Mat = Matrix4x4.TRS(position - transform.position, rotation, particleScale),
                    Color = Color.blue
                };

                particlesPositions[particlesCreatedSoFar + j] = position;
                properties[particlesCreatedSoFar + j] = props;
            }

            particlesCreatedSoFar += particlesInBox;
        }

        _meshPropertiesBuffer.SetData(properties);
        particleMaterial.SetBuffer(ShaderIDs.Properties, _meshPropertiesBuffer);

        simulator = new(particlesWidth, particlesHeight, particlesPositions, gridSize, gridResolution, particlesPerCell, flipness);
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

        simulator.Simulate(timeStep, worldMouseVelocity, camera.transform.position, worldSpaceMouseRay);

        Graphics.DrawMeshInstancedIndirect(_particleMesh, 0, particleMaterial, _bounds, _argsBuffer);
    }

    private void OnDestroy()
    {
        _meshPropertiesBuffer?.Release();
        _argsBuffer?.Release();
    }

    private void UpdateMouse(out Vector3 worldSpaceMouseRay, out Vector3 worldMouseVelocity)
    {
        // Calculate the field of view
        float fov = camera.fieldOfView * Mathf.Deg2Rad;

        // Get the mouse position in screen space
        Vector3 mousePosition = Input.mousePosition;

        // Convert the mouse position to normalized device coordinates (-1 to 1)
        Vector3 ndcMousePosition = new Vector3(
            (mousePosition.x / Screen.width) * 2 - 1,
            (mousePosition.y / Screen.height) * 2 - 1,
            -1.0f
        );

        // Calculate the view space mouse ray
        Vector3 viewSpaceMouseRay = new Vector3(
            ndcMousePosition.x * Mathf.Tan(fov / 2.0f) * (Screen.width / Screen.height),
            ndcMousePosition.y * Mathf.Tan(fov / 2.0f),
            -1.0f
        );

        // Calculate the mouse plane position
        Vector3 mousePlanePosition = viewSpaceMouseRay * cameraOrbit.distance;

        // Calculate the mouse velocity
        Vector3 mouseVelocity = mousePlanePosition - lastMousePlane;

        if (Input.GetMouseButton(0))
        {
            mouseVelocity = Vector3.zero;
        }

        lastMousePlane = mousePlanePosition;

        // Transform the view space mouse ray to world space
        worldSpaceMouseRay = camera.transform.TransformDirection(viewSpaceMouseRay);
        worldSpaceMouseRay.Normalize();

        // Get the camera's right and up vectors
        Vector3 cameraRight = camera.transform.right;
        Vector3 cameraUp = camera.transform.up;

        // Calculate the mouse velocity in world space
        worldMouseVelocity = mouseVelocity.x * cameraRight + mouseVelocity.y * cameraUp;
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
