using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.Rendering;

public class Simulator
{
    const int NumThreads = 8;

    public int particlesWidth, particlesHeight;
    private Vector2 particleResolution;
    private Vector2 invParticleResolution;

    public int gridWidth, gridHeight, gridDepth;
    private Vector3 gridSize;

    public int gridResolutionX, gridResolutionY, gridResolutionZ;
    private Vector3 gridResolution;
    private Vector3 invGridResolution;

    public int particleDensity = 0;

    // Parameters
    public float flipness = 0.99f;
    public int frameNumber = 0;

    // Objects
    public ComputeBuffer quadVertexBuffer, particleVertexBuffer;

    public RenderTexture particlePositionTexture, particlePositionTextureTemp;

    public RenderTexture particleVelocityTexture, particleVelocityTextureTemp;

    public RenderTexture particleRandomTexture;

    // Simulation Textures
    public RenderTexture velocityTexture, tempVelocityTexture, originalVelocityTexture, weightTexture;

    public RenderTexture markerTexture, divergenceTexture, pressureTexture, tempSimulationTexture;

    ComputeShader transferToGridShader, normalizeGridShader, addForcesShader, transferToParticlesShader, updateMeshPropertiesShader;

    #region Initializations

    public Simulator(int particlesWidth, int particlesHeight, Vector3[] particlePositions, Vector3Int gridSize, Vector3Int gridResolution, int particleDensity, float flipness)
    {
        InitShaders();

        this.flipness = flipness;

        this.particlesWidth = particlesWidth;
        this.particlesHeight = particlesHeight;

        particleResolution = new(particlesWidth, particlesHeight);
        invParticleResolution = new(1f / particlesWidth, 1f / particlesHeight);

        gridWidth = gridSize.x;
        gridHeight = gridSize.y;
        gridDepth = gridSize.z;

        this.gridSize = (Vector3)gridSize;

        gridResolutionX = gridResolution.x;
        gridResolutionY = gridResolution.y;
        gridResolutionZ = gridResolution.z;

        this.gridResolution = (Vector3)gridResolution;
        invGridResolution = new(1f / gridResolutionX, 1f / gridResolutionY, 1f / gridResolutionZ);

        this.particleDensity = particleDensity;

        InitializeBuffers();
        InitializeParticleTextures(particlePositions);
        InitializeSimulationTextures();
    }

    private void InitShaders()
    {
        transferToGridShader = Resources.Load<ComputeShader>("TransferToGrid");
        normalizeGridShader = Resources.Load<ComputeShader>("NormalizeGrid");
        addForcesShader = Resources.Load<ComputeShader>("AddForces");
        transferToParticlesShader = Resources.Load<ComputeShader>("TransferToParticles");
        updateMeshPropertiesShader = Resources.Load<ComputeShader>("UpdateMeshProperties");
    }

    private void InitializeSimulationTextures()
    {
        // Create simulation textures
        velocityTexture = CreateRenderTexture3D(gridResolutionX + 1, gridResolutionY + 1, gridResolutionZ + 1, RenderTextureFormat.ARGBHalf);

        tempVelocityTexture = CreateRenderTexture3D(gridResolutionX + 1, gridResolutionY + 1, gridResolutionZ + 1, RenderTextureFormat.ARGBHalf);

        originalVelocityTexture = CreateRenderTexture3D(gridResolutionX + 1, gridResolutionY + 1, gridResolutionZ + 1, RenderTextureFormat.ARGBHalf);

        weightTexture = CreateRenderTexture3D(gridResolutionX + 1, gridResolutionY + 1, gridResolutionZ + 1, RenderTextureFormat.ARGBHalf);

        markerTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.ARGBHalf);

        divergenceTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.ARGBHalf);

        pressureTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.ARGBHalf);

        tempSimulationTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.ARGBHalf);
    }

    private void InitializeParticleTextures(Vector3[] particlePositions)
    {
        // Generate initial particle positions
        Color[] particlePositionsData = new Color[particlesWidth * particlesHeight];
        Color[] particleRandoms = new Color[particlesWidth * particlesHeight];
        for (var i = 0; i < particlesWidth * particlesHeight; i++)
        {
            particlePositionsData[i] = new(
                particlePositions[i].x,
                particlePositions[i].y,
                particlePositions[i].z,
                0.0f
            );

            var theta = Random.value * 2f * Mathf.PI;
            var u = Random.value * 2f - 1f;
            particleRandoms[i] = new(
                Mathf.Sqrt(1f - u * u) * Mathf.Cos(theta),
                Mathf.Sqrt(1f - u * u) * Mathf.Sin(theta),
                u,
                0.0f
            );
        }

        Texture2D particlePositionsTexture2D = CreateTexture2D(particlesWidth, particlesHeight, particlePositionsData);
        particlePositionTexture = CreateRenderTexture2D(particlesWidth, particlesHeight, RenderTextureFormat.ARGBFloat);
        Graphics.Blit(particlePositionsTexture2D, particlePositionTexture);

        particlePositionTextureTemp = CreateRenderTexture2D(particlesWidth, particlesHeight, RenderTextureFormat.ARGBFloat);

        particleVelocityTexture = CreateRenderTexture2D(particlesWidth, particlesHeight, RenderTextureFormat.ARGBHalf);

        particleVelocityTextureTemp = CreateRenderTexture2D(particlesWidth, particlesHeight, RenderTextureFormat.ARGBHalf);

        Texture2D particleRandomsTexture2D = CreateTexture2D(particlesWidth, particlesHeight, particleRandoms);
        particleRandomTexture = CreateRenderTexture2D(particlesWidth, particlesHeight, RenderTextureFormat.ARGBFloat);
        Graphics.Blit(particleRandomsTexture2D, particleRandomTexture);
    }

    private void InitializeBuffers()
    {
        quadVertexBuffer = new ComputeBuffer(8, sizeof(float));
        float[] quadVertex = new float[] { -1.0f, -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, 1.0f, 1.0f };
        quadVertexBuffer.SetData(quadVertex);

        particleVertexBuffer = new ComputeBuffer(particlesHeight * particlesWidth, sizeof(float) * 2);
        Vector2[] particleTextureCoordinates = new Vector2[particlesWidth * particlesHeight];
        for (var y = 0; y < particlesHeight; y++)
        {
            for (var x = 0; x < particlesWidth; x++)
            {
                particleTextureCoordinates[y * particlesWidth + x] = new((x + 0.5f) / particlesWidth, (y + 0.5f) / particlesHeight);
            }
        }
        particleVertexBuffer.SetData(particleTextureCoordinates);
    }

    #endregion

    public void Simulate(float timeStep, Vector3 mouseVelocity, Vector3 mouseRayOrigin, Vector3 mouseRayDirection, ComputeBuffer _meshPropertiesBuffer)
    {
        frameNumber++;

        TransferToGrid();
        NormalizeGrid();
        // TODO: mark cells with fluid
        // TODO: save our original velocity grid
        AddForces(timeStep, mouseVelocity, mouseRayOrigin, mouseRayDirection);
        // TODO: enforce boundaries
        // TODO: compute divergence
        // TODO: compute pressure via jacobi
        // TODO: subtract pressure from velocity
        TransferToParticles();
        // TODO: advect particles
        UpdateMeshProperties(_meshPropertiesBuffer);
    }

    #region Helper Functions

    private void Swap(RenderTexture a, RenderTexture b)
    {
        (b, a) = (a, b);
    }

    private static RenderTexture CreateRenderTexture2D(int width, int height, RenderTextureFormat format)
    {
        var rt = new RenderTexture(width, height, 0, format)
        {
            enableRandomWrite = true,
            filterMode = FilterMode.Point,
            wrapMode = TextureWrapMode.Clamp
        };
        rt.Create();
        return rt;
    }

    private static RenderTexture CreateRenderTexture3D(int width, int height, int depth, RenderTextureFormat format)
    {
        var rt = new RenderTexture(width, height, 0, format)
        {
            dimension = TextureDimension.Tex3D,
            volumeDepth = depth,
            enableRandomWrite = true,
            filterMode = FilterMode.Bilinear,
            wrapMode = TextureWrapMode.Clamp
        };
        rt.Create();
        return rt;
    }

    private static Texture2D CreateTexture2D(int width, int height, Color[] colors)
    {
        var texture = new Texture2D(width, height, TextureFormat.RGBAFloat, false);
        texture.SetPixels(colors);
        texture.Apply();
        return texture;
    }

    #endregion

    ~Simulator()
    {
        quadVertexBuffer.Release();
        particleVertexBuffer.Release();

        particlePositionTexture.Release();
        particlePositionTextureTemp.Release();

        particleVelocityTexture.Release();
        particleVelocityTextureTemp.Release();

        particleRandomTexture.Release();

        velocityTexture.Release();
        tempVelocityTexture.Release();
        originalVelocityTexture.Release();
        weightTexture.Release();

        markerTexture.Release();
        divergenceTexture.Release();
        pressureTexture.Release();
        tempSimulationTexture.Release();
    }

    private void TransferToGrid()
    {
        // Set shader parameters
        transferToGridShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        transferToGridShader.SetVector(ShaderIDs.GridSize, new Vector3(gridWidth, gridHeight, gridDepth));
        transferToGridShader.SetVector(ShaderIDs.InvParticleResolution, invParticleResolution);
        transferToGridShader.SetVector(ShaderIDs.ParticleResolution, particleResolution);

        transferToGridShader.SetInt("_Accumulate", 0);

        var splatDepth = 5;
        for (int z = -(splatDepth - 1) / 2; z <= (splatDepth - 1) / 2; ++z)
        {
            transferToGridShader.SetInt("_ZOffset", z);

            // Set textures
            transferToGridShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);
            transferToGridShader.SetTexture(0, ShaderIDs.ParticleVelocityTexture, particleVelocityTexture);
            transferToGridShader.SetTexture(0, ShaderIDs.GridOutput, weightTexture);

            // Dispatch the compute shader
            int threadGroupsX = Mathf.CeilToInt((float)particlesWidth / NumThreads);
            int threadGroupsY = Mathf.CeilToInt((float)particlesHeight / NumThreads);
            transferToGridShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
        }

        transferToGridShader.SetInt("_Accumulate", 1);

        for (int z = -(splatDepth - 1) / 2; z <= (splatDepth - 1) / 2; ++z)
        {
            transferToGridShader.SetInt("_ZOffset", z);

            // Set textures
            transferToGridShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);
            transferToGridShader.SetTexture(0, ShaderIDs.ParticleVelocityTexture, particleVelocityTexture);
            transferToGridShader.SetTexture(0, ShaderIDs.GridOutput, tempVelocityTexture);

            // Dispatch the compute shader
            int threadGroupsX = Mathf.CeilToInt((float)particlesWidth / NumThreads);
            int threadGroupsY = Mathf.CeilToInt((float)particlesHeight / NumThreads);
            transferToGridShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
        }
    }

    private void NormalizeGrid()
    {
        // Set shader parameters
        normalizeGridShader.SetVector(ShaderIDs.InvGridResolution, invGridResolution);
        normalizeGridShader.SetVector(ShaderIDs.GridResolution, gridResolution);

        // Set textures
        normalizeGridShader.SetTexture(0, ShaderIDs.TempVelocityTexture, tempVelocityTexture);
        normalizeGridShader.SetTexture(0, ShaderIDs.VelocityTexture, velocityTexture);
        normalizeGridShader.SetTexture(0, ShaderIDs.WeightTexture, weightTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);
        normalizeGridShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);
    }

    private void AddForces(float timeStep, Vector3 mouseVelocity, Vector3 mouseRayOrigin, Vector3 mouseRayDirection)
    {
        // Set shader parameters
        addForcesShader.SetVector(ShaderIDs.InvGridResolution, invGridResolution);
        addForcesShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        addForcesShader.SetVector(ShaderIDs.GridSize, gridSize);
        addForcesShader.SetVector(ShaderIDs.MouseVelocity, mouseVelocity);
        addForcesShader.SetVector(ShaderIDs.MouseRayOrigin, mouseRayOrigin);
        addForcesShader.SetVector(ShaderIDs.MouseRayDirection, mouseRayDirection);
        addForcesShader.SetFloat(ShaderIDs.TimeStep, timeStep);

        // Set textures
        addForcesShader.SetTexture(0, ShaderIDs.TempVelocityTexture, tempVelocityTexture);
        addForcesShader.SetTexture(0, ShaderIDs.VelocityTexture, velocityTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);
        addForcesShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);

        Swap(velocityTexture, tempVelocityTexture);
    }

    private void TransferToParticles()
    {
        // Set shader parameters
        transferToParticlesShader.SetVector(ShaderIDs.InvParticleResolution, invParticleResolution);
        transferToParticlesShader.SetVector(ShaderIDs.ParticleResolution, particleResolution);
        transferToParticlesShader.SetVector(ShaderIDs.InvGridResolution, invGridResolution);
        transferToParticlesShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        transferToParticlesShader.SetVector(ShaderIDs.GridSize, gridSize);
        transferToParticlesShader.SetFloat(ShaderIDs.Flipness, flipness);

        // Set textures
        transferToParticlesShader.SetTexture(0, ShaderIDs.VelocityTexture, velocityTexture);
        transferToParticlesShader.SetTexture(0, ShaderIDs.OriginalVelocityTexture, originalVelocityTexture);
        transferToParticlesShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);
        transferToParticlesShader.SetTexture(0, ShaderIDs.ParticleVelocityTexture, particleVelocityTexture);
        transferToParticlesShader.SetTexture(0, ShaderIDs.ParticleVelocityTextureTemp, particleVelocityTextureTemp);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)particlesWidth / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)particlesHeight / NumThreads);
        transferToParticlesShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);

        Swap(particleVelocityTextureTemp, particleVelocityTexture);
    }

    private void UpdateMeshProperties(ComputeBuffer _meshPropertiesBuffer)
    {
        // Set shader parameters
        updateMeshPropertiesShader.SetVector(ShaderIDs.ParticleResolution, particleResolution);
        updateMeshPropertiesShader.SetVector(ShaderIDs.InvParticleResolution, invParticleResolution);
        updateMeshPropertiesShader.SetVector(ShaderIDs.GridSize, gridSize);
        updateMeshPropertiesShader.SetFloat(ShaderIDs.ParticleRadius, 7f / gridResolutionX);

        // Set textures
        updateMeshPropertiesShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);
        updateMeshPropertiesShader.SetBuffer(0, ShaderIDs.Properties, _meshPropertiesBuffer);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)particlesWidth / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)particlesHeight / NumThreads);
        updateMeshPropertiesShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
    }
}
