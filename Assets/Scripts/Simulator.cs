using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.Rendering;

public class Simulator
{
    const int NumThreads = 8;

    public int particlesWidth, particlesHeight;
    private Vector2 invParticleResolution;

    public int gridWidth, gridHeight, gridDepth;

    public int gridResolutionX, gridResolutionY, gridResolutionZ;
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

    ComputeShader transferToGridShader, normalizeGridShader, addForcesShader, transferToParticlesShader;

    #region Initializations

    public Simulator(int particlesWidth, int particlesHeight, Vector3[] particlePositions, Vector3Int gridSize, Vector3Int gridResolution, int particleDensity, float flipness)
    {
        InitShaders();

        this.flipness = flipness;

        this.particlesWidth = particlesWidth;
        this.particlesHeight = particlesHeight;

        invParticleResolution = new(1f / particlesWidth, 1f / particlesHeight);

        gridWidth = gridSize.x;
        gridHeight = gridSize.y;
        gridDepth = gridSize.z;

        gridResolutionX = gridResolution.x;
        gridResolutionY = gridResolution.y;
        gridResolutionZ = gridResolution.z;

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

    public void Simulate(float timeStep, Vector3 mouseVelocity, Vector3 mouseRayOrigin, Vector3 mouseRayDirection)
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
        // TODO: advect particles and update mesh properties
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
        transferToGridShader.SetVector("_GridResolution", new(gridResolutionX, gridResolutionY, gridResolutionZ));
        transferToGridShader.SetVector("_GridSize", new Vector3(gridWidth, gridHeight, gridDepth));
        transferToGridShader.SetVector("_InvParticlesResolution", invParticleResolution);

        transferToGridShader.SetInt("_Accumulate", 0);

        var splatDepth = 5;
        for (int z = -(splatDepth - 1) / 2; z <= (splatDepth - 1) / 2; ++z)
        {
            transferToGridShader.SetInt("_ZOffset", z);

            // Set textures
            transferToGridShader.SetTexture(0, "_ParticlePositionTexture", particlePositionTexture);
            transferToGridShader.SetTexture(0, "_ParticleVelocityTexture", particleVelocityTexture);
            transferToGridShader.SetTexture(0, "_GridOutput", weightTexture);

            // Dispatch the compute shader
            int threadGroupsX = Mathf.CeilToInt(particlesWidth / NumThreads);
            int threadGroupsY = Mathf.CeilToInt(particlesHeight / NumThreads);
            transferToGridShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
        }

        transferToGridShader.SetInt("_Accumulate", 1);

        for (int z = -(splatDepth - 1) / 2; z <= (splatDepth - 1) / 2; ++z)
        {
            transferToGridShader.SetInt("_ZOffset", z);

            // Set textures
            transferToGridShader.SetTexture(0, "_ParticlePositionTexture", particlePositionTexture);
            transferToGridShader.SetTexture(0, "_ParticleVelocityTexture", particleVelocityTexture);
            transferToGridShader.SetTexture(0, "_GridOutput", tempVelocityTexture);

            // Dispatch the compute shader
            int threadGroupsX = Mathf.CeilToInt(particlesWidth / NumThreads);
            int threadGroupsY = Mathf.CeilToInt(particlesHeight / NumThreads);
            transferToGridShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
        }
    }

    private void NormalizeGrid()
    {
        // Set shader parameters
        normalizeGridShader.SetVector("_InvGridResolution", invGridResolution);

        // Set textures
        normalizeGridShader.SetTexture(0, "_TempVelocityTexture", tempVelocityTexture);
        normalizeGridShader.SetTexture(0, "_VelocityTexture", velocityTexture);
        normalizeGridShader.SetTexture(0, "_WeightTexture", weightTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt(gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt(gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt(gridResolutionZ / NumThreads);
        normalizeGridShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);
    }

    private void AddForces(float timeStep, Vector3 mouseVelocity, Vector3 mouseRayOrigin, Vector3 mouseRayDirection)
    {
        // Set shader parameters
        addForcesShader.SetVector("_InvGridResolution", invGridResolution);
        addForcesShader.SetVector("_GridSize", new(gridWidth, gridHeight, gridDepth));
        addForcesShader.SetVector("_MouseVelocity", mouseVelocity);
        addForcesShader.SetVector("_MouseRayOrigin", mouseRayOrigin);
        addForcesShader.SetVector("_MouseRayDirection", mouseRayDirection);
        addForcesShader.SetFloat("_TimeStep", timeStep);

        // Set textures
        addForcesShader.SetTexture(0, "_TempVelocityTexture", tempVelocityTexture);
        addForcesShader.SetTexture(0, "_VelocityTexture", velocityTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt(gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt(gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt(gridResolutionZ / NumThreads);
        addForcesShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);

        Swap(velocityTexture, tempVelocityTexture);
    }

    private void TransferToParticles()
    {
        // Set shader parameters
        transferToParticlesShader.SetVector("_InvParticleResolution", invParticleResolution);
        transferToParticlesShader.SetVector("_InvGridResolution", invGridResolution);
        transferToParticlesShader.SetVector("_GridResolution", new(gridResolutionX, gridResolutionY, gridResolutionZ));
        transferToParticlesShader.SetVector("_GridSize", new(gridWidth, gridHeight, gridDepth));
        transferToParticlesShader.SetFloat("_Flipness", flipness);

        // Set textures
        transferToParticlesShader.SetTexture(0, "_VelocityTexture", velocityTexture);
        transferToParticlesShader.SetTexture(0, "_OriginalVelocityTexture", originalVelocityTexture);
        transferToParticlesShader.SetTexture(0, "_ParticlePositionTexture", particlePositionTexture);
        transferToParticlesShader.SetTexture(0, "_ParticleVelocityTexture", particleVelocityTexture);
        transferToParticlesShader.SetTexture(0, "_ParticleVelocityTextureTemp", particleVelocityTextureTemp);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt(particlesWidth / NumThreads);
        int threadGroupsY = Mathf.CeilToInt(particlesHeight / NumThreads);
        transferToParticlesShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);

        Swap(particleVelocityTextureTemp, particleVelocityTexture);
    }
}
