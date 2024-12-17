using System.Threading.Tasks;
using Unity.Mathematics;
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

    ComputeShader transferToGridShader, normalizeGridShader, addForcesShader, transferToParticlesShader, updateMeshPropertiesShader, advectShader, copyShader, markShader, enforceBoundariesShader, divergenceShader, jacobiShader;

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

    private void InitShaders()
    {
        transferToGridShader = Resources.Load<ComputeShader>("TransferToGrid");
        normalizeGridShader = Resources.Load<ComputeShader>("NormalizeGrid");
        addForcesShader = Resources.Load<ComputeShader>("AddForces");
        transferToParticlesShader = Resources.Load<ComputeShader>("TransferToParticles");
        updateMeshPropertiesShader = Resources.Load<ComputeShader>("UpdateMeshProperties");
        advectShader = Resources.Load<ComputeShader>("Advect");
        copyShader = Resources.Load<ComputeShader>("Copy");
        markShader = Resources.Load<ComputeShader>("Mark");
        enforceBoundariesShader = Resources.Load<ComputeShader>("EnforceBoundaries");
        divergenceShader = Resources.Load<ComputeShader>("Divergence");
        jacobiShader = Resources.Load<ComputeShader>("Jacobi");
    }

    private void InitializeSimulationTextures()
    {
        // Create simulation textures
        velocityTexture = CreateRenderTexture3D(gridResolutionX + 1, gridResolutionY + 1, gridResolutionZ + 1, RenderTextureFormat.ARGBHalf);

        tempVelocityTexture = CreateRenderTexture3D(gridResolutionX + 1, gridResolutionY + 1, gridResolutionZ + 1, RenderTextureFormat.ARGBHalf);

        originalVelocityTexture = CreateRenderTexture3D(gridResolutionX + 1, gridResolutionY + 1, gridResolutionZ + 1, RenderTextureFormat.ARGBHalf);

        weightTexture = CreateRenderTexture3D(gridResolutionX + 1, gridResolutionY + 1, gridResolutionZ + 1, RenderTextureFormat.ARGBHalf);

        markerTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RInt);

        divergenceTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RFloat);

        pressureTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RFloat);

        tempSimulationTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RFloat);
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

            // Generate a random normalized direction
            var randomDirection = UnityEngine.Random.insideUnitSphere.normalized;
            particleRandoms[i] = new Color(
                randomDirection.x,
                randomDirection.y,
                randomDirection.z,
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

    #region Helper Functions

    private void Swap(ref RenderTexture a, ref RenderTexture b)
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

    #region Simulation
    public void Simulate(float timeStep, Vector3 mouseVelocity, Vector3 mouseRayOrigin, Vector3 mouseRayDirection, ComputeBuffer _particleMeshPropertiesBuffer, ComputeBuffer _cellMeshPropertiesBuffer)
    {
        frameNumber++;

        TransferToGrid();
        NormalizeGrid();
        Mark();
        Copy();
        AddForces(timeStep, mouseVelocity, mouseRayOrigin, mouseRayDirection);
        EnforceBoundaries();
        Divergence();
        Jacobi();
        // TODO: subtract pressure from velocity
        TransferToParticles();
        Advect(timeStep);
        UpdateMeshProperties(_particleMeshPropertiesBuffer, _cellMeshPropertiesBuffer);
    }

    private void TransferToGrid()
    {
        int4[] init = new int4[gridResolutionX * gridResolutionY * gridResolutionZ];
        ComputeBuffer intWeightBuffer = new(gridResolutionX * gridResolutionY * gridResolutionZ, sizeof(int) * 4);
        intWeightBuffer.SetData(init);

        ComputeBuffer intTempVelocityBuffer = new(gridResolutionX * gridResolutionY * gridResolutionZ, sizeof(int) * 4);
        intTempVelocityBuffer.SetData(init);

        // Set shader parameters
        transferToGridShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        transferToGridShader.SetVector(ShaderIDs.GridSize, gridSize);
        transferToGridShader.SetVector(ShaderIDs.ParticleResolution, particleResolution);

        int threadGroupsX;
        int threadGroupsY;
        int threadGroupsZ;

        var splatDepth = 3;
        for (int x = -(splatDepth - 1) / 2; x <= (splatDepth - 1) / 2; ++x)
        {
            for (int y = -(splatDepth - 1) / 2; y <= (splatDepth - 1) / 2; ++y)
            {
                for (int z = -(splatDepth - 1) / 2; z <= (splatDepth - 1) / 2; ++z)
                {
                    transferToGridShader.SetVector(ShaderIDs.Offset, new(x, y, z));

                    // Set textures
                    transferToGridShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);
                    transferToGridShader.SetTexture(0, ShaderIDs.ParticleVelocityTexture, particleVelocityTexture);

                    transferToGridShader.SetBuffer(0, "_IntWeight", intWeightBuffer);
                    transferToGridShader.SetBuffer(0, "_IntTempVelocity", intTempVelocityBuffer);

                    // Dispatch the compute shader
                    threadGroupsX = Mathf.CeilToInt((float)particlesWidth / NumThreads);
                    threadGroupsY = Mathf.CeilToInt((float)particlesHeight / NumThreads);
                    transferToGridShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
                }
            }
        }

        transferToGridShader.SetTexture(1, ShaderIDs.WeightTexture, weightTexture);
        transferToGridShader.SetTexture(1, ShaderIDs.TempVelocityTexture, tempVelocityTexture);

        transferToGridShader.SetBuffer(1, "_IntWeight", intWeightBuffer);
        transferToGridShader.SetBuffer(1, "_IntTempVelocity", intTempVelocityBuffer);

        threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);
        transferToGridShader.Dispatch(1, threadGroupsX, threadGroupsY, threadGroupsZ);

        intWeightBuffer.Release();
        intTempVelocityBuffer.Release();
    }

    private void NormalizeGrid()
    {
        // Set shader parameters
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

        Swap(ref velocityTexture, ref tempVelocityTexture);
    }

    private void TransferToParticles()
    {
        // Set shader parameters
        transferToParticlesShader.SetVector(ShaderIDs.ParticleResolution, particleResolution);
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

        Swap(ref particleVelocityTextureTemp, ref particleVelocityTexture);
    }

    private void UpdateMeshProperties(ComputeBuffer _particleMeshPropertiesBuffer, ComputeBuffer _cellMeshPropertiesBuffer)
    {
        // Set shader parameters
        updateMeshPropertiesShader.SetVector(ShaderIDs.ParticleResolution, particleResolution);
        updateMeshPropertiesShader.SetVector(ShaderIDs.GridSize, gridSize);
        updateMeshPropertiesShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        updateMeshPropertiesShader.SetFloat(ShaderIDs.ParticleRadius, 7f / gridResolutionX);

        // Set textures
        updateMeshPropertiesShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);
        updateMeshPropertiesShader.SetBuffer(0, ShaderIDs.Properties, _particleMeshPropertiesBuffer);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)particlesWidth / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)particlesHeight / NumThreads);
        updateMeshPropertiesShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);

        // Set textures
        updateMeshPropertiesShader.SetTexture(1, ShaderIDs.MarkerTexture, markerTexture);
        updateMeshPropertiesShader.SetBuffer(1, ShaderIDs.Properties, _cellMeshPropertiesBuffer);

        // Dispatch the compute shader
        threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);
        updateMeshPropertiesShader.Dispatch(1, threadGroupsX, threadGroupsY, threadGroupsZ);
    }

    private void Advect(float timeStep)
    {
        // Set shader parameters
        advectShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        advectShader.SetVector(ShaderIDs.ParticleResolution, particleResolution);
        advectShader.SetVector(ShaderIDs.GridSize, gridSize);
        advectShader.SetInt(ShaderIDs.FrameNumber, frameNumber);
        advectShader.SetFloat(ShaderIDs.TimeStep, timeStep);

        // Set textures
        particlePositionTextureTemp.Release();
        particlePositionTextureTemp = CreateRenderTexture2D(particlesWidth, particlesHeight, RenderTextureFormat.ARGBFloat);
        advectShader.SetTexture(0, ShaderIDs.ParticlePositionTextureTemp, particlePositionTextureTemp);
        advectShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);
        advectShader.SetTexture(0, ShaderIDs.ParticleRandomTexture, particleRandomTexture);
        advectShader.SetTexture(0, ShaderIDs.VelocityTexture, velocityTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)particlesWidth / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)particlesHeight / NumThreads);
        advectShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);

        Swap(ref particlePositionTextureTemp, ref particlePositionTexture);
    }

    private void Copy()
    {
        // Set textures
        copyShader.SetTexture(0, ShaderIDs.VelocityTexture, velocityTexture);
        copyShader.SetTexture(0, ShaderIDs.OriginalVelocityTexture, originalVelocityTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);
        copyShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);
    }

    private void Mark()
    {
        // Set shader parameters
        markShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        markShader.SetVector(ShaderIDs.GridSize, gridSize);
        markShader.SetVector(ShaderIDs.ParticleResolution, particleResolution);

        // Set textures
        markerTexture.Release();
        markerTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RInt);
        markShader.SetTexture(0, ShaderIDs.MarkerTexture, markerTexture);
        markShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)particlesWidth / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)particlesHeight / NumThreads);
        markShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
    }

    private void EnforceBoundaries()
    {
        // Set shader parameters
        enforceBoundariesShader.SetVector(ShaderIDs.GridResolution, gridResolution);

        // Set textures
        enforceBoundariesShader.SetTexture(0, ShaderIDs.TempVelocityTexture, tempVelocityTexture);
        enforceBoundariesShader.SetTexture(0, ShaderIDs.VelocityTexture, velocityTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);
        enforceBoundariesShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);

        Swap(ref velocityTexture, ref tempVelocityTexture);
    }

    private void Divergence()
    {
        // Set shader parameters
        divergenceShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        divergenceShader.SetFloat(ShaderIDs.MaxDensity, particleDensity);

        // Set textures
        divergenceTexture.Release();
        divergenceTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RFloat);
        divergenceShader.SetTexture(0, ShaderIDs.DivergenceTexture, divergenceTexture);
        divergenceShader.SetTexture(0, ShaderIDs.VelocityTexture, velocityTexture);
        divergenceShader.SetTexture(0, ShaderIDs.MarkerTexture, markerTexture);
        divergenceShader.SetTexture(0, ShaderIDs.WeightTexture, weightTexture);

        // Dispatch the compute shader
        int threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);
        divergenceShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);
    }

    private void Jacobi()
    {
        // Set shader parameters
        jacobiShader.SetVector(ShaderIDs.GridResolution, gridResolution);
        jacobiShader.SetVector(ShaderIDs.InvGridResolution, invGridResolution);

        // Set textures
        pressureTexture.Release();
        pressureTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RFloat);
        jacobiShader.SetTexture(0, ShaderIDs.MarkerTexture, markerTexture);
        jacobiShader.SetTexture(0, ShaderIDs.DivergenceTexture, divergenceTexture);

        var pressureJacobiIterations = 50;
        for (var i = 0; i < pressureJacobiIterations; i++)
        {
            // Set textures
            jacobiShader.SetTexture(0, ShaderIDs.PressureTexture, pressureTexture);
            jacobiShader.SetTexture(0, ShaderIDs.TempSimulationTexture, tempSimulationTexture);

            // Dispatch the compute shader
            int threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
            int threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
            int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);
            jacobiShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);

            Swap(ref tempSimulationTexture, ref pressureTexture);
        }
    }

    #endregion
}
