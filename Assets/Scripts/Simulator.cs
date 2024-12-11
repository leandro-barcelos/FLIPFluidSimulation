using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.Rendering;

public class Simulator
{
    public int particlesWidth, particlesHeight;

    public int gridWidth, gridHeight, gridDepth;

    public int gridResolutionX, gridResolutionY, gridResolutionZ;

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

    #region Initializations

    public Simulator(int particlesWidth, int particlesHeight, Vector3[] particlePositions, Vector3Int gridSize, Vector3Int gridResolution, int particleDensity, float flipness)
    {
        this.flipness = flipness;

        this.particlesWidth = particlesWidth;
        this.particlesHeight = particlesHeight;

        gridWidth = gridSize.x;
        gridHeight = gridSize.y;
        gridDepth = gridSize.z;

        gridResolutionX = gridResolution.x;
        gridResolutionY = gridResolution.y;
        gridResolutionZ = gridResolution.z;

        this.particleDensity = particleDensity;

        InitializeBuffers();
        InitializeParticleTextures(particlePositions);
        InitializeSimulationTextures();
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

        // TODO: transfer to grid
        // TODO: normalize
        // TODO: mark cells with fluid
        // TODO: save our original velocity grid
        // TODO: add forces to velocity grid
        // TODO: enforce boundaries
        // TODO: compute divergence
        // TODO: compute pressure via jacobi
        // TODO: subtract pressure from velocity
        // TODO: transfer velocities back to particles
        // TODO: advect particles
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
            filterMode = FilterMode.Point,
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
}
