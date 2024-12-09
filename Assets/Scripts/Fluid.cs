using System.IO.Compression;
using UnityEngine;
using UnityEngine.Profiling;

public class Fluid
{
    const int NumThreads = 8;

    public enum CellType
    {
        Fluid,
        Solid,
        Air
    }

    public Vector3 simOrigin { get; set; }

    // Fluid
    public float density { get; set; }
    public Vector3Int fDimensions { get; set; }
    public float h { get; set; }
    public float fInvSpacing { get; set; }
    public int fNumCells { get; set; }

    public ComputeBuffer uBuffer { get; set; }
    public ComputeBuffer vBuffer { get; set; }
    public ComputeBuffer wBuffer { get; set; }
    public float[] du { get; set; }
    public float[] dv { get; set; }
    public float[] dw { get; set; }
    public float[] prevU { get; set; }
    public float[] prevV { get; set; }
    public float[] prevW { get; set; }
    public float[] p { get; set; }
    public float[] s { get; set; }
    public ComputeBuffer sBuffer { get; set; }
    public ComputeBuffer cellTypeBuffer { get; set; }
    public Color[] cellColor { get; set; }

    // Particles
    public int maxParticles { get; set; }
    public ComputeBuffer particlePosBuffer { get; set; }
    public ComputeBuffer particleVelBuffer { get; set; }
    public ComputeBuffer particleColorBuffer { get; set; }
    public float[] particleDensity { get; set; }

    public float particleRestDensity { get; set; } = 0.0f;
    public float particleRadius { get; set; }
    public float pInvSpacing { get; set; }
    public Vector3Int pDimensions { get; set; }
    public int pNumCells { get; set; }

    public ComputeBuffer numCellParticles { get; set; }
    public ComputeBuffer firstCellParticle { get; set; }
    public ComputeBuffer cellParticleIds { get; set; }

    public int numParticles { get; set; } = 0;

    // Shaders
    private ComputeShader integrateParticlesShader;
    private ComputeShader pushParticlesApartShader;
    private ComputeShader handleParticleCollisionsShader;
    private ComputeShader updateParticlePropertiesShader;
    private ComputeShader updateCellTypeShader;
    private ComputeShader setObstacleShader;
    private void LoadShaders()
    {
        integrateParticlesShader = Resources.Load<ComputeShader>("IntegrateParticle");
        pushParticlesApartShader = Resources.Load<ComputeShader>("PushParticlesApart");
        handleParticleCollisionsShader = Resources.Load<ComputeShader>("HandleParticleCollisions");
        updateParticlePropertiesShader = Resources.Load<ComputeShader>("UpdateParticleProperties");
        updateCellTypeShader = Resources.Load<ComputeShader>("UpdateCellType");
        setObstacleShader = Resources.Load<ComputeShader>("SetObstacle");
    }

    public Fluid(Vector3 simOrigin, float density, Vector3 dimentions, float spacing, float particleRadius, int maxParticles)
    {
        LoadShaders();

        this.simOrigin = simOrigin;

        this.density = density;
        fDimensions = Vector3Int.FloorToInt(dimentions / spacing) + Vector3Int.one;
        h = Mathf.Max(dimentions.x / fDimensions.x, dimentions.y / fDimensions.y, dimentions.z / fDimensions.z);
        fInvSpacing = 1f / h;
        fNumCells = fDimensions.x * fDimensions.y * fDimensions.z;

        uBuffer = new ComputeBuffer(fNumCells, sizeof(float));
        vBuffer = new ComputeBuffer(fNumCells, sizeof(float));
        wBuffer = new ComputeBuffer(fNumCells, sizeof(float));
        du = new float[fNumCells];
        dv = new float[fNumCells];
        dw = new float[fNumCells];
        prevU = new float[fNumCells];
        prevV = new float[fNumCells];
        prevW = new float[fNumCells];
        p = new float[fNumCells];
        s = new float[fNumCells];
        cellColor = new Color[fNumCells];

        this.maxParticles = maxParticles;

        Color[] particleColor = new Color[maxParticles];
        for (var i = 0; i < maxParticles; i++)
            particleColor[i] = Color.blue;

        particleDensity = new float[fNumCells];

        this.particleRadius = particleRadius;
        pInvSpacing = 1f / (2.2f * particleRadius);
        pDimensions = Vector3Int.FloorToInt(dimentions * pInvSpacing) + Vector3Int.one;
        pNumCells = pDimensions.x * pDimensions.y * pDimensions.z;

        // Initialize ComputeBuffers
        particlePosBuffer = new ComputeBuffer(maxParticles, sizeof(float) * 3);
        particleVelBuffer = new ComputeBuffer(maxParticles, sizeof(float) * 3);
        particleColorBuffer = new ComputeBuffer(maxParticles, sizeof(float) * 4);
        particleColorBuffer.SetData(particleColor);
        cellTypeBuffer = new ComputeBuffer(fNumCells, sizeof(int));
        sBuffer = new ComputeBuffer(fNumCells, sizeof(float));
        sBuffer.SetData(s);
        numCellParticles = new ComputeBuffer(pNumCells, sizeof(int));
        firstCellParticle = new ComputeBuffer(pNumCells + 1, sizeof(int));
        cellParticleIds = new ComputeBuffer(maxParticles, sizeof(int));
    }

    public void Simulate(float dt, float gravity, float flipRatio, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift, bool separateParticles, Vector3 obstaclePos, Vector3 obstacleVel, float obstacleRadius, ComputeBuffer meshPropertiesBuffer)
    {
        IntegrateParticles(dt, gravity);
        if (separateParticles)
            PushParticlesApart();
        HandleParticleCollisions(obstaclePos, obstacleVel, obstacleRadius); // TODO: handle obstacle collisions
        UpdateCellType();
        UpdateMeshProperties(meshPropertiesBuffer);
    }

    private void IntegrateParticles(float dt, float gravity)
    {
        integrateParticlesShader.SetFloat(ShaderIDs.TimeStep, dt);
        integrateParticlesShader.SetFloat(ShaderIDs.Gravity, gravity);
        integrateParticlesShader.SetInt(ShaderIDs.NumParticles, numParticles);

        integrateParticlesShader.SetBuffer(0, ShaderIDs.ParticlePos, particlePosBuffer);
        integrateParticlesShader.SetBuffer(0, ShaderIDs.ParticleVel, particleVelBuffer);

        integrateParticlesShader.Dispatch(0, Mathf.CeilToInt((float)numParticles / NumThreads), 1, 1);
    }

    private void PushParticlesApart()
    {
        int[] initComputeBuffer = new int[pNumCells];
        numCellParticles.SetData(initComputeBuffer);

        initComputeBuffer = new int[maxParticles];
        cellParticleIds.SetData(initComputeBuffer);

        float minDist = 2f * particleRadius;

        pushParticlesApartShader.SetInt(ShaderIDs.NumParticles, numParticles);
        pushParticlesApartShader.SetInt(ShaderIDs.PNumCells, pNumCells);
        pushParticlesApartShader.SetFloat(ShaderIDs.PInvSpacing, pInvSpacing);
        pushParticlesApartShader.SetFloat(ShaderIDs.ColorDiffusionCoeff, 0.001f);
        pushParticlesApartShader.SetFloat(ShaderIDs.MinDist, minDist);
        pushParticlesApartShader.SetFloat(ShaderIDs.MinDist2, minDist * minDist);
        pushParticlesApartShader.SetVector(ShaderIDs.PDimensions, (Vector3)pDimensions);
        pushParticlesApartShader.SetVector(ShaderIDs.SimOrigin, simOrigin);

        CountParticlesInCell(numCellParticles);

        PartialSums(numCellParticles, firstCellParticle);

        FillParticlesIntoCells(cellParticleIds, firstCellParticle);

        Push(cellParticleIds, firstCellParticle);
    }

    private void Push(ComputeBuffer cellParticlesIds, ComputeBuffer firstCellParticle)
    {
        int mainKernel = pushParticlesApartShader.FindKernel("CSMain");
        pushParticlesApartShader.SetBuffer(mainKernel, ShaderIDs.ParticlePos, particlePosBuffer);
        pushParticlesApartShader.SetBuffer(mainKernel, ShaderIDs.CellParticlesIds, cellParticlesIds);
        pushParticlesApartShader.SetBuffer(mainKernel, ShaderIDs.FirstCellParticle, firstCellParticle);
        pushParticlesApartShader.SetBuffer(mainKernel, ShaderIDs.ParticleColor, particleColorBuffer);

        pushParticlesApartShader.Dispatch(mainKernel, Mathf.CeilToInt((float)numParticles / NumThreads), 1, 1);
    }

    private void FillParticlesIntoCells(ComputeBuffer cellParticlesIds, ComputeBuffer firstCellParticle)
    {
        int countKernel = pushParticlesApartShader.FindKernel("FillParticlesIntoCells");
        pushParticlesApartShader.SetBuffer(countKernel, ShaderIDs.ParticlePos, particlePosBuffer);
        pushParticlesApartShader.SetBuffer(countKernel, ShaderIDs.CellParticlesIds, cellParticlesIds);
        pushParticlesApartShader.SetBuffer(countKernel, ShaderIDs.FirstCellParticle, firstCellParticle);

        pushParticlesApartShader.Dispatch(countKernel, Mathf.CeilToInt((float)numParticles / NumThreads), 1, 1);
    }

    private void PartialSums(ComputeBuffer numParticlesInCell, ComputeBuffer firstCellParticle)
    {
        int partialSumsKernel = pushParticlesApartShader.FindKernel("PartialSums");
        pushParticlesApartShader.SetBuffer(partialSumsKernel, ShaderIDs.FirstCellParticle, firstCellParticle);
        pushParticlesApartShader.SetBuffer(partialSumsKernel, ShaderIDs.NumParticlesInCell, numParticlesInCell);

        pushParticlesApartShader.Dispatch(partialSumsKernel, Mathf.CeilToInt((float)pNumCells / NumThreads), 1, 1);
    }

    private void CountParticlesInCell(ComputeBuffer numParticlesInCell)
    {
        int countKernel = pushParticlesApartShader.FindKernel("CountParticlesInCell");
        pushParticlesApartShader.SetBuffer(countKernel, ShaderIDs.ParticlePos, particlePosBuffer);
        pushParticlesApartShader.SetBuffer(countKernel, ShaderIDs.NumParticlesInCell, numParticlesInCell);

        pushParticlesApartShader.Dispatch(countKernel, Mathf.CeilToInt((float)numParticles / NumThreads), 1, 1);
    }

    private void HandleParticleCollisions(Vector3 obstaclePos, Vector3 obstacleVel, float obstacleRadius)
    {
        handleParticleCollisionsShader.SetFloat(ShaderIDs.ParticleRadius, particleRadius);
        handleParticleCollisionsShader.SetFloat(ShaderIDs.ObstacleRadius, obstacleRadius);
        handleParticleCollisionsShader.SetFloat(ShaderIDs.CellSpacing, h);
        handleParticleCollisionsShader.SetFloat(ShaderIDs.NumParticles, numParticles);
        // handleParticleCollisionsShader.SetFloat(ShaderIDs.CellInvSpacing, fInvSpacing);
        // handleParticleCollisionsShader.SetFloat(ShaderIDs.TimeStep, dt);

        handleParticleCollisionsShader.SetVector(ShaderIDs.FDimensions, (Vector3)fDimensions);
        handleParticleCollisionsShader.SetVector(ShaderIDs.SimOrigin, simOrigin);
        handleParticleCollisionsShader.SetVector(ShaderIDs.ObstaclePos, obstaclePos);
        handleParticleCollisionsShader.SetVector(ShaderIDs.ObstacleVel, obstacleVel);

        handleParticleCollisionsShader.SetBuffer(0, ShaderIDs.ParticlePos, particlePosBuffer);
        handleParticleCollisionsShader.SetBuffer(0, ShaderIDs.ParticleVel, particleVelBuffer);
        // handleParticleCollisionsShader.SetBuffer(0, ShaderIDs.CellType, cellTypeBuffer);
        // handleParticleCollisionsShader.SetBuffer(0, ShaderIDs.ObstaclesSDF, _obstaclesSDF);

        handleParticleCollisionsShader.Dispatch(0, Mathf.CeilToInt((float)numParticles / NumThreads), 1, 1);

        // _obstaclesSDF.Release();
    }

    private void UpdateMeshProperties(ComputeBuffer meshPropertiesBuffer)
    {

        updateParticlePropertiesShader.SetFloat(ShaderIDs.NumParticles, numParticles);

        updateParticlePropertiesShader.SetBuffer(0, ShaderIDs.ParticlePos, particlePosBuffer);
        updateParticlePropertiesShader.SetBuffer(0, ShaderIDs.ParticleColor, particleColorBuffer);
        updateParticlePropertiesShader.SetBuffer(0, ShaderIDs.Properties, meshPropertiesBuffer);

        updateParticlePropertiesShader.Dispatch(0,
            Mathf.CeilToInt((float)numParticles / NumThreads), 1, 1);
    }

    private void UpdateCellType()
    {
        updateCellTypeShader.SetInt(ShaderIDs.FNumCells, fNumCells);
        updateCellTypeShader.SetInt(ShaderIDs.NumParticles, numParticles);
        updateCellTypeShader.SetFloat(ShaderIDs.FInvSpacing, fInvSpacing);
        updateCellTypeShader.SetVector(ShaderIDs.SimOrigin, simOrigin);
        updateCellTypeShader.SetVector(ShaderIDs.FDimensions, (Vector3)fDimensions);

        int sKernel = updateCellTypeShader.FindKernel("BasedOnS");
        updateCellTypeShader.SetBuffer(sKernel, ShaderIDs.CellType, cellTypeBuffer);
        updateCellTypeShader.SetBuffer(sKernel, ShaderIDs.S, sBuffer);

        updateCellTypeShader.Dispatch(sKernel, Mathf.CeilToInt((float)fNumCells / NumThreads), 1, 1);

        int particleKernel = updateCellTypeShader.FindKernel("BasedOnParticles");
        updateCellTypeShader.SetBuffer(particleKernel, ShaderIDs.CellType, cellTypeBuffer);
        updateCellTypeShader.SetBuffer(particleKernel, ShaderIDs.ParticlePos, particlePosBuffer);

        updateCellTypeShader.Dispatch(particleKernel, Mathf.CeilToInt((float)numParticles / NumThreads), 1, 1);
    }

    public void SetObstacle(Vector3 obstaclePos, Vector3 obstacleVel, float obstacleRadius)
    {
        setObstacleShader.SetFloat(ShaderIDs.H, h);
        setObstacleShader.SetFloat(ShaderIDs.ObstacleRadius, obstacleRadius);
        setObstacleShader.SetVector(ShaderIDs.ObstaclePos, obstaclePos);
        setObstacleShader.SetVector(ShaderIDs.ObstacleVel, obstacleVel);

        setObstacleShader.SetBuffer(0, ShaderIDs.S, sBuffer);
        setObstacleShader.SetBuffer(0, ShaderIDs.U, uBuffer);
        setObstacleShader.SetBuffer(0, ShaderIDs.V, vBuffer);
        setObstacleShader.SetBuffer(0, ShaderIDs.W, wBuffer);

        setObstacleShader.Dispatch(0, Mathf.CeilToInt((float)fDimensions.x / NumThreads), Mathf.CeilToInt((float)fDimensions.y / NumThreads), Mathf.CeilToInt((float)fDimensions.z / NumThreads));
    }

    public void Destroy()
    {
        particlePosBuffer.Release();
        particleVelBuffer.Release();
        particleColorBuffer.Release();
        cellTypeBuffer.Release();
        sBuffer.Release();
        numCellParticles.Release();
        firstCellParticle.Release();
        cellParticleIds.Release();
        uBuffer.Release();
        vBuffer.Release();
        wBuffer.Release();
    }
}
