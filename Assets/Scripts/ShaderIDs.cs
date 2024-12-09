using UnityEngine;

public class ShaderIDs
{
    public static readonly int Properties = Shader.PropertyToID("_Properties");
    public static readonly int FDimensions = Shader.PropertyToID("_FDimensions");
    public static readonly int ParticlePos = Shader.PropertyToID("_ParticlePos");
    public static readonly int ParticleRadius = Shader.PropertyToID("_ParticleRadius");
    public static readonly int CellSpacing = Shader.PropertyToID("_CellSpacing");
    public static readonly int CellInvSpacing = Shader.PropertyToID("_CellInvSpacing");
    public static readonly int TimeStep = Shader.PropertyToID("_TimeStep");
    public static readonly int SimOrigin = Shader.PropertyToID("_SimOrigin");
    public static readonly int PDimensions = Shader.PropertyToID("_PDimensions");
    public static readonly int ParticleVel = Shader.PropertyToID("_ParticleVel");
    public static readonly int CellType = Shader.PropertyToID("_CellType");
    public static readonly int ObstaclesSDF = Shader.PropertyToID("_ObstaclesSDF");
    public static readonly int IntersectionCounts = Shader.PropertyToID("_IntersectionCounts");
    public static readonly int ClosestTriangles = Shader.PropertyToID("_ClosestTriangles");
    public static readonly int Transform = Shader.PropertyToID("_Transform");
    public static readonly int TriangleCount = Shader.PropertyToID("_TriangleCount");
    public static readonly int Vertices = Shader.PropertyToID("_Vertices");
    public static readonly int Triangles = Shader.PropertyToID("_Triangles");
    public static readonly int Direction = Shader.PropertyToID("_Direction");
    public static readonly int NumParticlesInCell = Shader.PropertyToID("_NumParticlesInCell");
    public static readonly int FirstCellParticle = Shader.PropertyToID("_FirstCellParticle");
    public static readonly int CellParticlesIds = Shader.PropertyToID("_CellParticlesIds");
    public static readonly int Gravity = Shader.PropertyToID("_Gravity");
    public static readonly int NumParticles = Shader.PropertyToID("_NumParticles");
    public static readonly int PInvSpacing = Shader.PropertyToID("_PInvSpacing");
    public static readonly int PNumCells = Shader.PropertyToID("_PNumCells");
    public static readonly int ParticleColor = Shader.PropertyToID("_ParticleColor");
    public static readonly int ColorDiffusionCoeff = Shader.PropertyToID("_ColorDiffusionCoeff");
    public static readonly int MinDist = Shader.PropertyToID("_MinDist");
    public static readonly int MinDist2 = Shader.PropertyToID("_MinDist2");
}
