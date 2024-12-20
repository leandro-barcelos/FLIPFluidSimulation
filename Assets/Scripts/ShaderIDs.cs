using UnityEngine;

public class ShaderIDs
{
    public static readonly int Properties = Shader.PropertyToID("_Properties");
    public static readonly int GridResolution = Shader.PropertyToID("_GridResolution");
    public static readonly int GridSize = Shader.PropertyToID("_GridSize");
    public static readonly int InvParticleResolution = Shader.PropertyToID("_InvParticleResolution");
    public static readonly int ParticleResolution = Shader.PropertyToID("_ParticleResolution");
    public static readonly int ParticlePositionTexture = Shader.PropertyToID("_ParticlePositionTexture");
    public static readonly int ParticleVelocityTexture = Shader.PropertyToID("_ParticleVelocityTexture");
    public static readonly int InvGridResolution = Shader.PropertyToID("_InvGridResolution");
    public static readonly int TempVelocityTexture = Shader.PropertyToID("_TempVelocityTexture");
    public static readonly int TempVelocityTextureX = Shader.PropertyToID("_TempVelocityTextureX");
    public static readonly int TempVelocityTextureY = Shader.PropertyToID("_TempVelocityTextureY");
    public static readonly int TempVelocityTextureZ = Shader.PropertyToID("_TempVelocityTextureZ");
    public static readonly int VelocityTexture = Shader.PropertyToID("_VelocityTexture");
    public static readonly int WeightTextureX = Shader.PropertyToID("_WeightTextureX");
    public static readonly int WeightTextureY = Shader.PropertyToID("_WeightTextureY");
    public static readonly int WeightTextureZ = Shader.PropertyToID("_WeightTextureZ");
    public static readonly int WeightTextureScalar = Shader.PropertyToID("_WeightTextureScalar");
    public static readonly int MouseVelocity = Shader.PropertyToID("_MouseVelocity");
    public static readonly int MouseRayOrigin = Shader.PropertyToID("_MouseRayOrigin");
    public static readonly int MouseRayDirection = Shader.PropertyToID("_MouseRayDirection");
    public static readonly int TimeStep = Shader.PropertyToID("_TimeStep");
    public static readonly int Flipness = Shader.PropertyToID("_Flipness");
    public static readonly int OriginalVelocityTexture = Shader.PropertyToID("_OriginalVelocityTexture");
    public static readonly int ParticleVelocityTextureTemp = Shader.PropertyToID("_ParticleVelocityTextureTemp");
    public static readonly int ParticleRadius = Shader.PropertyToID("_ParticleRadius");
    public static readonly int ParticlePositionTextureTemp = Shader.PropertyToID("_ParticlePositionTextureTemp");
    public static readonly int ParticleRandomTexture = Shader.PropertyToID("_ParticleRandomTexture");
    public static readonly int FrameNumber = Shader.PropertyToID("_FrameNumber");
    public static readonly int Offset = Shader.PropertyToID("_Offset");
    public static readonly int MarkerTexture = Shader.PropertyToID("_MarkerTexture");
    public static readonly int MaxDensity = Shader.PropertyToID("_MaxDensity");
    public static readonly int DivergenceTexture = Shader.PropertyToID("_DivergenceTexture");
    public static readonly int PressureTexture = Shader.PropertyToID("_PressureTexture");
    public static readonly int TempSimulationTexture = Shader.PropertyToID("_TempSimulationTexture");
}
