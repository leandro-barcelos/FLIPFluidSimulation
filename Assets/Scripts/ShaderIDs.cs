using UnityEngine;

public class ShaderIDs : MonoBehaviour
{
    public static readonly int Properties = Shader.PropertyToID("_Properties");
    public static readonly int TimeStep = Shader.PropertyToID("_TimeStep");
    public static readonly int Gravity = Shader.PropertyToID("_Gravity");
    public static readonly int ParticleSize = Shader.PropertyToID("_ParticleSize");
    public static readonly int ParticlePos = Shader.PropertyToID("_ParticlePos");
    public static readonly int ParticleVel = Shader.PropertyToID("_ParticleVel");
    public static readonly int MaxX = Shader.PropertyToID("_MaxX");
    public static readonly int MinX = Shader.PropertyToID("_MinX");
    public static readonly int MaxY = Shader.PropertyToID("_MaxY");
    public static readonly int MinY = Shader.PropertyToID("_MinY");
    public static readonly int MaxZ = Shader.PropertyToID("_MaxZ");
    public static readonly int MinZ = Shader.PropertyToID("_MinZ");
    public static readonly int GridVel = Shader.PropertyToID("_GridVel");
    public static readonly int GridSize = Shader.PropertyToID("_GridSize");
    public static readonly int CellSpacing = Shader.PropertyToID("_CellSpacing");
    public static readonly int Offset = Shader.PropertyToID("_Offset");
}
