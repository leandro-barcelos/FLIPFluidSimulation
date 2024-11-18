using UnityEngine;

public class ShaderIDs : MonoBehaviour
{
    public static readonly int Properties = Shader.PropertyToID("_Properties");
    public static readonly int TimeStep = Shader.PropertyToID("_TimeStep");
    public static readonly int Gravity = Shader.PropertyToID("_Gravity");
    public static readonly int Size = Shader.PropertyToID("_Size");
    public static readonly int ParticlePos = Shader.PropertyToID("_ParticlePos");
    public static readonly int ParticleVel = Shader.PropertyToID("_ParticleVel");
}
