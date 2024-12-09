using UnityEngine;

public class ObstacleMovement : MonoBehaviour
{
    private Vector3 direction = Vector3.zero;
    public float speed = 0.01f;
    public GameObject simulation;
    private FlipSimulation flipSimulation;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        flipSimulation = simulation.GetComponent<FlipSimulation>();
    }

    // Update is called once per frame
    void Update()
    {
        direction = Vector3.zero;

        if (Input.GetKey(KeyCode.A))
            direction += Vector3.left;
        if (Input.GetKey(KeyCode.D))
            direction += Vector3.right;

        if (Input.GetKey(KeyCode.W) && Input.GetKey(KeyCode.LeftShift))
            direction += Vector3.up;
        else if (Input.GetKey(KeyCode.W))
            direction += Vector3.forward;
        if (Input.GetKey(KeyCode.S) && Input.GetKey(KeyCode.LeftShift))
            direction += Vector3.down;
        else if (Input.GetKey(KeyCode.S))
            direction += Vector3.back;

        flipSimulation.SetObstacle(transform.position + direction * speed, false);
    }
}
