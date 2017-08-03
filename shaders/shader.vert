#version 450 core

in vec3 pos;
in int running;

out vec3 colour;

uniform mat4 MVP;
 
void main() {
    if(running == 0) {
        colour = vec3(1, 0, 0);
    } else {
        colour = vec3(0, 1, 0);
    }
    gl_Position = MVP * vec4(pos.x, pos.y, pos.z, 1.0);
    gl_PointSize = 1;
}
