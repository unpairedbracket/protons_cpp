#version 450 core

in vec3 colour;

out vec3 color;

void main() {
    color = colour;
    vec2 circCoord = 2.0 * gl_PointCoord - 1.0;
    if (dot(circCoord, circCoord) > 1.0) {
        //discard;
    }
}
