#version 130

precision mediump float;

in vec4 in_position;
in vec4 in_color;
in vec4 in_normal;
in vec2 in_texcoord;

out vec4 vColor;
out vec2 vTexCoord;

uniform mat4 MVP;

void main()
{
 
    gl_Position = MVP * in_position;
    vColor = in_color;
    vTexCoord = in_texcoord;

}
