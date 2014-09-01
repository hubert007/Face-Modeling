#version 130 
precision mediump float;

in  vec4 vColor;
in  vec2 vTexCoord;
uniform sampler2D face_tex;

out vec4 frag_color;

void main() 
{ 
    
    frag_color = vColor;
    if(vColor.y > 0.1){
//        frag_color = texture(face_tex, vTexCoord);
    }
} 
