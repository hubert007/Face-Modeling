// Homework #2 for Computer Graphic, UFL
// 
// Authoer: Ruijin Wu <ruijin@cise.ufl.edu>


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#ifdef __APPLE__
#include <GLEW/glew.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/freeglut.h>
#endif

#include "utils.h"
#include "debuglib.h"
#include "tga.h"


#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/transform.hpp"



// CONTSTANTS
const float PI = 3.1415926;
const int window_width = 500, window_height = 500;
const char * window_title = "HW2";

void draw_grid();

struct Vertex
{
    glm::vec4 pos;
    glm::vec4 color;
    glm::vec4 normal;

///////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    glm::vec2 tex_coord;
};


// runtime variable
char buf[1024];
vector<Vertex> grid_vertices;


GLuint vertex_shader, fragment_shader, program;



// function Prototype 
void graphics_init();

void initialize_points(){
    // create the grid

    for(int i = 0; i < 11 ; i++)
    {
        Vertex start,end;
        start.color = glm::vec4(1,1,1,1);
        end.color = glm::vec4(1,1,1,1);

        start.pos = glm::vec4(i - 5.0, 0, -5.0, 1.0); 
        end.pos = glm::vec4(i - 5.0, 0, 5.0, 1.0);
        grid_vertices.push_back(start);
        grid_vertices.push_back(end);
            
        start.pos = glm::vec4(-5.0, 0, i  - 5.0, 1.0); 
        end.pos = glm::vec4(5.0f, 0, i - 5.0, 1.0); 
        grid_vertices.push_back(start);
        grid_vertices.push_back(end);
    }
    // create the frame at origin

    Vertex vo, vx,vy,vz;
    vo.pos = glm::vec4(0,0.01,0,1);
    
    vx.pos = glm::vec4(3,0.01,0,1);
    vx.color = glm::vec4(1,0,0,1); // red
    
    vy.pos = glm::vec4(0,3.01,0,1);
    vy.color = glm::vec4(0,1,0,1); // green
    
    vz.pos = glm::vec4(0,0.01,3,1);
    vz.color = glm::vec4(0,0,1,1); // green

    vo.color = vx.color;
    grid_vertices.push_back(vo);
    grid_vertices.push_back(vx);

    vo.color = vy.color;
    grid_vertices.push_back(vo);
    grid_vertices.push_back(vy);
    
    vo.color = vz.color;
    grid_vertices.push_back(vo);
    grid_vertices.push_back(vz); 
}

Vertex control[12*12];
void initialize_controlpoints()
{
float bb = 2.75f;
float cc = 0.05f;
float dd = 0.05f;
for(int i=0;i<12;i++)
{
  float aa=-2.75f;
  float cc=0.05f;
  for(int j=0;j<12;j++)
  {
  control[j+i*12].pos=glm::vec4(aa,bb,0,1);
  control[j+i*12].color=glm::vec4(cc,dd,1,1);
  aa=aa+0.5f;
  cc=cc+0.05f;
  }
  bb=bb-0.25f;
  dd=dd+0.05f;
}

////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  float ff=0.0f;
for(int i=0;i<12;i++)
{
  float ee=0.0f;
  for(int j=0;j<12;j++)
  {
  control[j+i*12].tex_coord=glm::vec2(ee,ff);
  ee=ee+0.08f;
  }
  ff=ff+0.08f;
}
}



/*
void initialize_subdivision1(Vertex c1[], Vertex c2[], int np){
 //Pnew2i= (Poldi-1+6Poldi+Poldi+1)/8;
 //Pnew2i+1=(4Poldi+4Poldi+1)/8=(Poldi+Poldi+1)/2;
 //i=0
//np=24 for the first time
for(int i=0;i<np/2;i++)
{
  c2[np*i].pos[0]=(c1[np*i/2].pos[0]+6*c1[np*i/2].pos[0]+c1[np*i/2+1].pos[0])/8;
  c2[np*i].pos[1]=(c1[np*i/2].pos[1]+6*c1[np*i/2].pos[1]+c1[np*i/2+1].pos[1])/8;
  c2[np*i].pos[2]=(c1[np*i/2].pos[2]+6*c1[np*i/2].pos[2]+c1[np*i/2+1].pos[2])/8;

  c2[(i+1)*np-1].pos[0]=(c1[np/2*(i+1)-1].pos[0]+c1[np/2*(i+1)-1].pos[0])/2;
  c2[(i+1)*np-1].pos[1]=(c1[np/2*(i+1)-1].pos[1]+c1[np/2*(i+1)-1].pos[1])/2;
  c2[(i+1)*np-1].pos[2]=(c1[np/2*(i+1)-1].pos[2]+c1[np/2*(i+1)-1].pos[2])/2;


}

  for(int i=0;i<np/2;i++)
{
  for(int j=1; j<np-1;j=j+2)
 {
  c2[j+i*np].pos[0]=(c1[(j+np*i)/2-1].pos[0]+6*c1[(j+np*i)/2].pos[0]+c1[(j+np*i)/2+1].pos[0])/8;
  c2[j+i*np].pos[1]=(c1[(j+np*i)/2-1].pos[1]+6*c1[(j+np*i)/2].pos[1]+c1[(j+np*i)/2+1].pos[1])/8;
  c2[j+i*np].pos[2]=(c1[(j+np*i)/2-1].pos[2]+6*c1[(j+np*i)/2].pos[2]+c1[(j+np*i)/2+1].pos[2])/8;
 //cout<<j+i*np;
 //cout<<" ";
 }
  for(int j=2; j<np-1;j=j+2)
 {
  c2[j+i*np].pos[0]=(c1[(j+i*np-1)/2].pos[0]+c1[(j+i*np-1)/2+1].pos[0])/2;
  c2[j+i*np].pos[1]=(c1[(j+i*np-1)/2].pos[1]+c1[(j+i*np-1)/2+1].pos[1])/2;
  c2[j+i*np].pos[2]=(c1[(j+i*np-1)/2].pos[2]+c1[(j+i*np-1)/2+1].pos[2])/2;
 // cout<<j+i*np;
 // cout<<" ";
 }
}
  //i=np-1 e.g.143
for(int i=0; i<np*np/2; i++)
{
 c2[i].pos[3]=1.0f; 
 c2[i].color=glm::vec4(0.5f,0.5f,0.5f,1.0f);
}

}

void initialize_subdivision2(Vertex c1[], Vertex c2[], int np)
{
for(int i=0;i<np;i++)
{
  c2[i].pos[0]=(c1[i].pos[0]+6*c1[i].pos[0]+c1[np+i].pos[0])/8;
  c2[i].pos[1]=(c1[i].pos[1]+6*c1[i].pos[1]+c1[np+i].pos[1])/8;
  c2[i].pos[2]=(c1[i].pos[2]+6*c1[i].pos[2]+c1[np+i].pos[2])/8;

  c2[np*(np-1)+i].pos[0]=(c1[(np/2-1)*np+i].pos[0]+c1[(np/2-1)*np+i].pos[0])/2;
  c2[np*(np-1)+i].pos[1]=(c1[(np/2-1)*np+i].pos[1]+c1[(np/2-1)*np+i].pos[1])/2;
  c2[np*(np-1)+i].pos[2]=(c1[(np/2-1)*np+i].pos[2]+c1[(np/2-1)*np+i].pos[2])/2;
}

  for(int i=1;i<np-1;i=i+2)
{
  for(int j=0; j<np;j=j++)
 {
  c2[j+i*np].pos[0]=(c1[j+np*((i+1)/2-1)].pos[0]+6*c1[j+np*(i+1)/2].pos[0]+c1[j+np*((i+1)/2+1)].pos[0])/8;
  c2[j+i*np].pos[1]=(c1[j+np*((i+1)/2-1)].pos[1]+6*c1[j+np*(i+1)/2].pos[1]+c1[j+np*((i+1)/2+1)].pos[1])/8;
  c2[j+i*np].pos[2]=(c1[j+np*((i+1)/2-1)].pos[2]+6*c1[j+np*(i+1)/2].pos[2]+c1[j+np*((i+1)/2+1)].pos[2])/8;
 }
}
 for(int i=2;i<np-1;i=i+2)
{
  for(int j=0; j<np;j++)
 {
  c2[j+i*np].pos[0]=(c1[j+(i/2*np)].pos[0]+c1[j+(i/2+1)*np].pos[0])/2;
  c2[j+i*np].pos[1]=(c1[j+(i/2*np)].pos[1]+c1[j+(i/2+1)*np].pos[1])/2;
  c2[j+i*np].pos[2]=(c1[j+(i/2*np)].pos[2]+c1[j+(i/2+1)*np].pos[2])/2;
 }
}
  //i=np-1 e.g.143
for(int i=0; i<np*np; i++)
{
 c2[i].pos[3]=1.0f; 
 c2[i].color=glm::vec4(0.5f,0.5f,0.5f,1.0f);
}

}

Vertex subdivision11[2*144];
Vertex subdivision12[2*2*144];
Vertex subdivision21[2*4*144];
Vertex subdivision22[4*4*144];
Vertex subdivision31[4*8*144];
Vertex subdivision32[8*8*144];

void subdivision_level()
{
  initialize_subdivision1(control, subdivision11,24);
  initialize_subdivision2(subdivision11,subdivision12, 24);
  initialize_subdivision1(subdivision12,subdivision21, 48);
  initialize_subdivision2(subdivision21,subdivision22, 48);
  initialize_subdivision1(subdivision22,subdivision31, 96);
  initialize_subdivision2(subdivision31,subdivision32, 96);
}

Vertex surface[3*2*(96-1)*(96-1)];
void surface_triangles()
{ 
int np=96;
for(int i=0; i<8*8*12*12;i++)
{
 subdivision32[i].normal=glm::vec4(0,0,1,1.0f);
}
 for(int i=0;i<(np-1);i++)
{
 for(int j=0; j<np-1;j++)
{
 surface[(j+i*(np-1))*3]=subdivision32[(j+i*np)*3];
 surface[(j+i*(np-1))*3+1]=subdivision32[(j+i*np)*3+1];
 surface[(j+i*(np-1))*3+2]=subdivision32[(j+(i+1)*np)*3];
}
}
 for(int i=0;i<(np-1);i++)
{
 for(int j=0; j<np-1;j++)
{
 surface[(j+i*(np-1))*3+2*(96-1)*(96-1)]=subdivision32[(j+i*np)*3+1];
 surface[(j+i*(np-1))*3+1+2*(96-1)*(96-1)]=subdivision32[(j+i*np)*3];
 surface[(j+i*(np-1))*3+2+2*(96-1)*(96-1)]=subdivision32[(j+(i+1)*np)*3+1];
}
}


}


*/
void initialize_subdivision1(Vertex c1[], Vertex c2[], int np){
 //Pnew2i= (Poldi-1+6Poldi+Poldi+1)/8;
 //Pnew2i+1=(4Poldi+4Poldi+1)/8=(Poldi+Poldi+1)/2;
 //i=0
//np=24 for the first time
for(int i=0;i<np/2;i++)
{
 for(int j=0; j<np/2;j++)
{
  c2[(np-1)*i+2*j].pos=c1[np/2*i+j].pos;

///////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  c2[(np-1)*i+2*j].tex_coord=c1[np/2*i+j].tex_coord;
}
}

  for(int i=0;i<np/2;i++)
{
  for(int j=0; j<np/2-1;j++)
{
  c2[(np-1)*i+2*j+1].pos=(c1[np/2*i+j].pos+c1[np/2*i+j+1].pos)/2;
///////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  c2[(np-1)*i+2*j+1].tex_coord=(c1[np/2*i+j].tex_coord+c1[np/2*i+j+1].tex_coord)/2.0f;
}
} 

  //i=np-1 e.g.143
for(int i=0; i<(np-1)*(np/2); i++)
{
 c2[i].pos[3]=1.0f; 
 c2[i].color=glm::vec4(0.5f,0.5f,0.5f,1.0f);
}
}

void initialize_subdivision2(Vertex c1[], Vertex c2[], int np)
{

for(int i=0;i<np/2;i++)
{
 for(int j=0; j<np-1;j++)
{
  c2[(np-1)*2*i+j].pos=c1[(np-1)*i+j].pos;
///////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  c2[(np-1)*2*i+j].tex_coord=c1[(np-1)*i+j].tex_coord;
}
}

  for(int i=0;i<np/2-1;i++)
{
  for(int j=0; j<np-1;j++)
{
 c2[(np-1)*2*i+(np-1)+j].pos=(c1[(np-1)*i+j].pos+c1[(np-1)*(i+1)+j].pos)/2;
///////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
 c2[(np-1)*2*i+(np-1)+j].tex_coord=(c1[(np-1)*i+j].tex_coord+c1[(np-1)*(i+1)+j].tex_coord)/2.0f;
}
}

  //i=np-1 e.g.143
for(int i=0; i<(np-1)*(np-1); i++)
{
 c2[i].pos[3]=1.0f; 
 c2[i].color=glm::vec4(0.5f,0.5f,0.5f,1.0f);
}

}

Vertex subdivision11[12*23];
Vertex subdivision12[23*23];
Vertex subdivision21[23*45];
Vertex subdivision22[45*45];
Vertex subdivision31[89*45];
Vertex subdivision32[89*89];

void subdivision_level()
{
  initialize_subdivision1(control, subdivision11,24);
  initialize_subdivision2(subdivision11,subdivision12, 24);
  initialize_subdivision1(subdivision12,subdivision21, 46);
  initialize_subdivision2(subdivision21,subdivision22, 46);
  initialize_subdivision1(subdivision22,subdivision31, 90);
  initialize_subdivision2(subdivision31,subdivision32, 90);
}

Vertex surface[3*2*(89-1)*(89-1)];
void surface_triangles()
{ 
int np=89;
for(int i=0; i<3*2*(89-1)*(89-1);i++)
{
 surface[i].normal=glm::vec4(0,0,1,1.0f);
 surface[i].color=glm::vec4(0.5,0.5,0.5,1.0f);
}
 for(int i=0;i<np-1;i++)
{
 for(int j=0; j<np-1;j++)
{
 surface[(j*3+3*i*(np-1))].pos=subdivision32[j+i*np].pos;
 surface[(j*3+3*i*(np-1))+1].pos=subdivision32[j+i*np+1].pos;
 surface[(j*3+3*i*(np-1))+2].pos=subdivision32[j+(i+1)*np].pos;


///////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  surface[(j*3+3*i*(np-1))].tex_coord=subdivision32[j+i*np].tex_coord;
  surface[(j*3+3*i*(np-1))+1].tex_coord=subdivision32[j+i*np+1].tex_coord;
  surface[(j*3+3*i*(np-1))+2].tex_coord=subdivision32[j+(i+1)*np].tex_coord;
}
}
 for(int i=0;i<np-1;i++)
{
 for(int j=0; j<np-1;j++)
{
 surface[(j*3+3*i*(np-1))+3*(89-1)*(89-1)].pos=subdivision32[j+i*np+1].pos;
 surface[(j*3+3*i*(np-1))+1+3*(89-1)*(89-1)].pos=subdivision32[j+(i+1)*np].pos;
 surface[(j*3+3*i*(np-1))+2+3*(89-1)*(89-1)].pos=subdivision32[j+(i+1)*np+1].pos;


///////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT 
surface[(j*3+3*i*(np-1))+3*(89-1)*(89-1)].tex_coord=subdivision32[j+i*np+1].tex_coord;
 surface[(j*3+3*i*(np-1))+1+3*(89-1)*(89-1)].tex_coord=subdivision32[j+(i+1)*np].tex_coord;
 surface[(j*3+3*i*(np-1))+2+3*(89-1)*(89-1)].tex_coord=subdivision32[j+(i+1)*np+1].tex_coord;
}
}


}
int color_code_pick(int x,int y){
    
    // draw the scene
        glClearColor(0.0f,0.0f,0.0f,1.0f);
 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 	float cc = 0.08f;
	float dd = 0.08f;
	for(int i=0;i<12;i++)
	{
  	float cc=0.08f;
  	for(int j=0;j<12;j++)
  	{
  	control[j+i*12].color=glm::vec4(cc,dd,1,1);
  	cc=cc+0.08f;
  	}
  	dd=dd+0.08f;
	}
 	draw_grid();
 	glFlush();

  // read back the green channel of pixel under the cursor into data
  // data is in range [0,255]
  	GLubyte gdata;
        GLubyte rdata;
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	glReadPixels(x,viewport[3]-y,1,1,GL_RED,GL_UNSIGNED_BYTE,&rdata);
        glReadPixels(x,viewport[3]-y,1,1,GL_GREEN,GL_UNSIGNED_BYTE,&gdata);
        printf("color b: %d\n",rdata);
        printf("color d: %d\n",gdata);

    // TODO: identify the selected point and return its index
        if(gdata>0)
      	{
	float a;
	a=(float)rdata/(20.5);
	int b= ceil(a);
        float c;
        c=(float)gdata/(20.5);
        int d =ceil(c);
	return (d-1)*12+b-1;
	}
	else 
       {
	return -1;
	}
}
// MOUSE handling *******************************************
int last_x, last_y, last_z;
int selected_idx = -1;

void mouse(int button, int state, int x, int y ){
	if(state == GLUT_DOWN){
        selected_idx = color_code_pick(x,y);
	printf("pick point #%d\n", selected_idx);
	}else{
        
	}

    last_x = x;
    last_y = y;
    
}


void motion( int x, int y){
    GLint viewport[4];
    float dx,dy;
    
    glGetIntegerv(GL_VIEWPORT,viewport);
    
    dx = -(x - last_x)/(float)viewport[2] * 2;
    dy = -(y - last_y)/(float)viewport[3] * 2;
    
    last_x = x;
    last_y = y;
if(selected_idx>=0){
    control[selected_idx].pos[0]=control[selected_idx].pos[0]+dx;   
    control[selected_idx].pos[1] = control[selected_idx].pos[1]+dy;
    control[selected_idx].pos[2] = control[selected_idx].pos[2]+dx;
}
    subdivision_level();
    surface_triangles();
    glutPostRedisplay();
    
}


void lighting()
{
        glm::vec4 light = glm::vec4(0.0f,10.0f,0.0f,1.0f);
	glm::vec4 diffuse = glm::vec4(1,1,1,1);
        glm::vec4 diffusem=glm::vec4(0.8,0.8,0.8,1);
	glm::vec4 specular = glm::vec4(1,1,1,1);
        glm::vec4 specularm = glm::vec4(0.8,0.8,0.8,1);
        float diffusex=diffuse.x*diffusem.x;
        float diffusey=diffuse.y*diffusem.y;
        float diffusez=diffuse.z*diffusem.z;
        float specularx=specular.x*specularm.x;
        float speculary=specular.y*specularm.y;
        float specularz=specular.z*specularm.z;
   for(int i = 0; i<2*3*88*88; i++){

	float diffuse_intensity = max(glm::dot(glm::normalize(light - surface[i].pos), surface[i].normal),0.0f);
	glm::vec4 s = glm::normalize(light - surface[i].pos) + glm::normalize(glm::vec4(10.0f,10.0f,10.0f,1.0f) - surface[i].pos);
 	float specular_intensity = max(glm::dot(glm::normalize(s),surface[i].pos), 0.0f);
	float R = surface[i].color.x * ( diffuse_intensity * diffusex + specular_intensity * specularx);
	float G = surface[i].color.y * ( diffuse_intensity * diffusey + specular_intensity * speculary);
	float B = surface[i].color.z * ( diffuse_intensity * diffusez + specular_intensity * specularz);
	surface[i].color = glm::vec4(R,G,B,1.0f);
		}

}

void initial_color(){
for(int i = 0; i<2*3*88*88; i++){
surface[i].color=glm::vec4(0.5,0.5,0.5,1.0f);
}

}
// KEYBOARD handling *******************************************
int controlp=0;
int surfacep=0;
int light=0;
int texture=0;
void keyboard(unsigned char key, int x, int y)
{
    switch (key) {
    case 'r':
        graphics_init();
        glutPostRedisplay();
        break;
    case 27:
        exit(0);
        break;      
     break;
    case 'd':
        if(controlp==0)
        {  controlp=1;
           cout<<"show control points"<<endl;
        }
	else if (controlp==1)
 	{
	   controlp=0;
	   cout<<"hide control points"<<endl;
	}
     break;
    case 's':
        if(surfacep==0)
    	{
 	   surfacep=1;
           cout<<"show surface"<<endl;
     	}
	else if(surfacep==1)
	{
	   surfacep=0;
	   cout<<"hide surface"<<endl;
	}
     
     break;
     case 'l':
        if(light==0)
       {  
	  light=1;
          cout<<"show light"<<endl;
       }
        else if(light==1)
       {
	  light=0;
          cout<<"hide light"<<endl;
       }
     break;
     case 't':
        if(texture==0)
       {  
	  texture=1;
          cout<<"show texture"<<endl;
       }
        else if(texture==1)
       {
	  texture=0;
          cout<<"hide texture"<<endl;
       }
     break;
     }
    glutPostRedisplay();
}

float cx=4.1f; float cy=0.8f;
glm::mat4 newc;
void special_key(int key, int x, int y)
{
    // TODO: capture the arrow key here
//add move camera   
 if(key==GLUT_KEY_LEFT)
          {
              cy=cy+0.05f;
           
	  }
 	  else if(key==GLUT_KEY_RIGHT)
          {
	      cy=cy-0.05f;
	  }
 	  else if(key==GLUT_KEY_UP)
          {
              cx=cx+0.05f;
	  }
          else if(key==GLUT_KEY_DOWN)
          {
 	      cx=cx-0.05f;
	  }
        
    glutPostRedisplay();

}

// DISPLAY and RENDERING functions *************************

float rotate_x; float rotate_y; float rotate_z;
void draw_grid(){
    
    GLuint vao;
    GLuint grid_vbo;
    
    GLuint position_location;
    GLuint color_location;
    GLuint normal_location;

//////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    GLuint texcoord_location;
    GLuint MVP_location;

    GLuint control_vbo;
    GLuint surface_vbo;

    rotate_x=-cos(cx)*sin(cy)*10.0f;
    rotate_y=sin(cx)*10.0f;
    rotate_z=cos(cx)*cos(cy)*10.0f;
    

//add
//the camera is sitting at (10,10,10) and looking at (0,0,0), the up direction is (0,1,0)
    
    glm::vec3 POS=glm::vec3(rotate_x,rotate_y,rotate_z);
    glm::vec3 AT=glm::vec3(0.0f,0.0f,0.0f);
    glm::vec3 UP=glm::vec3(0.0f,1.0f,0.0f);
    glm::mat4 view = glm::lookAt(POS,AT,UP); 
    glm::mat4 projection=glm::perspective(45.0f,4.0f/3.0f,0.1f,100.0f);
    glm::mat4 model =glm::mat4(1.0f);
    glm::mat4 MVP = projection*view*model;
    

    // specify the shaders we want to use
    glUseProgram(program);

    // get the input variable location in this shader
    position_location = glGetAttribLocation(program, "in_position");
    color_location = glGetAttribLocation(program, "in_color");
    normal_location =glGetAttribLocation(program, "in_normal");

///////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    texcoord_location=glGetAttribLocation(program, "in_texcoord");


    MVP_location = glGetUniformLocation(program, "MVP");    

    // create and bind a VAO
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // enable the input locations we wan to use
    glEnableVertexAttribArray(position_location);
    glEnableVertexAttribArray(color_location);
    glEnableVertexAttribArray(normal_location);
//////////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    glEnableVertexAttribArray(texcoord_location);

//add   
    glUniformMatrix4fv(MVP_location, 1, GL_FALSE, &MVP[0][0]);

    // draw points

    // generate and bind a vertex buffer to hold the vertex data on GPU
    glGenBuffers(1, &grid_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, grid_vbo);

    // initialize the vertex buffer with the vertex data
    
    glBufferData(GL_ARRAY_BUFFER, grid_vertices.size() * sizeof(Vertex), &grid_vertices[0] , GL_STATIC_DRAW);

    // specify the input data format
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
    glVertexAttribPointer(texcoord_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));
     

/////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    glUniform1i(glGetUniformLocation(program, "face_tex"), 0);
    
   // draw points
    glPointSize(10);
    glDrawArrays(GL_LINES, 0, grid_vertices.size());

if(controlp==1)
{
    glGenBuffers(1, &control_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, control_vbo);
    glBufferData(GL_ARRAY_BUFFER, 12*12*sizeof(Vertex), &control[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
    glVertexAttribPointer(texcoord_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));
    glPointSize(5);
    glDrawArrays(GL_POINTS, 0, 12*12);
}

if(surfacep==1)
{
    glGenBuffers(1, &surface_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, surface_vbo);
    glBufferData(GL_ARRAY_BUFFER, 2*3*88*88*sizeof(Vertex), &surface[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
    glVertexAttribPointer(texcoord_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));
    glPointSize(5);
    glDrawArrays(GL_TRIANGLES, 0, 2*3*88*88);
}

    // unbind VAO and VBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Delete VAO and VBO
    glDeleteBuffers(1, &control_vbo);
    glDeleteBuffers(1, &surface_vbo);

    glDeleteBuffers(1, &grid_vbo);
    glDeleteVertexArrays(1, &vao);  
}

void display(){
    // Clear Viewport
    glClearColor(0.0f,0.0f,0.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for(int i =0; i<144; i++)
    {
    control[i].color=glm::vec4(0,1,0,1); 
    }
    
    if(light==1){
    lighting();
    }
    else if(light==0)
    {initial_color();
    }


    draw_grid();

    glFlush();
    glutSwapBuffers();

}

void reshape(int width, int height){
    // Clip the view port to match our ratio
    glViewport(0, 0, width, height);
    glutPostRedisplay();
}

void graphics_init(){
    cx=2.9f;
    cy=0.8f;

    // init vertex shader
    read_shader_source_code("vs.glsl", buf, 1024);
    cout << buf << endl;
    int vertex_shader_source_length = strlen(buf);
    const char *vertex_shader_sources[] = { buf };

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_sources, &vertex_shader_source_length);
    glCompileShader(vertex_shader);
    
    // init fragment shader 
    read_shader_source_code("fs.glsl", buf, 1024); 
    int fragment_shader_source_length = strlen(buf);
    const char *fragment_shader_sources[] = { buf };
    
    fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, fragment_shader_sources, &fragment_shader_source_length);
    glCompileShader(fragment_shader);

    // init program
    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);

    glLinkProgram(program);
  
    // enable depth test
    glEnable(GL_DEPTH_TEST);


///////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
   
    GLuint texture_id = load_texture_TGA("weilin.tga", NULL, NULL, GL_CLAMP, GL_CLAMP);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture_id);
 
}

void print_info()
{
    fprintf(
        stdout,
        "INFO: OpenGL Version: %s\n",
        glGetString(GL_VERSION)
        );
    fprintf(stdout, "INFO: Using GLEW %s\n", glewGetString(GLEW_VERSION));

    if(glewIsSupported("GL_ARB_debug_output"))
    {
        fprintf(stdout, "INFO: Support ARB_debug_output\n");
    }
    else
    {
        fprintf(stdout, "INFO: Not support ARB_debug_output\n");
    }
}

int main(int argc,char * argv[]){

    // Setup GLUT
    glutInit(&argc,argv);

    glutInitContextVersion (3, 0);
    glutInitContextFlags (GLUT_DEBUG);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB     | GLUT_DEPTH | GLUT_STENCIL);
    glutInitWindowSize (window_width,window_height);
    glutCreateWindow("First step - OpenGL 3");
  
    // Setup GLEW
    glewExperimental = true;
    GLenum err = glewInit();
    if(GLEW_OK != err)
    {
        /* Problem: glewInit failed, something is seriously wrong. */
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
    }

    print_info();

    if(glewIsSupported("GL_ARB_debug_output"))
        VSDebugLib::init();
    // Initialize OpenGL
    graphics_init ();

    initialize_points();
    initialize_controlpoints();
    subdivision_level();
    surface_triangles();

    // set up callbacks
    glutReshapeFunc(&reshape);
    glutDisplayFunc(&display);
    glutKeyboardFunc(&keyboard);
    glutMouseFunc(&mouse);
    glutMotionFunc(&motion);
    glutSpecialFunc(&special_key);

    // main loop
    glutMainLoop();
    return 0;
}
