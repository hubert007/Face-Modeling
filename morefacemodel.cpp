/*
 * Computer Graphics Fall 2013
 * 
 * Programmers:
 * - Melike Selin Aydin
 * - Wei Lin
 * - Icaro Alzuru
 *  
 * */   

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>


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


// CONSTANTS
const float PI = 3.1415926;
const int window_width = 500, window_height = 500;
const char * window_title = "HW3";
const int ptArrayLength = 60;

float camx = 50.0;
float camy = 50.0;
float camz = 50.0;
float degX = 0.0;
float degY = 0.0;
float upx = 0;
float upy = 1.0;
float upz = 0.0;
bool isDisc = false;
float isLighting = 0.1;
float isTexture = 0.0;

bool showPoints = true;
bool showSurface = true;
bool showMesh = false;
bool showAnimation = false;
int  showhair=0;
bool showEyes = false;


glm::vec4 light_loc = glm::vec4(10.0, 10.0, 10.0, 1.0);
glm::vec4 light_int = glm::vec4(1.0, 1.0, 1.0, 1.0);

struct Vertex
{
    glm::vec4 pos;
    glm::vec4 color;
    glm::vec4 normal;
    glm::vec2 tex_coord;
    int neighbor_count;
    Vertex(){ normal = glm::vec4(1.0, 1.0, 1.0, 1.0); 
    tex_coord = glm::vec2(0.0, 0.0); neighbor_count=0;};
};

void draw_grid();

// runtime variable
char buf[2048];
vector<Vertex> grid_vertices;
Vertex pt_array[ptArrayLength][ptArrayLength];

Vertex top[10][360];
Vertex top_triangles[3*2*360*19];

vector<Vertex> sub_surface;
vector<Vertex> temp;
vector<Vertex> mesh;
vector<Vertex> drawMesh;

vector<Vertex> eyeLeft;
vector<Vertex> eyeRight;
glm::vec4 animation[11][4];

GLuint vertex_shader, fragment_shader, program;


glm::mat4 view_matrix, projection_matrix, model_matrix, mvp_matrix, modelview_matrix;
glm::mat3 normal_mat;

// function Prototype 
void graphics_init();

int ind = 0;

void read_animation(){
	ifstream myfile;
	myfile.open("smile.txt");
	string line;

	float x, y, z;
	
	if (myfile.is_open())
	{		
		
		for (int i = 0; i < 11; i++)
		{
			getline(myfile, line);
			for (int j = 0; j < 4; j++)
			{
				myfile >> x >> y >> z;
				animation[i][j] = glm::vec4(x, y, z, 1.0);
			}		
		}
		myfile.close();
	}
	
	
}
void init_pt_array(){
	read_animation();
	float a, b, c;
	a = 0;
	b = -.2;
	c = 2.2;
	for (int i = 0; i < ptArrayLength; i++)
	{
		a = 0;
		for (int j = 0; j < ptArrayLength-1; j++)
		{
			pt_array[i][j].pos[0]=cos(a*PI/180)*c;
			pt_array[i][j].pos[1]=b;
			pt_array[i][j].pos[2]=sin(a*PI/180)*c;
			pt_array[i][j].pos[3] = 1.0;
	        a=a+(360/(ptArrayLength-1));
			pt_array[i][j].color = glm::vec4(((float)i/(float)ptArrayLength)* 0.8 + 0.2, ((float)j/(float)ptArrayLength)* 0.8 + 0.2, 0.1, 1.0);
			//~ pt_array[i][j].color = glm::vec4(0.1, 0.2, 1.0, 1.0);
			
			pt_array[i][j].tex_coord = glm::vec2(((float)i/(float)ptArrayLength)-.2, ((float)j/(float)ptArrayLength));
			//~ cout << " color at first: " << 255*(((float)i/(float)ptArrayLength)* 0.8 + 0.2) << " " << 255*(((float)j/(float)ptArrayLength)* 0.8 + 0.2) << endl;
		}
		if(int j=ptArrayLength-1)
		{
			pt_array[i][j].pos[0]=cos(0)*c;
			pt_array[i][j].pos[1]=b;
			pt_array[i][j].pos[2]=sin(0)*c;
			pt_array[i][j].pos[3] = 1.0;
			pt_array[i][j].color = glm::vec4(((float)i/(float)ptArrayLength)* 0.8 + 0.2, ((float)j/(float)ptArrayLength)* 0.8 + 0.2, 0.1, 1.0);
			//~ pt_array[i][j].color = glm::vec4(0.1, 0.2, 1.0, 1.0);
			
			pt_array[i][j].tex_coord = glm::vec2(((float)i/(float)ptArrayLength)-.2, ((float)j/(float)ptArrayLength));
		}
		
		b+=0.06;
	}	

		for (int i = 0; i <20; i++)
	{
		float c1=1.6;
		a = 0;
		for (int j = 0; j < ptArrayLength-1; j++)
		{
			pt_array[i][j].pos[0]=cos(a*PI/180)*c1;
			pt_array[i][j].pos[2]=sin(a*PI/180)*c1;
	        a=a+(360/(ptArrayLength-1));
		}
		if(int j=ptArrayLength-1)
		{
			pt_array[i][j].pos[0]=cos(0)*c1;
			pt_array[i][j].pos[2]=sin(0)*c1;
		}
	}	
    float c2=1.7;
		for (int i = 20; i <26; i++)
	{
		
		a = 0;
		for (int j = 0; j < ptArrayLength-1; j++)
		{
			pt_array[i][j].pos[0]=cos(a*PI/180)*c2;
			pt_array[i][j].pos[2]=sin(a*PI/180)*c2;
	        a=a+(360/(ptArrayLength-1));
		}
		if(int j=ptArrayLength-1)
		{
			pt_array[i][j].pos[0]=cos(0)*c2;
			pt_array[i][j].pos[2]=sin(0)*c2;
		}
		c2=c2+0.1;
	}	
	float c3=2.5;
		for (int i = 26; i <42; i++)
	{
		
		a = 0;
		for (int j = 0; j < ptArrayLength-1; j++)
		{
			pt_array[i][j].pos[0]=cos(a*PI/180)*c3;
			pt_array[i][j].pos[2]=sin(a*PI/180)*c3;
	        a=a+(360/(ptArrayLength-1));
		}
		if(int j=ptArrayLength-1)
		{
			pt_array[i][j].pos[0]=cos(0)*c3;
			pt_array[i][j].pos[2]=sin(0)*c3;
		}
		c3=c3+0.01;
	}	
float c4=2.6;
		for (int i = 42; i <60; i++)
	{
		
		a = 0;
		for (int j = 0; j < ptArrayLength-1; j++)
		{
			pt_array[i][j].pos[0]=cos(a*PI/180)*c4;
			pt_array[i][j].pos[2]=sin(a*PI/180)*c4;
	        a=a+(360/(ptArrayLength-1));
		}
		if(int j=ptArrayLength-1)
		{
			pt_array[i][j].pos[0]=cos(0)*c4;
			pt_array[i][j].pos[2]=sin(0)*c4;
		}
		c4=c4-0.024;
	}	
}

	

void createtop(){
	float a=0;
	float b=0;
	float c=0;
	for(int i=0;i<10;i++)
	{
		for(int j=0; j<360;j++)
		{
		top[i][j].pos[0]=cos(a*PI/180)*2.199*(1-c*c);
		top[i][j].pos[1]=3.3+b;
		top[i][j].pos[2]=sin(a*PI/180)*2.199*(1-c*c);
		top[i][j].pos[3]=1.0f;
		top[i][j].color=glm::vec4(0.0,1.0,0.0,1.0);
		a=a+1;
	}
	c=(c-0.100);
	b=b+0.1;
	}
 
}

void create_toptriangles()
{
	for(int i=0; i<19;i++)
	{
		for(int j=0; j<359;j++)
		{
			top_triangles[3*(i*360+j)]=top[i][j];
			top_triangles[3*(i*360+j)+1]=top[i][j+1];
			top_triangles[3*(i*360+j)+2]=top[i+1][j];
			top_triangles[3*360*19+3*(i*360+j)]=top[i][j+1];
			top_triangles[3*360*19+3*(i*360+j)+1]=top[i+1][j];
			top_triangles[3*360*19+3*(i*360+j)+2]=top[i+1][j+1];
			
		}
		if(int j=359)
		{
			top_triangles[3*(i*360+j)]=top[i][j];
			top_triangles[3*(i*360+j)+1]=top[i][0];
			top_triangles[3*(i*360+j)+2]=top[i+1][j];
			top_triangles[3*360*19+3*(i*360+j)]=top[i][0];
			top_triangles[3*360*19+3*(i*360+j)+1]=top[i+1][j];
			top_triangles[3*360*19+3*(i*360+j)+2]=top[i+1][0];
			
		}
	}
	for(int i=0; i<3*2*360*19;i++)
	{top_triangles[i].color=glm::vec4(214.0/255.0,135.0/255.0,88.0/255.0,1.0);
	}
}

Vertex toptop[360*3];
void create_toptopt()
{
	for(int i=0; i<359;i++)
	{
	toptop[3*i].pos=glm::vec4(0.0,4.5,0.0,1.0);
	toptop[3*i].color=glm::vec4(214.0/255.0,135.0/255.0,88.0/255.0,1.0);
	toptop[3*i+1]=top[19][i];
	toptop[3*i+2]=top[19][i+1];	
	}
	if(int i=359)
	{
	toptop[3*i].pos=glm::vec4(0.0,4.5,0.0,1.0);
	toptop[3*i].color=glm::vec4(214.0/255.0,135.0/255.0,88.0/255.0,1.0);
	toptop[3*i+1]=top[19][i];
	toptop[3*i+2]=top[19][0];		
	}
}


Vertex hair[360*20*3];
void initialize_hair()
{
for(int i=0; i<20;i++)
{
for(int j=0; j<359;j++)
{
hair[3*(i*360+j)].pos=top[i][j].pos;
hair[3*(i*360+j)+2].pos=top[i][j+1].pos;
float a=top[i][j].pos[0]*1.1;
float b=top[i][j].pos[2]*1.1;
hair[3*(i*360+j)+1].pos=glm::vec4(a,4.6f,b,1.0f);
}
if(int j=359)
{
hair[3*(i*360+j)].pos=top[i][j].pos;
hair[3*(i*360+j)+2].pos=top[i][0].pos;
float a=top[i][j].pos[0]*1.1;
float b=top[i][j].pos[2]*1.1;
hair[3*(i*360+j)+1].pos=glm::vec4(a,4.6f,b,1.0f);
}
}
for(int i=0; i<20*360*3;i++)
{
hair[i].color=glm::vec4(78.0/255.0,65.0/255.0,191.0/255.0,1.0f);
}
}


void initialize_rehair()
{
for(int i=0; i<20;i++)
{
for(int j=0; j<359;j++)
{
hair[3*(i*360+j)].pos=top[i][j].pos;
hair[3*(i*360+j)+2].pos=top[i][j+1].pos;
float a=top[i][j].pos[0]*1.5;
float b=top[i][j].pos[2]*1.5;
hair[3*(i*360+j)+1].pos=glm::vec4(a,5.6f,b,1.0f);
}
if(int j=359)
{
hair[3*(i*360+j)].pos=top[i][j].pos;
hair[3*(i*360+j)+2].pos=top[i][0].pos;
float a=top[i][j].pos[0]*1.5;
float b=top[i][j].pos[2]*1.5;
hair[3*(i*360+j)+1].pos=glm::vec4(a,5.6f,b,1.0f);
}
}
for(int i=0; i<20*360*3;i++)
{
hair[i].color=glm::vec4(78.0/255.0,65.0/255.0,191.0/255.0,1.0f);
}
}


void create_sub_surface(){
	temp.clear();
	sub_surface.clear();
	
	//n*2n matrix olusumu -> temp
	for (int i = 0; i < ptArrayLength; i++)
	{
		for (int j = 0; j < ptArrayLength; j++)
		{
			Vertex a, b;
			a.color = glm::vec4(0.1, 0.1, 0.1, 1.0);
			b.color = glm::vec4(0.1, 0.1, 0.1, 1.0);
			a.pos =(pt_array[i][max(j-1, 0)].pos + (pt_array[i][j].pos) + pt_array[i][min(j+1, ptArrayLength-1)].pos) / 3;			
			a.tex_coord =(pt_array[i][max(j-1, 0)].tex_coord + (1.0f * pt_array[i][j].tex_coord) + pt_array[i][min(j+1, ptArrayLength-1)].tex_coord) / 3.0f;
			b.pos = ((1 * pt_array[i][j].pos) + (1 * pt_array[i][min(j+1, ptArrayLength-1)].pos)) / 2;
			b.tex_coord = ((1.0f * pt_array[i][j].tex_coord) + (1.0f * pt_array[i][min(j+1, ptArrayLength-1)].tex_coord)) / 2.0f;
			temp.push_back(a);
			temp.push_back(b);
		}	
	}
	int sz = 2*ptArrayLength;
	for (int j = 0; j < sz; j++)
	{
		for (int i = 0; i < ptArrayLength; i++)
		{			
			Vertex a, b;
			a.color = glm::vec4(0.1, 0.1, 0.1, 1.0);
			b.color = glm::vec4(0.1, 0.1, 0.1, 1.0);
			a.pos = (temp[max(i-1, 0)*sz + j].pos + (1*temp[i*sz+j].pos) + temp[min(i+1, ptArrayLength-1)*sz + j].pos) / 3;
			a.tex_coord = (temp[max(i-1, 0)*sz + j].tex_coord + (1.0f*temp[i*sz+j].tex_coord) + temp[min(i+1, ptArrayLength-1)*sz + j].tex_coord) / 3.0f;
			b.pos = ((1*temp[i*sz+j].pos) + (1*temp[min(i+1, ptArrayLength-1)*sz+j].pos)) / 2;
			b.tex_coord = ((1.0f*temp[i*sz+j].tex_coord) + (1.0f*temp[min(i+1, ptArrayLength-1)*sz+j].tex_coord)) / 2.0f;
			sub_surface.push_back(a);
			sub_surface.push_back(b);
		}
		
	}
	//4n*4n matrix olusumu
	for (int i = 0; i < sz; i++)
	{
		for (int j = 0; j < sz; j++)
		{
			Vertex a, b;
			a.color = glm::vec4(0.1, 0.1, 0.1, 1.0);
			b.color = glm::vec4(0.1, 0.1, 0.1, 1.0);
			a.pos =(sub_surface[i*sz+max(j-1, 0)].pos + (1 * sub_surface[i*sz+j].pos) + sub_surface[i*sz+min(j+1, ptArrayLength-1)].pos) / 3;
			a.tex_coord =(sub_surface[i*sz+max(j-1, 0)].tex_coord + (1.0f * sub_surface[i*sz+j].tex_coord) + sub_surface[i*sz+min(j+1, ptArrayLength-1)].tex_coord) / 3.0f;
			b.pos = ((1 * sub_surface[i*sz+j].pos) + (1 * sub_surface[i*sz+min(j+1, ptArrayLength-1)].pos)) / 2;
			b.tex_coord = ((1.0f * sub_surface[i*sz+j].tex_coord) + (1.0f * sub_surface[i*sz+min(j+1, ptArrayLength-1)].tex_coord)) / 2.0f;
			temp.push_back(a);
			temp.push_back(b);
		}		
	}

	
	for (int i = 0; i < sub_surface.size(); i++)
	{
		sub_surface[i].pos[3] = 1.0;
	}
	
}

void create_rects(){
	temp.clear();
	int sz = round(sqrt(sub_surface.size()));
	for (int i = 0; i < sz; i++)
	{
		for (int j = 0; j < sz; j++)
		{
			temp.push_back(sub_surface[i*sz+j]);
			temp.push_back(sub_surface[(min(i+1, sz-1)*sz)+j]);
		}		
	}
	
}

int color_code_pick_green(int x,int y){
    
    // draw the scene
	glClearColor(0.0f,0.0f,0.0f,1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw_grid();
	glFlush();

    // read back the green channel of pixel under the cursor into data
    // data is in range [0,255]
	GLubyte data;
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	glReadPixels(x,viewport[3]-y,1,1,GL_GREEN,GL_UNSIGNED_BYTE,&data);
    printf("green color: %d\n", data);
    
    if(data == 0)
		return -1;
	else
	{
		float f_data = ((float)data) / 255.0;
		cout << "f_data: " << f_data << endl;
		int index = round((f_data - 0.2) / 0.8 * ptArrayLength);
		
		return index;
	}
}

int color_code_pick_red(int x,int y){
    
    // draw the scene
	glClearColor(0.0f,0.0f,0.0f,1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw_grid();
	glFlush();

    // read back the red channel of pixel under the cursor into data
    // data is in range [0,255]
	GLubyte data;
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	glReadPixels(x,viewport[3]-y,1,1,GL_RED,GL_UNSIGNED_BYTE,&data);
    printf("red color: %d\n", data);

    // identify the selected point and return its index
    // 
    
    if(data == 0)
		return -1;
	else
	{
		float f_data = (float)data / 255.0;
		int index = round((f_data - 0.2) / 0.8 * ptArrayLength);
		
		return index;
	}
}

glm::vec4 mycross(glm::vec4 a, glm::vec4 b){
	glm::vec4 ret(a[1]*b[2], a[2]*b[0], a[0]*b[1], 1.0);
	return ret;
}
glm::vec4 mynormalize(glm::vec4 a){
	float size = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	return (glm::vec4(a[0]/size, a[1]/size, a[2]/size, 1.0));
}

void read_mesh(){
	ifstream myfile;
	myfile.open("icaro.obj");
	
	if (myfile.is_open())
	{
		char temp;
		float t1, t2, t3;
		mesh.clear();
		myfile >> temp >> t1 >> t2 >> t3;
		while (temp == 'v')
		{
			Vertex a;
			a.pos = glm::vec4(t1, t2, t3, 1.0);
			a.color = glm::vec4(0.5, 0.5, 0.5, 1.0);
			mesh.push_back(a);
			myfile >> temp >> t1 >> t2 >> t3;
		}
		int i1, i2, i3;
		i1 = (int)t1 -1;
		i2 = (int)t2 -1;
		i3 = (int)t3 -1;
		while (myfile.good())
		{
			//do the stuff
			glm::vec4 tempVec1 = mesh[i2].pos - mesh[i1].pos;
			glm::vec4 tempVec2 = mesh[i3].pos - mesh[i1].pos;
			mesh[i1].normal = mesh[i1].normal + mycross(tempVec1, tempVec2);
			mesh[i1].neighbor_count++;
			
			tempVec1 =  mesh[i3].pos - mesh[i2].pos;
			tempVec2 =  mesh[i1].pos - mesh[i2].pos;
			mesh[i2].normal = mesh[i2].normal + mycross(tempVec1, tempVec2);
			mesh[i2].neighbor_count++;
			
			tempVec1 =  mesh[i2].pos - mesh[i3].pos;
			tempVec2 =  mesh[i1].pos - mesh[i3].pos;
			mesh[i3].normal = mesh[i3].normal + mycross(tempVec1, tempVec2);
			mesh[i3].neighbor_count++;
			
			
			//en sonunda (1, 1, 1, 0) i cikarip neighborcount a bolup normalize ediceksn.
			
			myfile >> temp >> i1 >> i2 >> i3;
			i1-=1;
			i2-=1;
			i3-=1;
		}

		for (int i = 0; i < mesh.size(); i++)
		{
			mesh[i].normal = mesh[i].normal - glm::vec4(1.0f, 1.0f, 1.0f, 0.0f);
			mesh[i].normal = mesh[i].normal / mesh[i].neighbor_count;
			mesh[i].normal = mynormalize(mesh[i].normal);
			mesh[i].normal[3] = 1.0;
		}
		
		
		myfile.close();
	}
}

void draw_mesh(){
	ifstream myfile;
	myfile.open("icaro.obj");
	if (myfile.is_open())
	{
		char temp;
		float t1, t2, t3;
		drawMesh.clear();
		myfile >> temp >> t1 >> t2 >> t3;
		while (temp == 'v')
		{
			myfile >> temp >> t1 >> t2 >> t3;
		}
		int i1, i2, i3;
		i1 = (int)t1 -1;
		i2 = (int)t2 -1;
		i3 = (int)t3 -1;
		while (myfile.good())
		{
			//do the stuff
			if (i1<0 || i2<0||i3<0)
			{
				cout<<"''''''''''''''''''''''''''''''''''''''''''''''''"<< temp << " " <<i1<<" " << i2<<" " <<i3<<endl;
			}
			
			drawMesh.push_back(mesh[i1]);
			drawMesh.push_back(mesh[i2]);
			drawMesh.push_back(mesh[i3]);
			myfile >> temp >> i1 >> i2 >> i3;
			i1-=1;
			i2-=1;
			i3-=1;
			
		}
		
		glm::mat4 model = glm::scale(glm::mat4(1.0f), glm::vec3(20.0f));
		for (int i = 0; i < drawMesh.size(); i++)
		{
			drawMesh[i].pos = drawMesh[i].pos * model;			
			
		}
		
		mesh.clear();
		myfile.close();
		//~ cout<<"'-------------------size"<< drawMesh.size()<<" " <<mesh.size()<<endl;
	}
}



float cx=0.1f; float cy=0.8f;
float rotate_x; float rotate_y; float rotate_z;

void set_mvp(){	
	rotate_x=-cos(cx)*sin(cy)*30.0;
    rotate_y=sin(cx)*30.0;
    rotate_z=cos(cx)*cos(cy)*30.0;

	glm::vec3 POS=glm::vec3(rotate_x,rotate_y,rotate_z);
	glm::vec3 AT=glm::vec3(0.0f,0.0f,0.0f);
	glm::vec3 UP=glm::vec3(0.0f,1.0f,0.0f);
	glm::mat4 view_matrix = glm::lookAt(POS,AT,UP); 
		
	projection_matrix = glm::ortho<GLfloat>(-5.0, 5.0, -5.0, 5.0,0.1, 100.0);
	// Model matrix : an identity matrix (model will be at the origin)
	model_matrix      = glm::mat4(1.0f);
	modelview_matrix = view_matrix*model_matrix;
	mvp_matrix = projection_matrix*view_matrix*model_matrix;

	normal_mat = glm::transpose(glm::inverse(glm::mat3(modelview_matrix)));
}


void initSphere( vector<Vertex> * sphere, float cx, float cy, float cz, float Radius, int Resolution ) {
	int w, h;
	
	// Auxiliary variables
	float X1, Y1, X2, Y2, Z1, Z2;
	float inc1, inc2, inc3, inc4, inc5, Radius1, Radius2;
	
	for( w = 0; w < Resolution; w++ ) {
		
		for( h = (-Resolution/2); h < (Resolution/2); h++ ) {
			inc1 = (w/(float)Resolution)*2*PI;
			inc2 = ((w+1)/(float)Resolution)*2*PI;
			inc3 = (h/(float)Resolution)*PI;
			inc4 = ((h+1)/(float)Resolution)*PI;

			X1 = sin(inc1);
			Y1 = cos(inc1);
			X2 = sin(inc2);
			Y2 = cos(inc2);
			
			// store the upper and lower radius
			Radius1 = Radius*cos(inc3);
			Radius2 = Radius*cos(inc4);

			Z1 = Radius*sin(inc3);
			Z2 = Radius*sin(inc4);

			Vertex v1,v2,v3,v4,v5,v6;
			// insert the triangle coordinates
			v1.pos 	  = glm::vec4( Radius1*X1 + cx, Z1 + cy, Radius1*Y1 + cz, 1.0 );
			v1.color  = glm::vec4( 0.0, 0.0, 1.0, 0.99 );
			v1.normal = glm::vec4(X1,Z1,Y1,1.0) / glm::length( glm::vec3(X1,Z1,Y1) );
			
			v2.pos 	  = glm::vec4( Radius1*X2 + cx, Z1 + cy, Radius1*Y2 + cz, 1.0 );
			v2.color  = glm::vec4( 0.0, 0.0, 1.0, 0.99 );
			v2.normal = glm::vec4(X2,Z1,Y2,1.0) / glm::length( glm::vec3(X2,Z1,Y2) ); 
			
			v3.pos 	  = glm::vec4( Radius2*X2 + cx, Z2 + cy, Radius2*Y2 + cz, 1.0 );
			v3.color  = glm::vec4( 0.0, 0.0, 1.0, 0.99 );
			v3.normal = glm::vec4(X2,Z2,Y2,1.0) / glm::length( glm::vec3(X2,Z2,Y2) ); 
			
			v4.pos 	  = glm::vec4( Radius1*X1 + cx, Z1 + cy, Radius1*Y1 + cz, 1.0 );
			v4.color  = glm::vec4( 0.0, 0.0, 1.0, 0.99 );
			v4.normal = glm::vec4(X1,Z1,Y1,1.0) / glm::length( glm::vec3(X1,Z1,Y1) );
			
			v5.pos 	  = glm::vec4( Radius2*X2 + cx, Z2 + cy, Radius2*Y2 + cz, 1.0 );
			v5.color  = glm::vec4( 0.0, 0.0, 1.0, 0.99 );
			v5.normal = glm::vec4(X2,Z2,Y2,1.0) / glm::length( glm::vec3(X2,Z2,Y2) );
			
			v6.pos 	  = glm::vec4( Radius2*X1 + cx, Z2 + cy, Radius2*Y1 + cz, 1.0 );
			v6.color  = glm::vec4( 0.0, 0.0, 1.0, 0.99 );		
			v6.normal = glm::vec4(X1,Z2,Y1,1.0) / glm::length( glm::vec3(X1,Z2,Y1) );
												
			sphere->push_back( v1 );
			sphere->push_back( v2 );
			sphere->push_back( v3 );
			sphere->push_back( v4 );
			sphere->push_back( v5 );
			sphere->push_back( v6 );
		}
	}	
}

void initialize_points(){
    init_pt_array();
    createtop();
    create_toptriangles();
  //  create_toptopt();
    initialize_hair();
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
    
    //Initialization of the points of the eyes
	initSphere( &eyeLeft, -2.2, 3.08, 0.505, 0.1, 10);
	initSphere( &eyeRight, -1.58, 3.12, 1.67, 0.1, 10);
}

// MOUSE handling *******************************************

int last_x, last_y;
int selected_idx = -1;
int selected_idy = -1;

void mouse(int button, int state, int x, int y ){
	if(state == GLUT_DOWN){
		selected_idx = color_code_pick_red(x,y);
		selected_idy = color_code_pick_green(x,y);
		printf("pick point #%d , %d\n", selected_idx, selected_idy);
	}

    last_x = x;
    last_y = y;
}


void motion( int x, int y){
    GLint viewport[4];
    float dx,dy;
    
    glGetIntegerv(GL_VIEWPORT,viewport);
    
    dx = (x - last_x)/(float)viewport[2] * 2;
    dy = -(y - last_y)/(float)viewport[3] * 2;

    last_x = x;
    last_y = y;
    
    if (selected_idx == -1 || selected_idy == -1)
	{
		last_x = x;
		last_y = y;	
	}
	else
	{
		if(selected_idx > 0 && selected_idx < ptArrayLength && selected_idy > 0 && selected_idy < ptArrayLength)
		{
		pt_array[selected_idx][selected_idy].pos[0] += dx;
		pt_array[selected_idx][selected_idy].pos[1] += dy;
		pt_array[selected_idx][selected_idy].pos[2] -= dx;
		
		printf("x, y, z %f %f %f  \n", pt_array[selected_idx][selected_idy].pos[0], pt_array[selected_idx][selected_idy].pos[1], pt_array[selected_idx][selected_idy].pos[2]);
		}
		
		last_x = x;
		last_y = y;	
	}
    glutPostRedisplay();    
}

// KEYBOARD handling *******************************************
void keyboard(unsigned char key, int x, int y)
{
    switch (key) {
		case 'r':
			graphics_init();
			camx = 10.0;
			camy = 10.0;
			camz = 10.0;
			set_mvp();
			glutPostRedisplay();
			break;
		case 27:
			exit(0);
			break;
		case 'L':
			if(!isLighting)
				isLighting = true;
			else
				isLighting = false;
			break;
		case 'l':
			if(isLighting != 0.0)
				isLighting = 0.0;
			else
				isLighting = 0.1;
			break;
		case '1':
			if(isLighting != 0.0)
				isLighting = 0.0;
			else
				isLighting = 0.1;
			break;
		case '2':
			if(!isDisc)
				isDisc = true;
			else
				isDisc = false;
			break;
		case 'd':
			if(showPoints)
				showPoints = false;
			else
				showPoints = true;
			break;
		case 's':
			if(showSurface)
				showSurface = false;
			else
				showSurface = true;
			break;
		case 't':
			if(isTexture == 0.0)
				isTexture = 0.1;
			else
				isTexture = 0.0;
			break;
		case 'f':
			if(showMesh)
				showMesh = false;
			else
				showMesh = true;
			break;
		case 'a':
			if (showAnimation)
				showAnimation = false;
			else
				showAnimation = true;
			break;
		case 'h':
			if (showhair==0)
				showhair = 1;
			else if( showhair==1)
			{
				showhair = 2;
			}
			else if(showhair==2)
			{
			    showhair==0;
			}
			break;	
		case 'e':
			if(showEyes)
				showEyes = false;
			else
				showEyes = true;
			break;
    }
    
    glutPostRedisplay();
}

void special_key(int key, int x, int y)
{
    //Capture the arrow key
	if(key==GLUT_KEY_LEFT) {
		cy=cy+0.05f;
	} else if(key==GLUT_KEY_RIGHT) {
		cy=cy-0.05f;
	} else if(key==GLUT_KEY_UP) {
		cx=cx+0.05f;
	} else if(key==GLUT_KEY_DOWN) {
		cx=cx-0.05f;
	}

    glutPostRedisplay();
}

// DISPLAY and RENDERING functions *************************
void draw_grid(){
    
    create_sub_surface();
    GLuint vao;
    GLuint grid_vbo;
    
    GLuint position_location;
    GLuint color_location;
    GLuint MVP_location;
    GLuint normal_location;
    GLuint light_pos_location;
    GLuint light_intensity_location;
    GLuint texture_location;
    
    GLuint camx_loc;
    GLuint camy_loc;
    GLuint camz_loc;
    GLuint isL_loc;
    
    GLuint isT_loc;
    
    GLuint top_vbo;
    GLuint toptop_vbo;
    GLuint hair_vbo;

    //~ GLenum ErrorValue = glGetError();
    //~ cout << "--------------------- " << ErrorValue << endl;
    // specify the shaders we want to use
    glUseProgram(program);
    //~ GLenum ErrorValue1 = glGetError();
    //~ cout << "---------------------1 " << ErrorValue1 << endl;

    // get the input variable location in this shader
    position_location = glGetAttribLocation(program, "in_position");
    color_location = glGetAttribLocation(program, "in_color");
    normal_location = glGetAttribLocation(program, "in_normal");
    texture_location = glGetAttribLocation(program, "in_texcoord");

    glUniformMatrix4fv(glGetUniformLocation(program, "modelview_mat"), 1, GL_FALSE, &modelview_matrix[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(program, "proj_mat"), 1, GL_FALSE, &projection_matrix[0][0]);
    glUniformMatrix3fv(glGetUniformLocation(program, "NormalMatrix"), 1, GL_FALSE, &normal_mat[0][0]);

	glUniform1i(glGetUniformLocation(program, "face_tex"), 0);

    isL_loc = glGetUniformLocation(program, "isL");
    glUniform1f(isL_loc, isLighting);
        
    isT_loc = glGetUniformLocation(program, "isT");
 
    // create and bind a VAO
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // enable the input locations we want to use
    glEnableVertexAttribArray(normal_location);
    glEnableVertexAttribArray(position_location);
    glEnableVertexAttribArray(color_location);
    glEnableVertexAttribArray(texture_location);
    
    //~ printf("LOCATIONS! %d %d %d    %d %d    %d %d  \n", position_location, color_location, normal_location,    texture_location,  glGetUniformLocation(program, "NormalMatrix"), glGetUniformLocation(program, "proj_mat"), glGetUniformLocation(program, "modelview_mat")); 

    // draw points
    glUniform1f(isT_loc, 0.0);

    // generate and bind a vertex buffer to hold the vertex data on GPU
    glGenBuffers(1, &grid_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, grid_vbo);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // initialize the vertex buffer with the vertex data
    glBufferData(GL_ARRAY_BUFFER, grid_vertices.size() * sizeof(Vertex), &grid_vertices[0] , GL_STATIC_DRAW);

    // specify the input data format
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
	glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
	glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));

    // draw points
    glPointSize(10);
    glDrawArrays(GL_LINES, 0, grid_vertices.size());
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    /**************Point array********/
    if (showPoints)
	{
		// draw points

		// generate and bind a vertex buffer to hold the vertex data on GPU
		glGenBuffers(1, &grid_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, grid_vbo);

		// initialize the vertex buffer with the vertex data
		glBufferData(GL_ARRAY_BUFFER, ptArrayLength * ptArrayLength * sizeof(Vertex), &pt_array[0] , GL_STATIC_DRAW);

		// specify the input data format
		glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
		glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));

		// draw points
		glPointSize(5);
		glDrawArrays(GL_POINTS, 0, ptArrayLength*ptArrayLength);
		
		
		glGenBuffers(1, &grid_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, grid_vbo);

		// initialize the vertex buffer with the vertex data
		glBufferData(GL_ARRAY_BUFFER, ptArrayLength * ptArrayLength * sizeof(Vertex), &top[0] , GL_STATIC_DRAW);

		// specify the input data format
		glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
		glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));

		// draw points
		glPointSize(5);
		glDrawArrays(GL_POINTS, 0, ptArrayLength*ptArrayLength);
		
	}
	

	if(isTexture)	
{	//draw top_triangles.................................
		glGenBuffers(1, &top_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, top_vbo);

		// initialize the vertex buffer with the vertex data
		
		glBufferData(GL_ARRAY_BUFFER, 3*2*360*19* sizeof(Vertex), &top_triangles[0] , GL_STATIC_DRAW);

		//~ printf("%d  %d  %d \n", position_location,color_location,normal_location);

		// specify the input data format
		glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
		glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));


		// draw points
		//glPointSize(10);
		glDrawArrays(GL_TRIANGLES, 0, 3*2*360*19);
		
	//draw toptop............................
    		glGenBuffers(1, &toptop_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, toptop_vbo);

		// initialize the vertex buffer with the vertex data
		
		glBufferData(GL_ARRAY_BUFFER, 3*360* sizeof(Vertex), &toptop[0] , GL_STATIC_DRAW);

		//~ printf("%d  %d  %d \n", position_location,color_location,normal_location);

		// specify the input data format
		glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
		glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));


		// draw points
		//glPointSize(10);
		glDrawArrays(GL_TRIANGLES, 0, 3*360);
    }
	if(showhair==1){
       //draw hair.....................
        glGenBuffers(1, &hair_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, hair_vbo);

		// initialize the vertex buffer with the vertex data
		
		glBufferData(GL_ARRAY_BUFFER, 3*20*360* sizeof(Vertex), &hair[0] , GL_STATIC_DRAW);

		//~ printf("%d  %d  %d \n", position_location,color_location,normal_location);

		// specify the input data format
		glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
		glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));

		// draw points
		//glPointSize(10);
		glDrawArrays(GL_TRIANGLES, 0, 20*3*360);
     }   
	 else 	if(showhair==2){
	   initialize_rehair();
       //draw hair.....................
        glGenBuffers(1, &hair_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, hair_vbo);

		// initialize the vertex buffer with the vertex data
		
		glBufferData(GL_ARRAY_BUFFER, 3*20*360* sizeof(Vertex), &hair[0] , GL_STATIC_DRAW);

		//~ printf("%d  %d  %d \n", position_location,color_location,normal_location);

		// specify the input data format
		glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
		glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));

		// draw points
		//glPointSize(10);
		glDrawArrays(GL_TRIANGLES, 0, 20*3*360);
     }   
    /**************Point array********/
	if (showMesh)
	{
		// draw points
		
		read_mesh();
		draw_mesh();
	
		glUniform1f(isT_loc, 0.0);
		// generate and bind a vertex buffer to hold the vertex data on GPU
		glGenBuffers(1, &grid_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, grid_vbo);

		// initialize the vertex buffer with the vertex data
		
		glBufferData(GL_ARRAY_BUFFER, drawMesh.size() * sizeof(Vertex), &drawMesh[0] , GL_STATIC_DRAW);

		//~ printf("%d  %d  %d \n", position_location,color_location,normal_location);

		// specify the input data format
		glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
		glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));

		// draw points
		glPointSize(10);
		glDrawArrays(GL_TRIANGLES, 0, drawMesh.size());		
	}
	
    /**************subsurface********/
    if (showSurface)
	{
		glUniform1f(isT_loc, isTexture);
		
		create_sub_surface();
		create_rects();

		// generate and bind a vertex buffer to hold the vertex data on GPU
		glGenBuffers(1, &grid_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, grid_vbo);

		// initialize the vertex buffer with the vertex data
		
		glBufferData(GL_ARRAY_BUFFER, temp.size() * sizeof(Vertex), &temp[0] , GL_STATIC_DRAW);

		//~ printf("%d  %d  %d \n", position_location,color_location,normal_location);

		// specify the input data format
		glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
		glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		glVertexAttribPointer(texture_location, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coord));

		// draw points
		glPointSize(5);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, temp.size());		
	}
	
	/*********************************************** sphere ********************************************/
	if (showEyes) {
		/********* Left Eye *******/
		// initialize the vertex buffer with the vertex data
		glBufferData(GL_ARRAY_BUFFER, eyeLeft.size() * sizeof(Vertex), &eyeLeft[0] , GL_STATIC_DRAW);
										
		// specify the input data format
		glVertexAttribPointer( position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer( color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*) offsetof(Vertex, color) );
		
		// draw triangles
		glDrawArrays( GL_TRIANGLES, 0, eyeLeft.size() - 2 );
		
		/********* Right Eye *******/
		// initialize the vertex buffer with the vertex data
		glBufferData(GL_ARRAY_BUFFER, eyeRight.size() * sizeof(Vertex), &eyeRight[0] , GL_STATIC_DRAW);
										
		// specify the input data format
		glVertexAttribPointer( position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
		glVertexAttribPointer( color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*) offsetof(Vertex, color) );
		
		// draw triangles
		glDrawArrays( GL_TRIANGLES, 0, eyeRight.size() - 2 );
	}
	
    glUniform1f(isT_loc, 0.0);
 
    // unbind VAO and VBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Delete VAO and VBO
    glDeleteBuffers(1, &top_vbo);
    glDeleteBuffers(1, &grid_vbo);
    glDeleteVertexArrays(1, &vao);
}


void display(){
    // Clear Viewport
    glClearColor(0.0f,0.0f,0.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    for (int i = 0; i < ptArrayLength; i++)
	{
		for (int j = 0; j < ptArrayLength; j++)
		{
			pt_array[i][j].color = glm::vec4(0.0, 1.0, 0.0, 1.0);
		}		
	}
	
	set_mvp();
    draw_grid();

    glFlush();
    glutSwapBuffers();

	for (int i = 0; i < ptArrayLength; i++)
	{
		for (int j = 0; j < ptArrayLength; j++)
		{
			pt_array[i][j].color = glm::vec4(((float)i/(float)ptArrayLength)* 0.8 + 0.2, ((float)j/(float)ptArrayLength)* 0.8 + 0.2, 0.1, 1.0);
		}		
	}
}

void timer(int value)
{
    // do the update you want, e.g. update the frame position
    if (showAnimation)
	{
	
		if (ind < 4)
		{
			pt_array[28][7].pos = animation[0][ind];
			pt_array[27][12].pos = animation[1][ind];
			pt_array[28][12].pos = animation[2][ind];
			pt_array[26][12].pos = animation[3][ind];
			pt_array[28][13].pos = animation[4][ind];
			pt_array[29][7].pos = animation[5][ind];
			pt_array[31][14].pos = animation[6][ind];
			pt_array[29][6].pos = animation[7][ind];
			pt_array[28][6].pos = animation[8][ind];
			pt_array[30][7].pos = animation[9][ind];
			pt_array[27][7].pos = animation[10][ind];
			
			create_sub_surface();
			create_rects();
			glutPostRedisplay();
			ind++;
		}
		
	}

    // reset timer to be trigger every 60 ms
    glutTimerFunc(60, &timer, 0);

    // redraw the screen
    glutPostRedisplay();
}

void reshape(int width, int height){
    // Clip the view port to match our ratio
    glViewport(0, 0, width, height);
    glutPostRedisplay();
}

void graphics_init(){

    // init vertex shader
    read_shader_source_code("vs.glsl", buf, 2048);
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

	GLuint texture_id = load_texture_TGA("icaro.tga", NULL, NULL, GL_CLAMP, GL_CLAMP );
	
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
    set_mvp();

    // set up callbacks
    glutReshapeFunc(&reshape);
    glutDisplayFunc(&display);
    glutKeyboardFunc(&keyboard);
    glutMouseFunc(&mouse);
    glutMotionFunc(&motion);
    glutSpecialFunc(&special_key);
	glutTimerFunc(60, &timer, 0);

    // main loop
    glutMainLoop();
    return 0;
}
