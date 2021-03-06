#ifdef _WIN32
#include <windows.h>
#endif

#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glu.h>
#include "\Code\vcpkg\Microsoft\vcpkg\installed\x64-windows\include\GL\glut.h"
#else
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#endif

#include "RayTracer.h"
#include <iostream>

//--------------------------------------------------------------------------------
// Change these to test ray tracer at different resolutions
int WIDTH = 640;
int HEIGHT = 480;
//--------------------------------------------------------------------------------

Render_World world;
Driver driver(world);

void Display()
{
    //glClear(GL_COLOR_BUFFER_BIT);
    //previously using gl_unsigned_int_8_8_8_8 
    if(world.camera.film.colors) glDrawPixels(WIDTH,HEIGHT,GL_RGBA, GL_UNSIGNED_BYTE,(GLvoid*)world.camera.film.colors);
    glFlush();
}

void Handle_Idle()
{
    driver.Render_More();
    Display();
}

void Initialize_Opengl_And_Glut(int argc, char** argv)
{
    glutInit(&argc,(char**)argv);
    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGBA);
    glutInitWindowSize(WIDTH,HEIGHT);
    glutCreateWindow("Ray Tracer");
    glutDisplayFunc(Display);
    glutIdleFunc(Handle_Idle);
    glClearColor(0,0,0,0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, WIDTH, 0.0, HEIGHT, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void Usage(const std::string& exec)
{
    std::cout<<"Usage: "<<exec<<" <test number>"<<"\n\t <test number> 1-4"<<std::endl;
    exit(1);
}

//--------------------------------------------------------------------------------
// Initialize_World
//    Sets up objects, materials, and lights for tests 1 - 4
//--------------------------------------------------------------------------------
void Initialize_World(Render_World& world,const int width,const int height,const int test_number)
{
    Vector_3D<double> color1(1,0,0);
    Vector_3D<double> color2(0.2,0.2,0.8);
    Vector_3D<double> plane_color(1,0.775,0.5431);

    Plane* plane=new Plane(Vector_3D<double>(),Vector_3D<double>(0,1,0));
    Sphere* sphere=new Sphere(Vector_3D<double>(0,1,0),1.0);
    Sphere* sphere2=new Sphere(Vector_3D<double>(-2,1,-1),1.0);


    if (test_number == 1)
    {
        plane->material_shader = new Flat_Shader(world, plane_color);
        sphere->material_shader = new Flat_Shader(world, color1);
        sphere2->material_shader = new Flat_Shader(world, color2);
        world.enable_shadows = false;
        world.recursion_depth_limit = 0;
        world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6), Vector_3D<double>(0, 1, 0), Vector_3D<double>(0, 1, 0));
        world.camera.Focus_Camera(1, (double)width / (double)height, (double)70 / (double)180 * PI);
    }
    else if (test_number == 2)
    {
        plane->material_shader = new Phong_Shader(world, plane_color, plane_color);
        sphere->material_shader = new Phong_Shader(world, color1, color1);
        sphere2->material_shader = new Phong_Shader(world, color2, color2);
        world.enable_shadows = false;
        world.recursion_depth_limit = 0;
        world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6), Vector_3D<double>(0, 1, 0), Vector_3D<double>(0, 1, 0));
        world.camera.Focus_Camera(1, (double)width / (double)height, (double)70 / (double)180 * PI);
    }
    else if (test_number == 3)
    {
        plane->material_shader = new Phong_Shader(world, plane_color, plane_color);
        sphere->material_shader = new Phong_Shader(world, color1, color1);
        sphere2->material_shader = new Phong_Shader(world, color2, color2);
        world.enable_shadows = true;
        world.recursion_depth_limit = 0;
        world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6), Vector_3D<double>(0, 1, 0), Vector_3D<double>(0, 1, 0));
        world.camera.Focus_Camera(1, (double)width / (double)height, (double)70 / (double)180 * PI);
    }
    else if (test_number == 4) 
    {
        plane->material_shader = new Phong_Shader(world, plane_color, plane_color);
        sphere->material_shader = new Reflective_Shader(world, color1, color1, Vector_3D<double>(1, 1, 1), 50, .5);
        sphere2->material_shader = new Reflective_Shader(world, color2, color2, Vector_3D<double>(1, 1, 1), 50, .5);
        world.enable_shadows = true;
        world.recursion_depth_limit = 4;
        world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6), Vector_3D<double>(-1, 1, 0), Vector_3D<double>(0, 1, 0));
        world.camera.Focus_Camera(1, (double)width / (double)height, (double)50 / (double)180 * PI);
    }
    else
    {
        std::cout << "Unrecognized test number" << std::endl;
        exit(1);
    }

    // Lights
    Light* point_light = new Point_Light(Vector_3D<double>(-2,7,-3),Vector_3D<double>(1,1,1),0.5);
    Light* point_light2 = new Point_Light(Vector_3D<double>(3,5,-3),Vector_3D<double>(1,1,1),0.5);

    world.objects.push_back(plane);
    world.objects.push_back(sphere);
    world.objects.push_back(sphere2);
    world.lights.push_back(point_light);
    world.lights.push_back(point_light2);

    // Camera
    world.camera.film.Set_Resolution(width,height);
}
//--------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    //VS being dumb, dont know how to add cmd args
    //if(argc!=2) Usage(argv[0]);
    //int test_number=atoi(argv[1]);
    //if(test_number<1 || test_number>4) Usage(argv[0]);
    int test_number = 4;
    Initialize_Opengl_And_Glut(argc,argv);
    Initialize_World(world,WIDTH,HEIGHT,test_number);

    glutMainLoop();

    return 0;
}