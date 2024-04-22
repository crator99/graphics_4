//
//  sphere_scene.c
//  Rasterizer
//
//
#include "color.h"
#include "ray.h"
#include "vec3.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <glad/glad.h>
#include <glad/glut.h>
#include <GLFW/glfw3.h>
unsigned char* Image = new unsigned char [512 * 512 * 3] ();
#define M_PI 3.14159265358979323846

int     gNumVertices    = 0;    // Number of 3D vertices.
int     gNumTriangles   = 0;    // Number of triangles.
int*    gIndexBuffer    = NULL; // Vertex indices for the triangles.

vec3* vertices = NULL;

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void MatrixMutiply(float result[4][4], const float mat1[4][4], const float mat2[4][4]){
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            result[i][j] = 0.0;
            for (int k = 0; k < 4; k++){
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

void MatrixMultiply1(float result[4][1], const float mat1[4][4], const float mat2[4][1]){
    for (int i = 0; i < 4; i++){
        result[i][0] = 0.0;
        for (int j = 0; j < 4; j++){
            result[i][0] += mat1[i][j] * mat2[j][0];
        }
    }
}

void MatrixNorm(float result[4][1]){
    float w = result[3][0];
    for(int i = 0; i <4; i++){
        result[i][0] /= w;
    }
}

vec3* create_scene(float Mvp[][4], float Morth[][4], float P[][4], float Mcam[][4], float Mm[][4])
{
    int width   = 32;
    int height  = 16;
    
    float theta, phi;
    int t;
    float result4x4[4][4];
    float result4x5[4][4];
    float result4x6[4][4];
    float result4x7[4][4];
    float result4x1[4][1];

    MatrixMutiply(result4x4, Mvp, Morth);
    MatrixMutiply(result4x5, result4x4, P);
    MatrixMutiply(result4x6, result4x5, Mcam);
    MatrixMutiply(result4x7, result4x6, Mm);

    gNumVertices    = (height - 2) * width + 2;
    gNumTriangles   = (height - 2) * (width - 1) * 2;

    // TODO: Allocate an array for gNumVertices vertices.
    vertices = new vec3[gNumVertices];

    gIndexBuffer    = new int[3*gNumTriangles];
    
    t = 0;
    for (int j = 1; j < height-1; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            theta = (float) j / (height-1) * M_PI;
            phi   = (float) i / (width-1)  * M_PI * 2;
            
            float   x   = sinf(theta) * cosf(phi);
            float   y   = cosf(theta);
            float   z   = -sinf(theta) * sinf(phi);

            float p0[4][1] = {{x}, {y}, {z}, {1}};

            MatrixMultiply1(result4x1, result4x7, p0);
            MatrixNorm(result4x1);
            
            // TODO: Set vertex t in the vertex array to {x, y, z}.
            vertices[t] = vec3(result4x1[0][0], result4x1[1][0], result4x1[2][0]);
            
            t++;
        }
    }
    
    // TODO: Set vertex t in the vertex array to {0, 1, 0}.
    // vertices[t] = float(0.0, 1.0, 0.0);
    float p0[4][1] = {{0.0}, {1.0}, {0.0}, {1.0}};

            // float p_vec[4][1] = 
    MatrixMultiply1(result4x1, result4x7, p0);
    MatrixNorm(result4x1);

            // std::cout << p_vec[1][4] << std::endl;
            
            // TODO: Set vertex t in the vertex array to {x, y, z}.
    vertices[t] = vec3(result4x1[0][0], result4x1[1][0], result4x1[2][0]);
            
    t++;
    
    // TODO: Set vertex t in the vertex array to {0, -1, 0}.
    // vertices[t] = vec3(0.0, -1.0, 0.0);
    float p1[4][1] = {{0.0}, {-1.0}, {0.0}, {1.0}};

            // float p_vec[4][1] = 
    MatrixMultiply1(result4x1, result4x7, p1);
    MatrixNorm(result4x1);

            // std::cout << p_vec[1][4] << std::endl;
            
            // TODO: Set vertex t in the vertex array to {x, y, z}.
    vertices[t] = vec3(result4x1[0][0], result4x1[1][0], result4x1[2][0]);
            
    t++;
    
    t = 0;
    for (int j = 0; j < height-3; ++j)
    {
        for (int i = 0; i < width-1; ++i)
        {
            gIndexBuffer[t++] = j*width + i;
            gIndexBuffer[t++] = (j+1)*width + (i+1);
            gIndexBuffer[t++] = j*width + (i+1);
            gIndexBuffer[t++] = j*width + i;
            gIndexBuffer[t++] = (j+1)*width + i;
            gIndexBuffer[t++] = (j+1)*width + (i+1);
        }
    }
    for (int i = 0; i < width-1; ++i)
    {
        gIndexBuffer[t++] = (height-2)*width;
        gIndexBuffer[t++] = i;
        gIndexBuffer[t++] = i + 1;
        gIndexBuffer[t++] = (height-2)*width + 1;
        gIndexBuffer[t++] = (height-3)*width + (i+1);
        gIndexBuffer[t++] = (height-3)*width + i;
    }
    
    // The index buffer has now been generated. Here's how to use to determine
    // the vertices of a triangle. Suppose you want to determine the vertices
    // of triangle i, with 0 <= i < gNumTriangles. Define:
    //
    // k0 = gIndexBuffer[3*i + 0]
    // k1 = gIndexBuffer[3*i + 1]
    // k2 = gIndexBuffer[3*i + 2]
    //
    // Now, the vertices of triangle i are at positions k0, k1, and k2 (in that
    // order) in the vertex array (which you should allocate yourself at line
    // 27).
    //
    // Note that this assumes 0-based indexing of arrays (as used in C/C++,
    // Java, etc.) If your language uses 1-based indexing, you will have to
    // add 1 to k0, k1, and k2.

    return vertices;
}

int main(int argc, char* argv[])
{
    float n_x = 512;
    float n_y = 512;
    float l = -0.1;
    float r = 0.1;
    float b = -0.1;
    float t = 0.1;
    float n = -0.1;
    float f = -1000;
    // float Mvp[4][4]= {(n_x/2, 0, 0, (n_x - 1)/2), (0, n_y/2, 0, (n_y-1)/2), (0, 0, 1, 0), (0, 0, 0, 1)};
    // std::cout << 1144 << std::endl;

    float Mvp[4][4]{
        {n_x/2, 0, 0, (n_x - 1)/2},
        {0, n_y/2, 0, (n_y-1)/2},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    float Morth[4][4] = {
        {2/(r-l), 0, 0, -(r+l)/(r-l)},
        {0, 2/(t-b), 0, -(t+b)/(t-b)},
        {0, 0, 2/(n-f), -(n+f)/(n-f)},
        {0, 0, 0, 1}
    };
    float P[4][4] = {
        {n, 0, 0, 0},
        {0, n, 0, 0},
        {0, 0, (n+f), -(n*f)},
        {0, 0, 1, 0}
    };
    float Mcam[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    float Mm[4][4] = {
        {2, 0, 0, 0},
        {0, 2, 0, 0},
        {0, 0, 2, -7},
        {0, 0, 0, 1}    
    };

    vec3 tri1;
    vec3 tri2;
    vec3 tri3;
    // float p0[1][4] = {(x), (y), (z), (1)};

    create_scene(Mvp, Morth, P, Mcam, Mm);

    // create_scene();

    // MatrixMultiply1(result4x1, result4x7, vertices);

    for (int i = 0; i < gNumTriangles; i++){
        int k0 = gIndexBuffer[3*i + 0];
        int k1 = gIndexBuffer[3*i + 1];
        int k2 = gIndexBuffer[3*i + 2];

        tri1 = vertices[k0];
        tri2 = vertices[k1];
        tri3 = vertices[k2];

        int x_min = floor(std::min(std::min(tri1.x(), tri2.x()), tri3.x()));
        int x_max = ceil(std::max(std::max(tri1.x(), tri2.x()), tri3.x()));
        int y_min = floor(std::min(std::min(tri1.y(), tri2.y()), tri3.y()));
        int y_max = ceil(std::max(std::max(tri1.y(), tri2.y()), tri3.y()));

        float beta = ((tri1.y() - tri3.y()) * x_min + (tri3.x() - tri1.x()) * y_min + tri1.x() * tri3.y() - tri3.x() * tri1.y()) / 
                    ((tri1.y() - tri3.y()) * tri2.x() + (tri3.x() - tri1.x()) * tri2.y() + tri1.x() * tri3.y() - tri3.x() * tri1.y());
        float gamma = ((tri1.y() - tri2.y()) * x_min + (tri2.x() - tri1.x()) * y_min + tri1.x() * tri2.y() - tri2.x() * tri1.y()) /
                    ((tri1.y() - tri2.y()) * tri3.x() + (tri2.x() - tri1.x()) * tri3.y() + tri1.x() * tri2.y() - tri2.x() * tri1.y());

        int n = (x_max - x_min) + 1;

        float beta_x = (tri1.y() - tri3.y()) /
                    ((tri1.y() - tri3.y()) * tri2.x() + (tri3.x() - tri1.x()) * tri2.y() + tri1.x() * tri3.y() - tri3.x() * tri1.y());
        float beta_y = (tri3.x() - tri1.x()) /
                    ((tri1.y() - tri3.y()) * tri2.x() + (tri3.x() - tri1.x()) * tri2.y() + tri1.x() * tri3.y() - tri3.x() * tri1.y());
        float gamma_x = (tri1.y() - tri2.y()) /
                    ((tri1.y() - tri2.y()) * tri3.x() + (tri2.x() - tri1.x()) * tri3.y() + tri1.x() * tri2.y() - tri2.x() * tri1.y());
        float gamma_y = (tri2.x() - tri1.x()) /
                    ((tri1.y() - tri2.y()) * tri3.x() + (tri2.x() - tri1.x()) * tri3.y() + tri1.x() * tri2.y() - tri2.x() * tri1.y());

        for (int y = y_min; y <= y_max; y++){
            for (int x = x_min; x <= x_max; x++){
                if (beta > 0 && gamma > 0 && (beta + gamma) < 1){
                    int index = (y * n_x + x) * 3;
                    Image[index] = 255;
                    Image[index + 1] = 255;
                    Image[index + 2] = 255;
                }
                beta += beta_x;
                gamma += gamma_x;
            }
            beta += beta_y - n * beta_x;
            gamma += gamma_y - n * gamma_x;
        }
    }

    glfwInit();

    if (!glfwInit())
    {
        std::cout << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(n_x, n_y, "HW4_Q1", NULL, NULL);
    
    if(!window){
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glViewport(0, 0, n_x, n_y);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    while(!glfwWindowShouldClose(window))
    {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        glad_glDrawPixels(n_x, n_y, GL_RGB, GL_UNSIGNED_BYTE, Image);
    
        glfwSwapBuffers(window);
        glfwPollEvents();    
    };

    glfwTerminate();
    return 0;
}
