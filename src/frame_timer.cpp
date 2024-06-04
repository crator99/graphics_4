/*
   How to use this code:

   Call init_timer before starting rendering, i.e., before calling
   glutMainLoop. Then, make sure your display function is organized
   roughly as the example below.
*/

#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <float.h>
#include <iostream>
#include <math.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "load_mesh.cpp"

float                 gTotalTimeElapsed    = 0;
int                gTotalFrames      = 0;
GLuint                gTimer;

GLuint VAO, VBO[2], EBO;

void init_timer()
{
   glGenQueries(1, &gTimer);
}

void start_timing()
{
   glBeginQuery(GL_TIME_ELAPSED, gTimer);
}

float stop_timing()
{
   glEndQuery(GL_TIME_ELAPSED);

   GLint available = GL_FALSE;
   while (available == GL_FALSE)
      glGetQueryObjectiv(gTimer, GL_QUERY_RESULT_AVAILABLE, &available);

   GLint result;
   glGetQueryObjectiv(gTimer, GL_QUERY_RESULT, &result);

   float timeElapsed = result / (1000.0f * 1000.0f * 1000.0f);
   return timeElapsed;
}

/*
   Your display function should look roughly like the following.
*/

void display_immediate_mode() {
    // 화면 및 깊이 버퍼 초기화
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    start_timing();

    glBegin(GL_TRIANGLES);
    for (const auto& triangle : gTriangles) {
        for (int i = 0; i < 3; ++i) {
            int index = triangle.indices[i];
            const auto& normal = gNormals[index];
            const auto& position = gPositions[index];
            glNormal3f(normal.x, normal.y, normal.z);
            glVertex3f(position.x, position.y, position.z);
        }
    }
    glEnd();

    float timeElapsed = stop_timing();
    gTotalFrames++;
    gTotalTimeElapsed += timeElapsed;
    float fps = gTotalFrames / gTotalTimeElapsed;
    std::string title = "Immediate Mode: " + std::to_string(fps) + " FPS";
}

void display_vertex_array_mode() {
    // 화면 및 깊이 버퍼 초기화
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    start_timing();

    GLuint positionBuffer;
    GLuint normalBuffer;
    GLuint indexBuffer;

    glGenBuffers(1, &positionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
    glBufferData(GL_ARRAY_BUFFER, gPositions.size() * sizeof(Vector3), gPositions.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &normalBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
    glBufferData(GL_ARRAY_BUFFER, gNormals.size() * sizeof(Vector3), gNormals.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &indexBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, gTriangles.size() * sizeof(Triangle), gTriangles.data(), GL_STATIC_DRAW);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
    glVertexPointer(3, GL_FLOAT, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
    glNormalPointer(GL_FLOAT, 0, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glDrawElements(GL_TRIANGLES, gTriangles.size() * 3, GL_UNSIGNED_INT, 0);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);

    float timeElapsed = stop_timing();
    gTotalFrames++;
    gTotalTimeElapsed += timeElapsed;
    float fps = gTotalFrames / gTotalTimeElapsed;
    std::string title = "Vertex Array Mode: " + std::to_string(fps) + " FPS";
}