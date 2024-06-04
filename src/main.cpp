//
//  sphere_scene.c
//  Rasterizer
//
//
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "frame_timer.cpp"

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void setup() {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.1, -1.0, -1.5);
    glScalef(10, 10, 10);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-0.1, 0.1, -0.1, 0.1, 0.1, 1000);

    glViewport(0, 0, 512, 512);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    float ka[] = { 1, 1, 1, 1 };
    float kd[] = { 1, 1, 1, 1 };
    float ks[] = { 0, 0, 0, 1 };
    float p = 0;

    glMaterialfv(GL_FRONT, GL_AMBIENT, ka);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, kd);
    glMaterialfv(GL_FRONT, GL_SPECULAR, ks);
    glMaterialf(GL_FRONT, GL_SHININESS, p);

    float Ia[] = {0.2, 0.2, 0.2, 0 };
    float l[] = {1, 1, 1, 0 };
    float la[] = {0, 0, 0, 0 };
    float ld[] = {10, 10, 10, 0 };
    float ls[] = {0, 0, 0, 0 };

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Ia);
    glLightfv(GL_LIGHT0, GL_POSITION, l);
    glLightfv(GL_LIGHT0, GL_AMBIENT, la);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, ld);
    glLightfv(GL_LIGHT0, GL_SPECULAR, ls);
}

void setWindowTitle(GLFWwindow* window, const std::string& title) {
    glfwSetWindowTitle(window, title.c_str());
}

int main(int argc, char* argv[])
{
    glfwInit();
        
    GLFWwindow* window1 = glfwCreateWindow(512, 512, "Intermiate Mode", NULL, NULL);
    if (window1 == NULL) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window1);
    glfwSetFramebufferSizeCallback(window1, framebuffer_size_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    GLFWwindow* window2 = glfwCreateWindow(512, 512, "Vertex Array", NULL, NULL);
    if (window2 == NULL) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwSetFramebufferSizeCallback(window2, framebuffer_size_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD for 2nd window" << std::endl;
        return -1;
    }

    glEnable(GL_DEPTH_TEST);
     
    load_mesh("bunny.obj");
    setup();

    init_timer();
  
    while (!glfwWindowShouldClose(window1) && !glfwWindowShouldClose(window2)) {
        glfwMakeContextCurrent(window1);
        display_immediate_mode();
        //immediate mode 디스플레이 호출
        
        float fps1 = gTotalFrames / gTotalTimeElapsed;
        setWindowTitle(window1, "Immediate Mode: " + std::to_string(fps1) + " FPS");
        glfwSwapBuffers(window1);
      
        glfwPollEvents();
      
        glfwMakeContextCurrent(window2);
        glEnable(GL_DEPTH_TEST);
        
        setup();
        display_vertex_array_mode();
        //vertex array 디스플레이 호출
        
        float fps2 = gTotalFrames / gTotalTimeElapsed;
        setWindowTitle(window2, "Vertex Array Mode: " + std::to_string(fps2) + " FPS");
        glfwSwapBuffers(window2);

        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}