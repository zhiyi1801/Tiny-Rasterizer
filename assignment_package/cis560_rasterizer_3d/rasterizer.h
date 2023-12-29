#pragma once
#include "polygon.h"
#include "camera.h"
#include <QImage>

enum lightMode
{
    Plain,
    Lambertian
};

class Rasterizer
{
private:
    //This is the set of Polygons loaded from a JSON scene file
    std::vector<Polygon> m_polygons;
    Camera myCamera;
public:
    Rasterizer(const std::vector<Polygon>& polygons, Camera = Camera(glm::vec4(0,0,10,1), glm::vec4(0,0,-1,0)));
    QImage RenderScene();
    void ClearScene();
    void moveCamera(glm::vec4 dPos);
    void addCameraFOV(float dFOV);
    void addRotation(glm::vec4 axis, float degree);
    std::vector<float> depthBuffer;
    std::vector<glm::vec3> AABuffer;
    static int AntiAliasing;
    static bool LineMode;
    static void setAntiAliasing(int AA);
    lightMode lMode;
};

void getBoundingBox2D(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3, BoundingBox2D &box, glm::vec2 xlimit, glm::vec2 ylimit);
void getBoundingBox2D(const glm::vec4& p1, const glm::vec4& p2, const glm::vec4& p3, BoundingBox2D &box, glm::vec2 xlimit, glm::vec2 ylimit);

// for p1,p2 (x,y,z in the camera space, z in the screen space)
void Bresenham(const glm::vec4 &p1, const glm::vec4 &p2, std::vector<float> &depthBuffer, std::vector<glm::vec3> &AABuffer, const glm::vec3 &color, const int width, const int height);

void drawTriangle(const lightMode lMode, const Camera &myCam, const std::array<glm::vec4,3> normals, const std::array<Vertex,3>v, const std::array<glm::vec4,3> &screenp, const std::array<glm::vec4, 3> &camerap, QImage &result, std::vector<glm::vec3> &AABuffer, std::vector<float> &depthBuffer, const QImage *currentTexture);

void drawTriangle(const lightMode lMode, const Camera &myCam, const std::array<glm::vec4,3> normals, const std::array<glm::vec2,3>uvMap, const std::array<glm::vec4,3> &screenp, const std::array<glm::vec4, 3> &camerap, QImage &result, std::vector<glm::vec3> &AABuffer, std::vector<float> &depthBuffer, const QImage *currentTexture);

float interpolateZ(glm::vec3 p1, glm::vec3 p2, float targetZ);

float interpolateZ(glm::vec4 p1, glm::vec4 p2, float targetZ);
