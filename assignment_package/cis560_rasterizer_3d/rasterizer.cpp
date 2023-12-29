#include "rasterizer.h"
#include "macro.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <QElapsedTimer>
#include <sys/time.h>

int Rasterizer::AntiAliasing = 1;
bool Rasterizer::LineMode = false;
Rasterizer::Rasterizer(const std::vector<Polygon>& polygons, Camera cam)
    : m_polygons(polygons), myCamera(cam), depthBuffer(WIDTH*HEIGHT,10), AABuffer(WIDTH*HEIGHT),lMode(Plain)
{}

QImage Rasterizer::RenderScene()
{
    const int AA = LineMode ? 1 : this->AntiAliasing;
    const int AA2 = AA*AA;
    QImage result(AA*WIDTH, AA*HEIGHT, QImage::Format_RGB32);
    depthBuffer.resize(result.width()*result.height());
    for(int i = 0; i < (int)depthBuffer.size();++i) depthBuffer[i] = 10;
    AABuffer.resize(result.width()*result.height());
    for(int i = 0; i < (int)AABuffer.size();++i) AABuffer[i] = glm::vec3(0);
    // Fill the image with black pixels.
    // Note that qRgb creates a QColor,
    // and takes in values [0, 255] rather than [0, 1].
    result.fill(qRgb(0.f, 0.f, 0.f));
    // TODO: Complete the various components of code that make up this function.
    // It should return the rasterized image of the current scene.

    // Make liberal use of helper functions; writing your rasterizer as one
    // long RenderScene function will make it (a) hard to debug and
    // (b) hard to write without copy-pasting. Also, Adam will be sad when
    // he reads your code.

    // Also! As per the style requirements for this assignment, make sure you
    // use std::arrays to store things like your line segments, Triangles, and
    // vertex coordinates. This lets you easily use loops to perform operations
    // on your scene components, rather than copy-pasting operations three times
    // each!
#if MODE2D
    depthBuffer.fill(-1.0f);
    BoundingBox2D boundingBox;
    glm::vec3 barycentricv;
    Camera myCamera;
    myCamera.getViewMat();
    for(const Polygon& currentPolygon:m_polygons)
    {
        for (const Triangle &t : currentPolygon.m_tris)
        {
            Vertex v1 = currentPolygon.m_verts[t.m_indices[0]];
            Vertex v2 = currentPolygon.m_verts[t.m_indices[1]];
            Vertex v3 = currentPolygon.m_verts[t.m_indices[2]];

            glm::vec4 p1 = v1.m_pos;
            glm::vec4 p2 = v2.m_pos;
            glm::vec4 p3 = v3.m_pos;

            glm::vec3 color1 = v1.m_color;
            glm::vec3 color2 = v2.m_color;
            glm::vec3 color3 = v3.m_color;

            std::vector<LineSegment2D> l;
            l.push_back(LineSegment2D(p1,p2));
            l.push_back(LineSegment2D(p2,p3));
            l.push_back(LineSegment2D(p3,p1));

            getBoundingBox2D(p1,p2,p3,boundingBox);
            if(boundingBox.xmax <= boundingBox.xmin || boundingBox.ymax <= boundingBox.ymin) continue;
            for(unsigned int yi = std::floor(boundingBox.ymin); yi < std::ceil(boundingBox.ymax); ++yi)
            {
                float hitxmin = WIDTH + 1, hitxmax = -1;
                float hitx[3];
                bool ifHit[3];
                ifHit[0] = l[0].getIntersection(yi,hitx);
                ifHit[1] = l[1].getIntersection(yi,hitx+1);
                ifHit[2] = l[2].getIntersection(yi,hitx+2);
                for(int j = 0; j < 3; ++j)
                {
                    if(ifHit[j])
                    {
                        hitx[j] = std::clamp(hitx[j], 0.0f, float(WIDTH));
                        hitxmin = std::min(hitxmin, hitx[j]);
                        hitxmax = std::max(hitxmax, hitx[j]);
                    }
                }
                for(int xi = std::round(hitxmin); xi <= std::round(hitxmax); ++xi)
                {
                    barycentricv = GetBarycentric(v1,v2,v3,glm::vec3(xi,yi,0.0));

                    if(barycentricv.x < 0 || barycentricv.x + barycentricv.y > 1) continue;

                    glm::vec3 pixelColor = barycentricv.x * color1 + barycentricv.y * color2 + barycentricv.z * color3;
                    float pixelZ = barycentricv.x * p1.z + barycentricv.y * p2.z + barycentricv.z * p3.z;
                    if(depthBuffer[xi+yi*WIDTH] == -1 || depthBuffer[xi+yi*WIDTH] > pixelZ)
                    {
                        result.setPixel(xi, yi, qRgb(pixelColor.x, pixelColor.y, pixelColor.z));
                        depthBuffer[xi+yi*WIDTH] = pixelZ;
                    }
                }
            }
        }
    }
    return result;
#else
    BoundingBox2D boundingBox;
    glm::vec3 barycentricv;
    glm::mat4 modelMat;
    glm::mat4 viewMat = myCamera.getViewMat();
    glm::mat4 projectionMat = myCamera.getProjectionMat();
    glm::mat4 NDC2Screen = Camera::getNDC2ScreenMat();
    for(const Polygon& currentPolygon:m_polygons)
    {
        QImage *currentTexture = currentPolygon.mp_texture;
        for (const Triangle &t : currentPolygon.m_tris)
        {
//            Vertex v1 = currentPolygon.m_verts[t.m_indices[0]];
//            Vertex v2 = currentPolygon.m_verts[t.m_indices[1]];
//            Vertex v3 = currentPolygon.m_verts[t.m_indices[2]];

//            glm::vec4 worldp1 = modelMat * v1.m_pos;
//            glm::vec4 worldp2 = modelMat * v2.m_pos;
//            glm::vec4 worldp3 = modelMat * v3.m_pos;

//            glm::vec4 camerap1 = viewMat * worldp1;
//            glm::vec4 camerap2 = viewMat * worldp2;
//            glm::vec4 camerap3 = viewMat * worldp3;

//            glm::vec4 projp1 = projectionMat * camerap1;
//            glm::vec4 projp2 = projectionMat * camerap2;
//            glm::vec4 projp3 = projectionMat * camerap3;

//            glm::vec4 p1 = projp1 / projp1[3];
//            glm::vec4 p2 = projp2 / projp2[3];
//            glm::vec4 p3 = projp3 / projp3[3];

//            p1[0] = (1+p1[0]) * result.width() /2;
//            p1[1] = (1-p1[1]) * result.height()/2;
//            p2[0] = (1+p2[0]) * result.width() /2;
//            p2[1] = (1-p2[1]) * result.height() /2;
//            p3[0] = (1+p3[0]) * result.width() /2;
//            p3[1] = (1-p3[1]) * result.height() /2;
//            std::array<LineSegment2D,3> l{LineSegment2D(p1,p2),LineSegment2D(p2,p3),LineSegment2D(p3,p1)};

            std::array<Vertex,3> v{currentPolygon.m_verts[t.m_indices[0]],currentPolygon.m_verts[t.m_indices[1]],currentPolygon.m_verts[t.m_indices[2]]};
            std::array<glm::vec2,4> uvMap{v[0].m_uv, v[1].m_uv, v[2].m_uv}, clipuvMap;
            std::array<glm::vec4,4> normals{v[0].m_normal, v[1].m_normal, v[2].m_normal}, clipNormals;
            std::array<glm::vec4,4> worldp{modelMat * v[0].m_pos, modelMat * v[1].m_pos, modelMat * v[2].m_pos};
            std::array<glm::vec4,4> camerap{viewMat * worldp[0], viewMat * worldp[1], viewMat * worldp[2]}, clipCamerap;
            std::vector<int> behindNearPlane, frontOfNearPlane;
            std::array<glm::vec4,4>projp;
            std::array<glm::vec4,4>screenp;
            for(int i = 0; i < 3; ++i)
            {
                if(camerap[i].z > -myCamera.getNear())
                {
                    behindNearPlane.push_back(i);
                }
                else
                {
                    frontOfNearPlane.push_back(i);
                }
            }
            int behindNearPlaneCount = behindNearPlane.size();
            int pointCount = ((behindNearPlaneCount == 1) ? 4 : 3);
            if(behindNearPlaneCount == 3) continue;

            else if(behindNearPlaneCount == 2)
            {
                float u1 = interpolateZ(camerap[frontOfNearPlane[0]], camerap[behindNearPlane[0]], -myCamera.getNear());
                float u2 = interpolateZ(camerap[frontOfNearPlane[0]], camerap[behindNearPlane[1]], -myCamera.getNear());
                glm::vec4 tempp1 = (1-u1)*camerap[frontOfNearPlane[0]] + u1*camerap[behindNearPlane[0]], tempp2 = (1-u2)*camerap[frontOfNearPlane[0]] + u2*camerap[behindNearPlane[1]];
                glm::vec2 uv1 = (1-u1)*uvMap[frontOfNearPlane[0]] + u1*uvMap[behindNearPlane[0]], uv2 = (1-u2)*uvMap[frontOfNearPlane[0]] + u2*uvMap[behindNearPlane[1]];
                glm::vec4 normal1 = (1-u1)*normals[frontOfNearPlane[0]] + u1*normals[behindNearPlane[0]], normal2 = (1-u2)*normals[frontOfNearPlane[0]] + u2*normals[behindNearPlane[1]];
                clipCamerap = std::array<glm::vec4, 4>{camerap[frontOfNearPlane[0]], tempp1, tempp2};
                projp = std::array<glm::vec4,4>{projectionMat * clipCamerap[0], projectionMat * clipCamerap[1], projectionMat * clipCamerap[2]};
                clipuvMap = std::array<glm::vec2,4>{uvMap[frontOfNearPlane[0]], uv1, uv2};
                clipNormals = std::array<glm::vec4,4>{normals[frontOfNearPlane[0]], normal1, normal2};
            }

            else if(behindNearPlaneCount == 1)
            {
                float u1 = interpolateZ(camerap[frontOfNearPlane[0]], camerap[behindNearPlane[0]], -myCamera.getNear());
                float u2 = interpolateZ(camerap[frontOfNearPlane[1]], camerap[behindNearPlane[0]], -myCamera.getNear());
                glm::vec4 tempp1 = (1-u1)*camerap[frontOfNearPlane[0]] + u1*camerap[behindNearPlane[0]], tempp2 = (1-u2)*camerap[frontOfNearPlane[1]] + u2*camerap[behindNearPlane[0]];
                glm::vec2 uv1 = (1-u1)*uvMap[frontOfNearPlane[0]] + u1*uvMap[behindNearPlane[0]], uv2 = (1-u2)*uvMap[frontOfNearPlane[1]] + u2*uvMap[behindNearPlane[0]];
                glm::vec4 normal1 = (1-u1)*normals[frontOfNearPlane[0]] + u1*normals[behindNearPlane[0]], normal2 = (1-u2)*normals[frontOfNearPlane[1]] + u2*normals[behindNearPlane[0]];
                clipCamerap = std::array<glm::vec4, 4>{camerap[frontOfNearPlane[0]], tempp1, tempp2, camerap[frontOfNearPlane[1]]};
                projp = std::array<glm::vec4,4>{projectionMat * clipCamerap[0], projectionMat * clipCamerap[1], projectionMat * clipCamerap[2], projectionMat * clipCamerap[3]};
                clipuvMap = std::array<glm::vec2,4>{uvMap[frontOfNearPlane[0]], uv1, uv2, uvMap[frontOfNearPlane[1]]};
                clipNormals = std::array<glm::vec4,4>{normals[frontOfNearPlane[0]], normal1, normal2, normals[frontOfNearPlane[1]]};
            }

            else
            {
                projp = std::array<glm::vec4,4>{projectionMat * camerap[0], projectionMat * camerap[1], projectionMat * camerap[2]};
                clipuvMap = std::array<glm::vec2,4>(uvMap);
                clipCamerap = std::array<glm::vec4, 4>(camerap);
                clipNormals = std::array<glm::vec4, 4>(normals);
            }

            for(int i = 0; i < pointCount; ++i)
            {
                screenp[i] = projp[i]/projp[i][3];
            }
            for(int i = 0; i < pointCount; ++i)
            {
                screenp[i][0] = (1 + screenp[i][0]) * result.width()/2;
                screenp[i][1] = (1 - screenp[i][1]) * result.height()/2;
            }
            // the solide mode
            if(!this->LineMode)
            {
                std::array<glm::vec2,3> tempuvMap{clipuvMap[0], clipuvMap[1], clipuvMap[2]};
                std::array<glm::vec4,3> tempScreenp{screenp[0], screenp[1], screenp[2]};
                std::array<glm::vec4,3> tempCamerap{clipCamerap[0], clipCamerap[1], clipCamerap[2]};
                std::array<glm::vec4,3> tempNormals{};
                for(int drawi = 0; drawi < pointCount-2; ++drawi)
                {
                    tempuvMap = std::array<glm::vec2,3>{clipuvMap[0], clipuvMap[drawi+1], clipuvMap[drawi+2]};
                    tempScreenp = std::array<glm::vec4,3>{screenp[0], screenp[drawi+1], screenp[drawi+2]};
                    tempCamerap = std::array<glm::vec4,3>{clipCamerap[0], clipCamerap[drawi+1], clipCamerap[drawi+2]};
                    tempNormals = std::array<glm::vec4,3>{clipNormals[0], clipNormals[drawi+1], clipNormals[drawi+2]};
                    drawTriangle(lMode,myCamera, tempNormals, tempuvMap, tempScreenp, tempCamerap, result, AABuffer, depthBuffer,currentTexture);
                }
            }
            // the wireframe mode
            else
            {
                if((screenp[0].z > 1&&screenp[1].z > 1&&screenp[2].z > 1)||(screenp[0].z < -1&&screenp[1].z < -1&&screenp[2].z < -1)) continue;
                Bresenham(glm::vec4(screenp[0].x,screenp[0].y,camerap[0].z,screenp[0].z), glm::vec4(screenp[1].x,screenp[1].y,camerap[1].z,screenp[1].z), depthBuffer, AABuffer, glm::vec3(255.0f,255.0f,255.0f), result.width(), result.height());
                Bresenham(glm::vec4(screenp[0].x,screenp[0].y,camerap[0].z,screenp[0].z), glm::vec4(screenp[2].x,screenp[2].y,camerap[2].z,screenp[2].z), depthBuffer, AABuffer, glm::vec3(255.0f,255.0f,255.0f), result.width(), result.height());
                Bresenham(glm::vec4(screenp[2].x,screenp[2].y,camerap[2].z,screenp[2].z), glm::vec4(screenp[1].x,screenp[1].y,camerap[1].z,screenp[1].z), depthBuffer, AABuffer, glm::vec3(255.0f,255.0f,255.0f), result.width(), result.height());
            }
        }
    }
    result = QImage(WIDTH, HEIGHT, QImage::Format_RGB32);
    result.fill(0);
    for(int yi = 0; yi < result.height(); ++yi)
    {
        for(int xi = 0; xi < result.width(); ++xi)
        {
            glm::vec3 color(0);
            for(int i = 0; i < AA2; ++i)
            {
                int dx = i % AA;
                int dy = std::floor(i / AA);
                color += AABuffer[AA*xi + (AA*yi+dy)*(AA*result.width()) + dx];
            }
            color/=AA2;
            result.setPixel(xi , yi, qRgb(color.x, color.y, color.z));
        }
    }
    return result;
#endif
}


void Rasterizer::ClearScene()
{
    m_polygons.clear();
}

void Rasterizer::moveCamera(glm::vec4 dPos)
{
    this->myCamera.movePosition(dPos);
}

void Rasterizer::addCameraFOV(float dFOV)
{
    this->myCamera.addFOV(dFOV);
}

void Rasterizer::addRotation(glm::vec4 axis, float degree)
{
    myCamera.addRotation(axis, degree);
}

void getBoundingBox2D(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3, BoundingBox2D &box, glm::vec2 xlimit, glm::vec2 ylimit)
{
    float xmax,xmin,ymax,ymin;
    xmax = std::max(p1.x,p2.x);
    xmax = std::max(xmax, p3.x);
    xmax = std::min(xmax, xlimit[1]);

    ymax = std::max(p1.y,p2.y);
    ymax = std::max(ymax, p3.y);
    ymax = std::min(ymax, ylimit[1]);

    xmin = std::min(p1.x,p2.x);
    xmin = std::min(xmin, p3.x);
    xmin = std::max(xmin, xlimit[0]);

    ymin = std::min(p1.y,p2.y);
    ymin = std::min(ymin, p3.y);
    ymin = std::max(ymin, ylimit[0]);

    box.xmin = xmin;
    box.xmax = xmax;
    box.ymax = ymax;
    box.ymin = ymin;
    return;
}

void getBoundingBox2D(const glm::vec4& p1, const glm::vec4& p2, const glm::vec4& p3, BoundingBox2D &box,glm::vec2 xlimit, glm::vec2 ylimit)
{
    float xmax,xmin,ymax,ymin;
    xmax = std::max(p1.x,p2.x);
    xmax = std::max(xmax, p3.x);
    xmax = std::min(xmax, xlimit[1]);

    ymax = std::max(p1.y,p2.y);
    ymax = std::max(ymax, p3.y);
    ymax = std::min(ymax, ylimit[1]);

    xmin = std::min(p1.x,p2.x);
    xmin = std::min(xmin, p3.x);
    xmin = std::max(xmin, xlimit[0]);

    ymin = std::min(p1.y,p2.y);
    ymin = std::min(ymin, p3.y);
    ymin = std::max(ymin, ylimit[0]);

    box.xmin = xmin;
    box.xmax = xmax;
    box.ymax = ymax;
    box.ymin = ymin;
    return;
}

void Rasterizer::setAntiAliasing(int AA)
{
    Rasterizer::AntiAliasing = AA;
}

// for p1,p2 (x,y,z in the camera space, z in the screen space)
void Bresenham(const glm::vec4 &p1, const glm::vec4 &p2, std::vector<float> &depthBuffer, std::vector<glm::vec3> &AABuffer, const glm::vec3 &color, const int width, const int height)
{
    glm::vec2 tempp1(std::round(p1.x), std::round(p1.y)), tempp2(std::round(p2.x), std::round(p2.y));
    bool reverse = false;
    float dy = tempp2.y - tempp1.y;
    float dx = tempp2.x - tempp1.x;
    if(std::abs(dy) > std::abs(dx))
    {
        reverse = true;
        std::swap(tempp1.x, tempp1.y);
        std::swap(tempp2.x, tempp2.y);
        std::swap(dy,dx);
    }
    int error = 0;
    int derror = 2 * abs(dy);
    dx = std::abs(dx);
    if(tempp1.x > tempp2.x) std::swap(tempp1, tempp2);
    int yi = tempp1.y;
    for(int xi = tempp1.x; xi <= tempp2.x; ++xi)
    {
        error += derror;
        if(!(xi > width || xi < 0 || yi > height || yi < 0))
        {
            float u = (xi - tempp1.x)/(tempp2.x - tempp1.x);
            float cameraZ = 1/((1-u)/p1.z + u/p2.z);
            float pixelZ = (1-u)*p1[3]/p1.z + u*p2[3]/p2.z;
            pixelZ *= cameraZ;
            int x = xi, y = yi;
            if(reverse) std::swap(x,y);
            if(pixelZ < depthBuffer[x + y * width])
            {
                depthBuffer[x + y * width] = pixelZ;
                AABuffer[x + y * width] = color;
            }
        }
        if(error > dx)
        {
            yi += (tempp2.y > tempp1.y ? 1:-1);
            error -= dx*2;
        }
    }
}

void drawTriangle(const lightMode lMode, const Camera &myCam, const std::array<glm::vec4,3> normals, const std::array<Vertex,3>v, const std::array<glm::vec4,3> &screenp, const std::array<glm::vec4, 3> &camerap, QImage &result, std::vector<glm::vec3> &AABuffer, std::vector<float> &depthBuffer, const QImage *currentTexture)
{
    if(screenp[0].x < 0 && screenp[1].x < 0 && screenp[2].x < 0) return;
    if(screenp[0].y < 0 && screenp[1].y < 0 && screenp[2].y < 0) return;
    if(screenp[0].x > result.width() && screenp[1].x > result.width() && screenp[2].x > result.width()) return;
    if(screenp[0].y > result.height() && screenp[1].y > result.height() && screenp[2].y > result.height()) return;
    BoundingBox2D boundingBox;
    getBoundingBox2D(screenp[0],screenp[1],screenp[2],boundingBox, glm::vec2(0,result.width()-1), glm::vec2(0,result.height()-1));
    if(boundingBox.xmax <= boundingBox.xmin || boundingBox.ymax <= boundingBox.ymin) return;
    std::array<LineSegment2D, 3> l{{LineSegment2D(screenp[0],screenp[1]),LineSegment2D(screenp[1],screenp[2]),LineSegment2D(screenp[0],screenp[2])}};
    for(int yi = std::floor(boundingBox.ymin); yi < std::ceil(boundingBox.ymax); ++yi)
    {
        float hitxmin = result.width() + 1, hitxmax = -1;
        std::array<float,3> hitx;
        std::array<bool,3> ifHit;
        ifHit[0] = l[0].getIntersection(yi,&hitx[0]);
        ifHit[1] = l[1].getIntersection(yi,&hitx[1]);
        ifHit[2] = l[2].getIntersection(yi,&hitx[2]);
        for(int j = 0; j < 3; ++j)
        {
            if(ifHit[j] && l[j].dy == 0)
            {
                hitxmin = std::clamp(std::min(l[j].pStart.x, l[j].pEnd.x),0.0f,float(result.width() - 1));
                hitxmax = std::clamp(std::max(l[j].pStart.x, l[j].pEnd.x),0.0f,float(result.width() - 1));
            }
            else if(ifHit[j])
            {
                hitx[j] = std::clamp(hitx[j], 0.0f, float(result.width() - 1));
                hitxmin = std::min(hitxmin, hitx[j]);
                hitxmax = std::max(hitxmax, hitx[j]);
            }
        }
        if(hitxmin > hitxmax) continue;
        for(int xi = std::round(hitxmin); xi <= std::round(hitxmax); ++xi)
        {
            glm::vec3 barycentricv = GetBarycentric(screenp[0],screenp[1],screenp[2],glm::vec3(xi,yi,0));
            if(barycentricv.z < 0 || barycentricv.y < 0 || barycentricv.y + barycentricv.z > 1) continue;

            float cameraZInv = barycentricv.x/camerap[0].z + barycentricv.y/camerap[1].z + barycentricv.z/camerap[2].z;
            float cameraZ = 1/cameraZInv;

            float pixelZ = barycentricv.x * screenp[0].z / camerap[0].z + barycentricv.y * screenp[1].z/ camerap[1].z + barycentricv.z * screenp[2].z / camerap[2].z;
            pixelZ *= cameraZ;

            if(pixelZ > depthBuffer[xi + yi* result.width()] || pixelZ >= 1 || pixelZ <= -1) continue;

            glm::vec2 pixeluv = barycentricv.x * v[0].m_uv / camerap[0].z + barycentricv.y * v[1].m_uv / camerap[1].z + barycentricv.z * v[2].m_uv / camerap[2].z;
            pixeluv *= cameraZ;

            glm::vec4 pixelNormal = barycentricv.x * normals[0] / camerap[0].z + barycentricv.y * normals[1] / camerap[1].z + barycentricv.z * normals[2] / camerap[2].z;
            pixelNormal *= cameraZ;
            pixelNormal = glm::normalize(pixelNormal);

            glm::vec3 pixelColor = GetImageColor(pixeluv, currentTexture);
            if(lMode == Lambertian)
            {
                pixelColor *= (glm::clamp(glm::dot(-myCam.getForward(), pixelNormal),.0f,1.0f) * 0.7f + 0.3f);
            }

            if(depthBuffer[xi+yi*result.width()] > pixelZ && (pixelZ <=1 && pixelZ > -1))
            {
                AABuffer[xi+yi*result.width()] = glm::vec3(pixelColor.x, pixelColor.y, pixelColor.z);
                //result.setPixel(xi, yi, qRgb(pixelColor.x, pixelColor.y, pixelColor.z));
                depthBuffer[xi+yi*result.width()] = pixelZ;
            }
        }
    }
}

void drawTriangle(const lightMode lMode, const Camera &myCam, const std::array<glm::vec4,3> normals, const std::array<glm::vec2,3>uvMap, const std::array<glm::vec4,3> &screenp, const std::array<glm::vec4, 3> &camerap, QImage &result, std::vector<glm::vec3> &AABuffer, std::vector<float> &depthBuffer, const QImage *currentTexture)
{
    if(screenp[0].x < 0 && screenp[1].x < 0 && screenp[2].x < 0) return;
    if(screenp[0].y < 0 && screenp[1].y < 0 && screenp[2].y < 0) return;
    if(screenp[0].x > result.width() && screenp[1].x > result.width() && screenp[2].x > result.width()) return;
    if(screenp[0].y > result.height() && screenp[1].y > result.height() && screenp[2].y > result.height()) return;
    BoundingBox2D boundingBox;
    getBoundingBox2D(screenp[0],screenp[1],screenp[2],boundingBox, glm::vec2(0,result.width()-1), glm::vec2(0,result.height()-1));
    if(boundingBox.xmax <= boundingBox.xmin || boundingBox.ymax <= boundingBox.ymin) return;
    std::array<LineSegment2D, 3> l{{LineSegment2D(screenp[0],screenp[1]),LineSegment2D(screenp[1],screenp[2]),LineSegment2D(screenp[0],screenp[2])}};
    for(int yi = std::floor(boundingBox.ymin); yi < std::ceil(boundingBox.ymax); ++yi)
    {
        float hitxmin = result.width() + 1, hitxmax = -1;
        std::array<float,3> hitx;
        std::array<bool,3> ifHit;
        ifHit[0] = l[0].getIntersection(yi,&hitx[0]);
        ifHit[1] = l[1].getIntersection(yi,&hitx[1]);
        ifHit[2] = l[2].getIntersection(yi,&hitx[2]);
        for(int j = 0; j < 3; ++j)
        {
            if(ifHit[j] && l[j].dy == 0)
            {
                hitxmin = std::clamp(std::min(l[j].pStart.x, l[j].pEnd.x),0.0f,float(result.width() - 1));
                hitxmax = std::clamp(std::max(l[j].pStart.x, l[j].pEnd.x),0.0f,float(result.width() - 1));
            }
            else if(ifHit[j])
            {
                hitx[j] = std::clamp(hitx[j], 0.0f, float(result.width() - 1));
                hitxmin = std::min(hitxmin, hitx[j]);
                hitxmax = std::max(hitxmax, hitx[j]);
            }
        }
        if(hitxmin > hitxmax) continue;
        for(int xi = std::round(hitxmin); xi <= std::round(hitxmax); ++xi)
        {
            glm::vec3 barycentricv = GetBarycentric(screenp[0],screenp[1],screenp[2],glm::vec3(xi,yi,0));
            if(barycentricv.z < 0 || barycentricv.y < 0 || barycentricv.y + barycentricv.z > 1) continue;

            float cameraZInv = barycentricv.x/camerap[0].z + barycentricv.y/camerap[1].z + barycentricv.z/camerap[2].z;
            float cameraZ = 1/cameraZInv;

            float pixelZ = barycentricv.x * screenp[0].z / camerap[0].z + barycentricv.y * screenp[1].z/ camerap[1].z + barycentricv.z * screenp[2].z / camerap[2].z;
            pixelZ *= cameraZ;

            if(pixelZ > depthBuffer[xi + yi* result.width()] || pixelZ > 1 || pixelZ < -1) continue;

            glm::vec2 pixeluv = barycentricv.x * uvMap[0] / camerap[0].z + barycentricv.y * uvMap[1] / camerap[1].z + barycentricv.z * uvMap[2] / camerap[2].z;
            pixeluv *= cameraZ;

            glm::vec4 pixelNormal = barycentricv.x * normals[0] / camerap[0].z + barycentricv.y * normals[1] / camerap[1].z + barycentricv.z * normals[2] / camerap[2].z;
            pixelNormal *= cameraZ;
            pixelNormal = glm::normalize(pixelNormal);

            glm::vec3 pixelColor = GetImageColor(pixeluv, currentTexture);
            if(lMode == Lambertian)
            {
                pixelColor *= (glm::clamp(glm::dot(-myCam.getForward(), pixelNormal),.0f,1.0f) * 0.7f + 0.3f);
            }

            if(depthBuffer[xi+yi*result.width()] > pixelZ && (pixelZ < 1 && pixelZ > -1))
            {
                AABuffer[xi+yi*result.width()] = glm::vec3(pixelColor.x, pixelColor.y, pixelColor.z);
                //result.setPixel(xi, yi, qRgb(pixelColor.x, pixelColor.y, pixelColor.z));
                depthBuffer[xi+yi*result.width()] = pixelZ;
            }
        }
    }
}

float interpolateZ(glm::vec3 p1, glm::vec3 p2, float targetZ)
{
    if(p2.z == p1.z) return -1;
    else return (targetZ-p1.z)/(p2.z-p1.z);
}

float interpolateZ(glm::vec4 p1, glm::vec4 p2, float targetZ)
{
    if(p2.z == p1.z) return -1;
    else return (targetZ-p1.z)/(p2.z-p1.z);
}

