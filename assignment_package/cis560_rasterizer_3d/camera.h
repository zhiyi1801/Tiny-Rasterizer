#ifndef CAMERA_H
#define CAMERA_H
#include<iostream>
#include<glm/glm.hpp>
#include<glm/gtx/transform.hpp>
#include "macro.h"

class Camera
{
private:
    glm::vec4 position{0,0,10,1};
    glm::vec4 forwardDir{0,0,-1,0}, upDir{0,1,0,0}, rightDir{1,0,0,0};
    glm::vec4 zAxis, yAxis, xAxis;
    float fov = 60;
    float near = 0.01f,far = 100.0f;
    float aspectRatio = 1.0f;
public:
    Camera(glm::vec4 pos = glm::vec4(0,0,10,1), glm::vec4 fDir = glm::vec4{0,0,-1,0},glm::vec4 uDir = glm::vec4{0,1,0,0},glm::vec4 rDir = glm::vec4{1,0,0,0}, float inputFOV = 60, float inputNear = 0.01f, float inputFar = 100.0f, float inputAspectRatio = 1.0f);

    glm::mat4 getViewMat();
    glm::mat4 getProjectionMat();
    static glm::mat4 getNDC2ScreenMat();

    glm::mat4 getAxes();
    void movePosition(glm::vec4 dpos);
    void addFOV(float dFOV);
    void addRotation(glm::vec4 axis, float degree);
    float getNear()const;
    float getFar()const;
    glm::vec4 getForward()const;
};

#endif // CAMERA_H
