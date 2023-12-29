#include "camera.h"

Camera::Camera(glm::vec4 pos, glm::vec4 fDir,glm::vec4 uDir,glm::vec4 rDir, float inputFOV, float inputNear, float inputFar, float inputAspectRatio)
    :position(pos), forwardDir(fDir), upDir(uDir), rightDir(rDir), fov(inputFOV), near(inputNear), far(inputFar), aspectRatio(inputAspectRatio)
{
    zAxis = glm::normalize(-1.0f * forwardDir);
    yAxis = upDir - glm::dot(upDir,zAxis) * zAxis;
    yAxis = glm::normalize(yAxis);
    xAxis = glm::vec4(glm::cross(glm::vec3(yAxis), glm::vec3(zAxis)),0);
    xAxis = glm::normalize(xAxis);
}


// matrix translate a point from world coordinate system to camera system
glm::mat4 Camera::getViewMat()
{
    glm::mat3 axesMatT((glm::vec3(xAxis)), (glm::vec3(yAxis)), (glm::vec3(zAxis)));
    axesMatT = glm::transpose(axesMatT);

    glm::mat4 viewMat = glm::mat4(axesMatT);

    glm::mat4 translateView(glm::vec4(1, 0, 0, 0),
                            glm::vec4(0, 1, 0, 0),
                            glm::vec4(0, 0, 1, 0),
                            glm::vec4(-glm::vec3(position), 1));

    viewMat = viewMat * translateView;
    return viewMat;
}

glm::mat4 Camera::getProjectionMat()
{
    glm::mat4 projectionMat;
    float t = tan(MY_PI * (fov/2)/180) * (near);
    float b = -t;
    float r = t/aspectRatio;
    float l = -r;
    projectionMat[0] = glm::vec4(2*near/(r-l), 0, 0, 0);
    projectionMat[1] = glm::vec4(0, 2*near/(t-b), 0, 0);
    projectionMat[2] = glm::vec4((l+r)/(r-l), (t+b)/(t-b), (near+far)/(near-far), -1);
    projectionMat[3] = glm::vec4(0,0,2*near*far/(near-far),0);
    return projectionMat;
}


glm::mat4 Camera::getNDC2ScreenMat()
{
    glm::mat4 NDC2ScreenMat;
    NDC2ScreenMat[0] = glm::vec4(WIDTH/2,0,0,0);
    NDC2ScreenMat[1] = glm::vec4(0,-HEIGHT/2,0,0);
    NDC2ScreenMat[3] = glm::vec4(WIDTH/2,HEIGHT/2,0,1.0f);
    return NDC2ScreenMat;
}

void Camera::movePosition(glm::vec4 dpos)
{
    glm::mat4 axes(xAxis,yAxis,-zAxis,glm::vec4(0,0,0,1));
    this->position += axes * dpos;
}

void Camera::addFOV(float dFOV)
{
    this->fov += dFOV;
    this->fov = std::max(this->fov, 30.0f);
    this->fov = std::min(this->fov, 120.0f);
}

void Camera::addRotation(glm::vec4 axis, float degree)
{
    glm::mat4 rotationMat(1); // Creates a identity matrix
    rotationMat = glm::rotate(rotationMat, float(degree), glm::vec3(axis));
    glm::mat4 axes(xAxis,yAxis,zAxis,glm::vec4(0,0,0,1));
    axes *= rotationMat;
    rightDir = xAxis = axes[0];
    upDir = yAxis = axes[1];
    zAxis = axes[2];
    forwardDir = -zAxis;
}

float Camera::getNear()const
{
    return this->near;
}

float Camera::getFar()const
{
    return this->far;
}

glm::vec4 Camera::getForward()const
{
    return this->forwardDir;
}
