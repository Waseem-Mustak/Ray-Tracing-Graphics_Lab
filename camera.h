#ifndef CAMERA_H
#define CAMERA_H

#include "vector3D.h"
#include <cmath>
using namespace std;

// Camera variables
float camPos[3] = {4, 3, 4};
float lookAt[3] = {0, 0, 0};
float up[3] = {0, 1, 0};

float rotSpeed = 0.08;    // Slower for look rotation (1-6 keys)
float moveSpeed = 0.08;   // Faster for movement (arrow keys)
float fastMoveSpeed = 1.0; // Even faster for special keys

// Helper function for cross product
void cross(float a[3], float b[3], float result[3]) 
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

// Helper function for normalization
void normalize(float v[3]) 
{
    float len = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if(len == 0) return;
    v[0] /= len; v[1] /= len; v[2] /= len;
}

// Rotate view around an axis
void rotateView(float angle, float axis[3]) 
{
    float viewDir[3] = {lookAt[0]-camPos[0], lookAt[1]-camPos[1], lookAt[2]-camPos[2]};
    float c = cos(angle), s = sin(angle);
    float dot = axis[0]*viewDir[0] + axis[1]*viewDir[1] + axis[2]*viewDir[2];

    float newDir[3] = {
        viewDir[0]*c + (axis[1]*viewDir[2] - axis[2]*viewDir[1])*s + axis[0]*dot*(1-c),
        viewDir[1]*c + (axis[2]*viewDir[0] - axis[0]*viewDir[2])*s + axis[1]*dot*(1-c),
        viewDir[2]*c + (axis[0]*viewDir[1] - axis[1]*viewDir[0])*s + axis[2]*dot*(1-c)
    };

    lookAt[0] = camPos[0] + newDir[0];
    lookAt[1] = camPos[1] + newDir[1];
    lookAt[2] = camPos[2] + newDir[2];
}

// Rotate view and update UP vector
void rotateViewAndUpdateUP(float angle, float axis[3]) 
{
    float viewDir[3] = {lookAt[0]-camPos[0], lookAt[1]-camPos[1], lookAt[2]-camPos[2]};
    float c = cos(angle), s = sin(angle);
    float dot = axis[0]*viewDir[0] + axis[1]*viewDir[1] + axis[2]*viewDir[2];

    float newDir[3] = {
        viewDir[0]*c + (axis[1]*viewDir[2] - axis[2]*viewDir[1])*s + axis[0]*dot*(1-c),
        viewDir[1]*c + (axis[2]*viewDir[0] - axis[0]*viewDir[2])*s + axis[1]*dot*(1-c),
        viewDir[2]*c + (axis[0]*viewDir[1] - axis[1]*viewDir[0])*s + axis[2]*dot*(1-c)
    };

    lookAt[0] = camPos[0] + newDir[0];
    lookAt[1] = camPos[1] + newDir[1];
    lookAt[2] = camPos[2] + newDir[2];

    float dot1=(axis[0]*up[0] + axis[1]*up[1] + axis[2]*up[2]);
    float newUp[3] = {
        up[0]*c + (axis[1]*up[2] - axis[2]*up[1])*s + axis[0]*dot1*(1-c),
        up[1]*c + (axis[2]*up[0] - axis[0]*up[2])*s + axis[1]*dot1*(1-c),
        up[2]*c + (axis[0]*up[1] - axis[1]*up[0])*s + axis[2]*dot1*(1-c)
    };
    for(int i=0; i<3; i++) up[i] = newUp[i];
    normalize(up);
}

// Change camera position and look at position together
void changeCamPosAndLookAtPos(float dx, float dy, float dz) 
{
    camPos[0] += dx;
    camPos[1] += dy; 
    camPos[2] += dz; 
    lookAt[0] += dx;
    lookAt[1] += dy;
    lookAt[2] += dz;
}

// Change only camera position
void changeCamPos(float dx, float dy, float dz) 
{
    camPos[0] += dx;
    camPos[1] += dy; 
    camPos[2] += dz; 
}

// Get normalized view direction vector
void getViewDirection(float viewDir[3]) 
{
    viewDir[0] = lookAt[0] - camPos[0];
    viewDir[1] = lookAt[1] - camPos[1];
    viewDir[2] = lookAt[2] - camPos[2];
    normalize(viewDir);
}

// Get normalized right vector
void getRightVector(float right[3]) 
{
    float viewDir[3];
    getViewDirection(viewDir);
    cross(viewDir, up, right);
    normalize(right);
}

// Camera movement functions for keyboard input
void rotateCameraLeft() 
{
    rotateView(rotSpeed, up);
}

void rotateCameraRight() 
{
    rotateView(-rotSpeed, up);
}

void rotateCameraUp() 
{
    float right[3];
    getRightVector(right);
    rotateViewAndUpdateUP(rotSpeed, right);
}

void rotateCameraDown() 
{
    float right[3];
    getRightVector(right);
    rotateViewAndUpdateUP(-rotSpeed, right);
}

void rollCameraLeft() 
{
    float viewDir[3];
    getViewDirection(viewDir);
    rotateViewAndUpdateUP(rotSpeed, viewDir);
}

void rollCameraRight() 
{
    float viewDir[3];
    getViewDirection(viewDir);
    rotateViewAndUpdateUP(-rotSpeed, viewDir);
}

void moveCameraUp() 
{
    changeCamPos(0, 0, moveSpeed);
}

void moveCameraDown() 
{
    changeCamPos(0, 0, -moveSpeed);
}

void moveCameraForward() 
{
    float viewDir[3];
    getViewDirection(viewDir);
    changeCamPosAndLookAtPos(viewDir[0]*moveSpeed, viewDir[1]*moveSpeed, viewDir[2]*moveSpeed);
}

void moveCameraBackward() 
{
    float viewDir[3];
    getViewDirection(viewDir);
    changeCamPosAndLookAtPos(-viewDir[0]*moveSpeed, -viewDir[1]*moveSpeed, -viewDir[2]*moveSpeed);
}

void moveCameraLeft() 
{
    float right[3];
    getRightVector(right);
    changeCamPosAndLookAtPos(-right[0]*moveSpeed, -right[1]*moveSpeed, -right[2]*moveSpeed);
}

void moveCameraRight() 
{
    float right[3];
    getRightVector(right);
    changeCamPosAndLookAtPos(right[0]*moveSpeed, right[1]*moveSpeed, right[2]*moveSpeed);
}

void moveCameraUpGlobal() 
{
    changeCamPosAndLookAtPos(up[0]*moveSpeed, up[1]*moveSpeed, up[2]*moveSpeed);
}

void moveCameraDownGlobal() 
{
    changeCamPosAndLookAtPos(-up[0]*moveSpeed, -up[1]*moveSpeed, -up[2]*moveSpeed);
}

// Fast movement functions for special keys
void moveCameraForwardFast() 
{
    float viewDir[3];
    getViewDirection(viewDir);
    changeCamPosAndLookAtPos(viewDir[0]*fastMoveSpeed, viewDir[1]*fastMoveSpeed, viewDir[2]*fastMoveSpeed);
}

void moveCameraBackwardFast() 
{
    float viewDir[3];
    getViewDirection(viewDir);
    changeCamPosAndLookAtPos(-viewDir[0]*fastMoveSpeed, -viewDir[1]*fastMoveSpeed, -viewDir[2]*fastMoveSpeed);
}

void moveCameraLeftFast() 
{
    float right[3];
    getRightVector(right);
    changeCamPosAndLookAtPos(-right[0]*fastMoveSpeed, -right[1]*fastMoveSpeed, -right[2]*fastMoveSpeed);
}

void moveCameraRightFast() 
{
    float right[3];
    getRightVector(right);
    changeCamPosAndLookAtPos(right[0]*fastMoveSpeed, right[1]*fastMoveSpeed, right[2]*fastMoveSpeed);
}

void moveCameraUpGlobalFast() 
{
    changeCamPosAndLookAtPos(up[0]*fastMoveSpeed, up[1]*fastMoveSpeed, up[2]*fastMoveSpeed);
}

void moveCameraDownGlobalFast() 
{
    changeCamPosAndLookAtPos(-up[0]*fastMoveSpeed, -up[1]*fastMoveSpeed, -up[2]*fastMoveSpeed);
}

// Simple Camera struct - just holds vectors without redundant wrapper methods
#include "vector3D.h"

struct Camera {
    Vector3D pos, l, u, r;
    double rotationAngle;
    
    Camera() {
        updateFromGlobalCamera();
        rotationAngle = 3.14159/180.0;
    }
    
    // Update local vectors from global camera variables
    void updateFromGlobalCamera() {
        pos = makeVector(camPos[0], camPos[1], camPos[2]);
        Vector3D lookAtPos = makeVector(lookAt[0], lookAt[1], lookAt[2]);
        l = subtract(lookAtPos, pos);
        normalize(l);
        u = makeVector(up[0], up[1], up[2]);
        normalize(u);
        r = cross(l, u);
        normalize(r);
    }
    
    void setCoords(Vector3D a, Vector3D b, Vector3D c, double ang=3.14159/180.0) {
        pos = a;
        l = b;
        normalize(l);
        u = c;
        normalize(u);
        r = cross(l, u);
        normalize(r);
        rotationAngle = ang;
        // Update global variables
        camPos[0] = pos.x; camPos[1] = pos.y; camPos[2] = pos.z;
        lookAt[0] = pos.x + l.x; lookAt[1] = pos.y + l.y; lookAt[2] = pos.z + l.z;
        ::up[0] = u.x; ::up[1] = u.y; ::up[2] = u.z;
    }
};

#endif
