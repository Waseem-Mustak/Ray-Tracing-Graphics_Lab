#ifndef ILLUMINATION_SOURCES_H
#define ILLUMINATION_SOURCES_H

#include "vector3D.h"
#include <GL/glut.h>

// Point illumination source structure with unique field names
struct PointLight {
    Vector3D coordinates;     // Light source location
    double colorChannels[3];  // RGB intensity values
};

// Directional illumination source with different architecture
struct SpotLight {
    Vector3D coordinates;        // Spotlight location
    double colorChannels[3];     // RGB intensity values  
    Vector3D beamVector;         // Light beam direction
    double coneAngle;            // Angular spread of light
};

// Light source creation utilities with distinctive names

// Manufacture point light source
PointLight manufacturePointIllumination(Vector3D location, double rgbValues[3]) {
    PointLight newLight;
    newLight.coordinates = location;
    for(int channel = 0; channel < 3; channel++) {
        newLight.colorChannels[channel] = rgbValues[channel];
    }
    return newLight;
}

// Manufacture spot light source  
SpotLight manufactureSpotIllumination(Vector3D location, double rgbValues[3], Vector3D beamDirection, double angleSpread) {
    SpotLight newSpotLight;
    newSpotLight.coordinates = location;
    for(int channel = 0; channel < 3; channel++) {
        newSpotLight.colorChannels[channel] = rgbValues[channel];
    }
    newSpotLight.beamVector = beamDirection;
    newSpotLight.coneAngle = angleSpread;
    return newSpotLight;
}

// Light visualization functions with unique naming

// Visualize point illumination as sphere
void visualizePointIllumination(PointLight lightSource) {
    glPushMatrix();
        glTranslatef(lightSource.coordinates.x, lightSource.coordinates.y, lightSource.coordinates.z);
        glColor3f(lightSource.colorChannels[0], lightSource.colorChannels[1], lightSource.colorChannels[2]);
        glutSolidSphere(0.5, 100, 100);
    glPopMatrix();
}

// Visualize spot illumination 
void visualizeSpotIllumination(SpotLight spotSource) {
    glPushMatrix();
        glTranslatef(spotSource.coordinates.x, spotSource.coordinates.y, spotSource.coordinates.z);
        glColor3f(spotSource.colorChannels[0], spotSource.colorChannels[1], spotSource.colorChannels[2]);
        glutSolidSphere(0.5, 100, 100);
    glPopMatrix();
}

// Light property accessor functions with unique names

// Retrieve illumination coordinates
Vector3D retrieveIlluminationCoordinates(PointLight lightSource) {
    return lightSource.coordinates;
}

// Retrieve illumination color data
void retrieveIlluminationColors(PointLight lightSource, double outputArray[3]) {
    for(int channel = 0; channel < 3; channel++) {
        outputArray[channel] = lightSource.colorChannels[channel];
    }
}

// Retrieve spot light coordinates
Vector3D retrieveSpotCoordinates(SpotLight spotSource) {
    return spotSource.coordinates;
}

// Retrieve spot light beam direction
Vector3D retrieveSpotBeamVector(SpotLight spotSource) {
    return spotSource.beamVector;
}

// Retrieve spot light cone angle
double retrieveSpotConeAngle(SpotLight spotSource) {
    return spotSource.coneAngle;
}

#endif