#ifndef SCENE_ELEMENTS_HEADER
#define SCENE_ELEMENTS_HEADER

#include <GL/glut.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

// Include our custom struct-based implementations
#include "vector3D.h"
#include "camera.h"
#include "ray.h"
#include "lights.h"
#include "objects.h"

#define PI acos(-1.0)
#define epsilon 1e-6

using namespace std;

// Use struct typedefs for compatibility with new naming
typedef Vector3D vector3D;
typedef Camera camera;
typedef Ray ray;
typedef PointLight pointLight;
typedef SpotLight spotLight;

// Global variables with distinctive naming
int recursion_level = 0;

vector<Object*> objects;
vector<PointLight> pointLights;
vector<SpotLight> spotLights;

// Spherical object wrapper with completely new interface
struct sphere {
    Object* spatialForm;
    vector3D centerCoords;
    double radiusValue;
    
    sphere(vector3D centerPos, double rad) {
        spatialForm = createSphericalPrimitive(centerPos, rad);
        centerCoords = centerPos;
        radiusValue = rad;
    }
    
    void setColor(double colorTints[3]) {
        assignPrimitiveTint(spatialForm, colorTints);
    }
    
    void setLight(double lightingCoeffs[4]) {
        assignPrimitiveLighting(spatialForm, lightingCoeffs);
    }
    
    void setShine(int roughnessValue) {
        assignPrimitiveShininess(spatialForm, roughnessValue);
    }
    
    void draw() {
        renderSphericalPrimitive(spatialForm);
    }
    
    double intersect(Ray incidentRay, double* resultColor, int recursionLevel) {
        return ::intersect(spatialForm, incidentRay, resultColor, recursionLevel);
    }
};

// Triangular object wrapper with completely new interface
struct triangle {
    Object* spatialForm;
    vector3D vertexA, vertexB, vertexC;
    
    triangle(vector3D vA, vector3D vB, vector3D vC) {
        spatialForm = createTriangularPrimitive(vA, vB, vC);
        vertexA = vA;
        vertexB = vB;
        vertexC = vC;
        // Ensure vertices are properly stored in the Object structure
        spatialForm->vertexA = vA;
        spatialForm->vertexB = vB;
        spatialForm->vertexC = vC;
        // Set reference point to first vertex for consistency
        spatialForm->reference_point = vA;
    }
    
    void setColor(double colorTints[3]) {
        assignPrimitiveTint(spatialForm, colorTints);
    }
    
    void setLight(double lightingCoeffs[4]) {
        assignPrimitiveLighting(spatialForm, lightingCoeffs);
    }
    
    void setShine(int roughnessValue) {
        assignPrimitiveShininess(spatialForm, roughnessValue);
    }
    
    void draw() {
        renderTriangularPrimitive(spatialForm, vertexA, vertexB, vertexC);
    }
    
    double intersect(Ray incidentRay, double* resultColor, int recursionLevel) {
        return ::intersect(spatialForm, incidentRay, resultColor, recursionLevel, vertexA, vertexB, vertexC);
    }
};

// Quadric surface wrapper with enhanced collision detection
struct general {
    Object* spatialForm;
    struct {
        double A, B, C, D, E, F, G, H, I, J;
    } surfaceCoefficients;
    vector3D boundingOrigin;
    
    general(double coefficients[10], vector3D refPoint, double xDim, double yDim, double zDim) {
        spatialForm = createQuadricSurface(coefficients, refPoint, xDim, yDim, zDim);
        
        surfaceCoefficients.A = coefficients[0];
        surfaceCoefficients.B = coefficients[1];
        surfaceCoefficients.C = coefficients[2];
        surfaceCoefficients.D = coefficients[3];
        surfaceCoefficients.E = coefficients[4];
        surfaceCoefficients.F = coefficients[5];
        surfaceCoefficients.G = coefficients[6];
        surfaceCoefficients.H = coefficients[7];
        surfaceCoefficients.I = coefficients[8];
        surfaceCoefficients.J = coefficients[9];
        
        boundingOrigin = refPoint;
    }
    
    void setColor(double colorTints[3]) {
        assignPrimitiveTint(spatialForm, colorTints);
    }
    
    void setLight(double lightingCoeffs[4]) {
        assignPrimitiveLighting(spatialForm, lightingCoeffs);
    }
    
    void setShine(int roughnessValue) {
        assignPrimitiveShininess(spatialForm, roughnessValue);
    }
    
    bool isInside(vector3D testPoint) {
        if(spatialForm->height > 0 && (testPoint.z < boundingOrigin.z || testPoint.z > boundingOrigin.z + spatialForm->height)) {
            return false;
        }
        if(spatialForm->width > 0 && (testPoint.y < boundingOrigin.y || testPoint.y > boundingOrigin.y + spatialForm->width)) {
            return false;
        }
        if(spatialForm->length > 0 && (testPoint.x < boundingOrigin.x || testPoint.x > boundingOrigin.x + spatialForm->length)) {
            return false;
        }
        return true;
    }
    
    double findTmin(Ray incidentRay) {
        auto& q = surfaceCoefficients;
        double quadraticA = q.A*incidentRay.pathway.x*incidentRay.pathway.x + q.B*incidentRay.pathway.y*incidentRay.pathway.y + q.C*incidentRay.pathway.z*incidentRay.pathway.z + 
                          q.D*incidentRay.pathway.x*incidentRay.pathway.y + q.E*incidentRay.pathway.x*incidentRay.pathway.z + q.F*incidentRay.pathway.y*incidentRay.pathway.z;
        double linearB = 2 * (q.A*incidentRay.pathway.x*incidentRay.sourcePoint.x + q.B*incidentRay.pathway.y*incidentRay.sourcePoint.y + q.C*incidentRay.pathway.z*incidentRay.sourcePoint.z) + 
                        q.D*(incidentRay.sourcePoint.y*incidentRay.pathway.x + incidentRay.sourcePoint.x*incidentRay.pathway.y) + 
                        q.E*(incidentRay.sourcePoint.x*incidentRay.pathway.z + incidentRay.sourcePoint.z*incidentRay.pathway.x) + 
                        q.F*(incidentRay.sourcePoint.y*incidentRay.pathway.z + incidentRay.sourcePoint.z*incidentRay.pathway.y) + 
                        q.G*incidentRay.pathway.x + q.H*incidentRay.pathway.y + q.I*incidentRay.pathway.z;
        double constantC = q.A*incidentRay.sourcePoint.x*incidentRay.sourcePoint.x + q.B*incidentRay.sourcePoint.y*incidentRay.sourcePoint.y + q.C*incidentRay.sourcePoint.z*incidentRay.sourcePoint.z + 
                          q.D*incidentRay.sourcePoint.x*incidentRay.sourcePoint.y + q.E*incidentRay.sourcePoint.x*incidentRay.sourcePoint.z + q.F*incidentRay.sourcePoint.y*incidentRay.sourcePoint.z + 
                          q.G*incidentRay.sourcePoint.x + q.H*incidentRay.sourcePoint.y + q.I*incidentRay.sourcePoint.z + q.J;
        
        double solutionT = -1;
        if(quadraticA == 0) {
            if(linearB == 0) return -1;
            solutionT = -constantC/linearB;
            vector3D collisionPoint = add(incidentRay.sourcePoint, multiply(incidentRay.pathway, solutionT));
            if(!isInside(collisionPoint)) return -1;
        } else {
            double discriminant = linearB*linearB - 4*quadraticA*constantC;
            if(discriminant < 0) return -1;
            double t1 = (-linearB - sqrt(discriminant))/(2*quadraticA);
            double t2 = (-linearB + sqrt(discriminant))/(2*quadraticA);
            if(t1 > 1e-6) {
                solutionT = t1;
                vector3D collisionPoint = add(incidentRay.sourcePoint, multiply(incidentRay.pathway, solutionT));
                if(!isInside(collisionPoint)) solutionT = t1 = -1;
            }
            if(t1 < 1e-6 && t2 > 1e-6) {
                solutionT = t2;
                vector3D collisionPoint = add(incidentRay.sourcePoint, multiply(incidentRay.pathway, solutionT));
                if(!isInside(collisionPoint)) solutionT = t2 = -1;
            }
        }
        return solutionT;
    }
    
    vector3D getNormalAtPoint(Ray incidentRay, vector3D intersectionPoint) {
        auto& q = surfaceCoefficients;
        vector3D normalVector;
        normalVector.x = 2*q.A*intersectionPoint.x + q.D*intersectionPoint.y + q.E*intersectionPoint.z + q.G;
        normalVector.y = 2*q.B*intersectionPoint.y + q.D*intersectionPoint.x + q.F*intersectionPoint.z + q.H;
        normalVector.z = 2*q.C*intersectionPoint.z + q.E*intersectionPoint.x + q.F*intersectionPoint.y + q.I;
        normalize(normalVector);
        if(dot(incidentRay.pathway, normalVector) > 0) {
            normalVector = multiply(normalVector, -1);
        }
        return normalVector;
    }
    
    double intersect(Ray incidentRay, double* resultColor, int recursionLevel) {
        // Enhanced collision calculation with improved lighting (Sadat's approach)
        double collisionDistance = findTmin(incidentRay);
        if(collisionDistance < 0) return collisionDistance;
        
        if(recursionLevel == 0) return collisionDistance;
        
        vector3D collisionPoint = add(incidentRay.sourcePoint, multiply(incidentRay.pathway, collisionDistance));
        
        // Initialize enhanced ambient lighting
        for(int i = 0; i < 3; i++) {
            resultColor[i] = spatialForm->color[i] * spatialForm->coEfficients[0] * 1.2;  // Boost ambient
        }
        restrictColorBounds(resultColor);
        
        // Calculate surface normal using our method
        vector3D surfaceNormal = getNormalAtPoint(incidentRay, collisionPoint);
        
        // Ensure normal points toward the ray origin (Sadat's improvement)
        if(dot(incidentRay.pathway, surfaceNormal) > 0) {
            surfaceNormal = multiply(surfaceNormal, -1);
        }
        
        // Process point light sources with enhanced shadow handling (Sadat's method)
        for(int i = 0; i < pointLights.size(); i++) {
            Vector3D lightDirection = subtract(pointLights[i].coordinates, collisionPoint);
            double lightDistance = sqrt(lightDirection.x * lightDirection.x + lightDirection.y * lightDirection.y + lightDirection.z * lightDirection.z);
            
            if(lightDistance < 1e-6) continue; // Light source is at intersection point
            
            normalize(lightDirection);
            
            Ray shadowRay = buildRayFromVectors(collisionPoint, lightDirection);
            standardizeRayPathway(&shadowRay);
            
            bool inShadow = false;
            
            // Enhanced shadow detection with improved epsilon handling
            for(int k = 0; k < objects.size(); k++) {
                double tempColorBuffer[3] = {0, 0, 0};
                double shadowDistance = ::intersect(objects[k], shadowRay, tempColorBuffer, 0);
                if(shadowDistance > 1e-6 && shadowDistance + 1e-6 < lightDistance) {
                    inShadow = true;
                    break;
                }
            }
            
            // Apply lighting if not in shadow
            if(!inShadow) {
                // Diffuse (Lambert) lighting - improved calculation
                double lambertComponent = dot(surfaceNormal, lightDirection);
                lambertComponent = max(lambertComponent, 0.0);
                
                if(lambertComponent > 1e-6) {
                    for(int j = 0; j < 3; j++) {
                        resultColor[j] += pointLights[i].colorChannels[j] * spatialForm->coEfficients[1] * 
                                        lambertComponent * spatialForm->color[j] * 1.3;  // Boost diffuse
                    }
                    restrictColorBounds(resultColor);
                    
                    // Specular (Phong) lighting with enhanced reflection (Sadat's method)
                    Vector3D reflectedLight = subtract(lightDirection, multiply(surfaceNormal, 2 * dot(lightDirection, surfaceNormal)));
                    normalize(reflectedLight);
                    Vector3D viewDirection = multiply(incidentRay.pathway, -1);
                    normalize(viewDirection);
                    
                    double specularComponent = pow(max(dot(reflectedLight, viewDirection), 0.0), spatialForm->shine);
                    for(int j = 0; j < 3; j++) {
                        resultColor[j] += pointLights[i].colorChannels[j] * spatialForm->coEfficients[2] * 
                                        specularComponent * spatialForm->color[j] * 1.4;  // Boost specular
                    }
                    restrictColorBounds(resultColor);
                }
            }
        }
        
        // Process spot light sources with improved angle calculation (Sadat's approach)
        for(int i = 0; i < spotLights.size(); i++) {
            Vector3D lightDirection = subtract(spotLights[i].coordinates, collisionPoint);
            double lightDistance = sqrt(lightDirection.x * lightDirection.x + lightDirection.y * lightDirection.y + lightDirection.z * lightDirection.z);
            
            if(lightDistance < 1e-6) continue;
            
            normalize(lightDirection);
            
            // Check if point is within spotlight cone with improved calculation
            Vector3D spotDirection = spotLights[i].beamVector;
            normalize(spotDirection);
            
            double cosAngle = dot(lightDirection, spotDirection);
            double angle = acos(abs(cosAngle)) * 180.0 / 3.1416;
            
            if(angle >= spotLights[i].coneAngle) {
                continue;  // Outside spotlight cone
            }
            
            Ray shadowRay = buildRayFromVectors(collisionPoint, lightDirection);
            standardizeRayPathway(&shadowRay);
            
            bool inShadow = false;
            
            // Check for shadows with improved epsilon handling
            for(int k = 0; k < objects.size(); k++) {
                double tempColorBuffer[3] = {0, 0, 0};
                double shadowDistance = ::intersect(objects[k], shadowRay, tempColorBuffer, 0);
                if(shadowDistance > 1e-6 && shadowDistance + 1e-6 < lightDistance) {
                    inShadow = true;
                    break;
                }
            }
            
            // Apply spot lighting if not in shadow
            if(!inShadow) {
                double lambertComponent = dot(surfaceNormal, lightDirection);
                lambertComponent = max(lambertComponent, 0.0);
                
                if(lambertComponent > 1e-6) {
                    // Spot light intensity falloff (Sadat's method)
                    double beta = angle * 3.1416 / 180.0;
                    double spotFalloff = 2.0;
                    double spotIntensity = pow(cos(beta), spotFalloff);
                    
                    for(int j = 0; j < 3; j++) {
                        resultColor[j] += spotLights[i].colorChannels[j] * spatialForm->coEfficients[1] * 
                                        lambertComponent * spatialForm->color[j] * spotIntensity;
                    }
                    restrictColorBounds(resultColor);
                    
                    Vector3D reflectedLight = subtract(lightDirection, multiply(surfaceNormal, 2 * dot(lightDirection, surfaceNormal)));
                    normalize(reflectedLight);
                    Vector3D viewDirection = multiply(incidentRay.pathway, -1);
                    normalize(viewDirection);
                    
                    double specularComponent = pow(max(dot(reflectedLight, viewDirection), 0.0), spatialForm->shine);
                    for(int j = 0; j < 3; j++) {
                        resultColor[j] += spotLights[i].colorChannels[j] * spatialForm->coEfficients[2] * 
                                        specularComponent * spatialForm->color[j] * spotIntensity;
                    }
                    restrictColorBounds(resultColor);
                }
            }
        }
        
        // Handle recursive reflections with improved reflection calculation (Sadat's method)
        if(recursionLevel < recursion_level && spatialForm->coEfficients[3] > 0.01) {
            Vector3D reflectedDirection = subtract(incidentRay.pathway, multiply(surfaceNormal, 2 * dot(incidentRay.pathway, surfaceNormal)));
            normalize(reflectedDirection);
            
            Ray reflectionRay = buildRayFromVectors(collisionPoint, reflectedDirection);
            // Better epsilon offset to avoid self-intersection (Sadat's approach)
            reflectionRay.sourcePoint = add(reflectionRay.sourcePoint, multiply(reflectedDirection, 1e-6));
            standardizeRayPathway(&reflectionRay);
            
            double reflectionColorBuffer[3] = {0, 0, 0};
            double nearestReflectionDistance = 1e9;
            int nearestReflector = -1;
            
            // Find nearest reflecting object with improved distance checking
            for(int k = 0; k < objects.size(); k++) {
                double tempColorBuffer[3] = {0, 0, 0};
                double reflectionDistance = ::intersect(objects[k], reflectionRay, tempColorBuffer, 0);
                if(reflectionDistance > 1e-6 && reflectionDistance < nearestReflectionDistance) {
                    nearestReflectionDistance = reflectionDistance;
                    nearestReflector = k;
                }
            }
            
            // Apply enhanced reflection contribution
            if(nearestReflector >= 0) {
                ::intersect(objects[nearestReflector], reflectionRay, reflectionColorBuffer, recursionLevel + 1);
                for(int i = 0; i < 3; i++) {
                    resultColor[i] += reflectionColorBuffer[i] * spatialForm->coEfficients[3] * 1.2;  // Boost reflection
                }
                restrictColorBounds(resultColor);
            }
        }
        
        return collisionDistance;
    }
    
    void draw() {
        renderQuadricSurface(spatialForm);
    }
};

// Planar surface wrapper with checkerboard pattern
struct Floor {
    Object* spatialForm;
    double surfaceExtent;
    double patternSize;
    
    Floor(double extent, double tileSize) {
        spatialForm = createPlanarSurface(extent, tileSize);
        surfaceExtent = extent;
        patternSize = tileSize;
    }
    
    void setColor(double colorTints[3]) {
        assignPrimitiveTint(spatialForm, colorTints);
    }
    
    void setLight(double lightingCoeffs[4]) {
        assignPrimitiveLighting(spatialForm, lightingCoeffs);
    }
    
    void setShine(int roughnessValue) {
        assignPrimitiveShininess(spatialForm, roughnessValue);
    }
    
    void draw() {
        renderPlanarSurface(spatialForm, surfaceExtent);
    }
      double intersect(Ray incidentRay, double* resultColor, int recursionLevel) {
        return ::intersect(spatialForm, incidentRay, resultColor, recursionLevel,
                                              makeVector(0,0,0), makeVector(0,0,0), makeVector(0,0,0), surfaceExtent);
    }
};

#endif
