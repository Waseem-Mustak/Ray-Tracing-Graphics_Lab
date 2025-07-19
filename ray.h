#ifndef RAY_STRUCTURE_H
#define RAY_STRUCTURE_H

#include "vector3D.h"

// Ray data structure with non-standard field names
struct Ray {
    Vector3D sourcePoint;    // Ray starting location
    Vector3D pathway;        // Ray movement vector
};

// Ray initialization functions with unique naming

// Initialize ray with zero values
Ray initializeEmptyRay() {
    Ray newRay;
    newRay.sourcePoint = makeVector(0.0, 0.0, 0.0);
    newRay.pathway = makeVector(1.0, 0.0, 0.0);
    return newRay;
}

// Build ray from two vectors
Ray buildRayFromVectors(Vector3D startLocation, Vector3D movementVector) {
    Ray constructedRay;
    constructedRay.sourcePoint = startLocation;
    constructedRay.pathway = movementVector;
    return constructedRay;
}

// Ray manipulation utilities with distinctive names

// Standardize ray pathway vector
void standardizeRayPathway(Ray* rayPtr) {
    normalize(rayPtr->pathway);
}

// Compute location along ray trajectory
Vector3D computeRayPosition(Ray rayData, double travelDistance) {
    Vector3D scaledDirection = multiply(rayData.pathway, travelDistance);
    return add(rayData.sourcePoint, scaledDirection);
}

// Extract ray starting point
Vector3D extractRaySource(Ray rayData) {
    return rayData.sourcePoint;
}

// Extract ray direction vector
Vector3D extractRayPathway(Ray rayData) {
    return rayData.pathway;
}

// Update ray source point
void updateRaySource(Ray* rayPtr, Vector3D newSource) {
    rayPtr->sourcePoint = newSource;
}

// Update ray direction
void updateRayPathway(Ray* rayPtr, Vector3D newPathway) {
    rayPtr->pathway = newPathway;
}

#endif