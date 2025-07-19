#ifndef SPATIAL_PRIMITIVES_HEADER
#define SPATIAL_PRIMITIVES_HEADER

#include "vector3D.h"
#include "ray.h"
#include "lights.h"
#include <GL/glut.h>
#include <vector>
#include <climits>

// Add texture support
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Texture variables
unsigned char* textureData = nullptr;
int textureWidth = 0;
int textureHeight = 0;
int textureChannels = 0;
bool useTexture = false;  // For easy switching with button 9

// Color structure for texture sampling
struct Color {
    double r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(double red, double green, double blue) : r(red), g(green), b(blue) {}
};

// Texture sampling function - follows specification exactly
Color sampleTexture(double u, double v) {
    if (!textureData || textureWidth <= 0 || textureHeight <= 0) {
        return Color(0.5, 0.5, 0.5); // Gray fallback
    }

    // Clamp u and v to [0,1] - specification requirement
    u = std::max(0.0, std::min(1.0, u));
    v = std::max(0.0, std::min(1.0, v));

    // Normalized -> pixel coords - specification implementation
    int pixel_x = (int)(u * (textureWidth - 1));
    int pixel_y = (int)((1.0 - v) * (textureHeight - 1)); // Flip Y - specification

    // Safety clamp - specification requirement
    pixel_x = std::max(0, std::min(textureWidth - 1, pixel_x));
    pixel_y = std::max(0, std::min(textureHeight - 1, pixel_y));

    // Compute array index - specification implementation
    int index = (pixel_y * textureWidth + pixel_x) * textureChannels;
    int max_index = textureWidth * textureHeight * textureChannels;
    
    if (index < 0 || index + 2 >= max_index) {
        return Color(1.0, 0.0, 1.0); // Magenta = error - specification
    }

    Color color;
    color.r = textureData[index] / 255.0;

    if (textureChannels >= 2) {
        color.g = textureData[index + 1] / 255.0;
    } else {
        color.g = color.r; // Grayscale
    }

    if (textureChannels >= 3) {
        color.b = textureData[index + 2] / 255.0;
    } else {
        color.b = color.r; // Grayscale
    }

    return color;
}

// Load texture function
bool loadTexture(const char* filename) {
    if (textureData) {
        stbi_image_free(textureData);
        textureData = nullptr;
    }
    
    textureData = stbi_load(filename, &textureWidth, &textureHeight, &textureChannels, 0);
    
    if (!textureData) {
        printf("Failed to load texture: %s\n", filename);
        return false;
    }
    
    printf("Loaded texture: %dx%d with %d channels\n", textureWidth, textureHeight, textureChannels);
    return true;
}

// Toggle texture mode function - for easy switching as per specification
void toggleTextureMode() {
    useTexture = !useTexture;
    printf("Texture mode: %s\n", useTexture ? "ON (using texture)" : "OFF (using checkerboard)");
}

// Primary spatial primitive structure with distinctive field names
struct Object {
    Vector3D reference_point;     // Core reference coordinates
    double height, width, length;    // Shape-defining parameters
    double color[3];       // RGB material coloration
    double coEfficients[4];      // PHONG illumination coefficients [ambient, diffuse, specular, reflection]
    int shine;         // Material shininess factor
    int type;        // Shape classification identifier
    // Triangle vertex storage for proper reflection calculation
    Vector3D vertexA, vertexB, vertexC; // For triangle primitives only
};

// Primitive category constants
#define SPHERE 1
#define TRIANGLE 2
#define GENERAL 3
#define FLOOR 4

// Global primitive collections with distinctive names
extern std::vector<Object*> objects;
extern std::vector<PointLight> pointLights;
extern std::vector<SpotLight> spotLights;
extern int recursion_level;

// Initialize spatial primitive with default values
Object constructDefaultPrimitive() {
    Object primitive;
    primitive.reference_point = makeVector(0, 0, 0);
    primitive.height = 0.0;
    primitive.width = 0.0;
    primitive.length = 0.0;
    for(int i = 0; i < 3; i++) {
        primitive.color[i] = 1.0;  // Default white
    }
    for(int i = 0; i < 4; i++) {
        primitive.coEfficients[i] = 0.25;  // Default lighting
    }
    primitive.shine = 1;
    primitive.type = 0;
    return primitive;
}

// Color utility functions with enhanced intensity
void restrictColorBounds(double* colorArray) {
    for(int i = 0; i < 3; i++) {
        if(colorArray[i] > 1.5) colorArray[i] = 1.5;  // Allow higher intensity
        if(colorArray[i] < 0.0) colorArray[i] = 0.0;
    }
}

// Spherical primitive factory function
Object* createSphericalPrimitive(Vector3D centerPos, double sphereRadius) {
    Object* sphere = new Object();
    *sphere = constructDefaultPrimitive();
    sphere->type = SPHERE;
    sphere->reference_point = centerPos;
    sphere->length = sphereRadius;
    return sphere;
}

// Configure primitive material tint
void assignPrimitiveTint(Object* primitive, double colorValues[3]) {
    for(int i = 0; i < 3; i++) {
        primitive->color[i] = colorValues[i];
    }
}

// Configure primitive lighting model
void assignPrimitiveLighting(Object* primitive, double lightCoeffs[4]) {
    for(int i = 0; i < 4; i++) {
        primitive->coEfficients[i] = lightCoeffs[i];
    }
}

// Configure primitive surface roughness
void assignPrimitiveShininess(Object* primitive, int shineValue) {
    primitive->shine = shineValue;
}

// Spherical intersection computation with improved precision
double computeSphericalIntersection(Object* spherePrimitive, Ray incidentRay) {
    Vector3D rayToCenter = subtract(incidentRay.sourcePoint, spherePrimitive->reference_point);
    double rayMagnitude = dot(incidentRay.pathway, incidentRay.pathway);
    double rayDotOffset = 2.0 * dot(rayToCenter, incidentRay.pathway);
    double offsetMagnitude = dot(rayToCenter, rayToCenter) - spherePrimitive->length * spherePrimitive->length;
    
    double discriminant = rayDotOffset * rayDotOffset - 4 * rayMagnitude * offsetMagnitude;
    if(discriminant < 0) return -1;
    
    double sqrtDiscriminant = sqrt(discriminant);
    double t1 = (-rayDotOffset - sqrtDiscriminant) / (2 * rayMagnitude);
    double t2 = (-rayDotOffset + sqrtDiscriminant) / (2 * rayMagnitude);
    
    if(t1 > 1e-6) return t1;
    if(t2 > 1e-6) return t2;
    return -1;
}

// Spherical surface normal computation
Vector3D computeSphericalNormal(Object* spherePrimitive, Ray incidentRay, Vector3D intersectionPoint) {
    Vector3D surfaceNormal = subtract(intersectionPoint, spherePrimitive->reference_point);
    normalize(surfaceNormal);
    if(dot(incidentRay.pathway, surfaceNormal) > 0) {
        surfaceNormal = multiply(surfaceNormal, -1);
    }
    return surfaceNormal;
}

// Render spherical primitive
void renderSphericalPrimitive(Object* spherePrimitive) {
    glPushMatrix();
        glTranslatef(spherePrimitive->reference_point.x, spherePrimitive->reference_point.y, spherePrimitive->reference_point.z);
        glColor3f(spherePrimitive->color[0], spherePrimitive->color[1], spherePrimitive->color[2]);
        glutSolidSphere(spherePrimitive->length, 50, 50);
    glPopMatrix();
}

// Triangle primitive factory function
Object* constructTriangularPrimitive(Vector3D vertexA, Vector3D vertexB, Vector3D vertexC) {
    Object* triangle = new Object();
    *triangle = constructDefaultPrimitive();
    triangle->type = TRIANGLE;
    triangle->reference_point = vertexA; // Store first vertex as reference
    // Store all vertices for proper intersection and reflection calculation
    triangle->vertexA = vertexA;
    triangle->vertexB = vertexB;
    triangle->vertexC = vertexC;
    return triangle;
}

// Triangle primitive factory function
Object* createTriangularPrimitive(Vector3D vertexA, Vector3D vertexB, Vector3D vertexC) {
    Object* triangle = new Object();
    *triangle = constructDefaultPrimitive();
    triangle->type = TRIANGLE;
    triangle->reference_point = vertexA; // Store first vertex as reference
    // Store all vertices for proper intersection and reflection calculation
    triangle->vertexA = vertexA;
    triangle->vertexB = vertexB;
    triangle->vertexC = vertexC;
    return triangle;
}

// Determinant calculation helper
double computeDeterminant(Vector3D v1, Vector3D v2, Vector3D v3) {
    return v1.x * (v2.y * v3.z - v2.z * v3.y) - 
           v1.y * (v2.x * v3.z - v2.z * v3.x) + 
           v1.z * (v2.x * v3.y - v2.y * v3.x);
}

// Triangular intersection using barycentric coordinates
double computeTriangularIntersection(Object* trianglePrimitive, Ray incidentRay, Vector3D vertA, Vector3D vertB, Vector3D vertC) {
    Vector3D edgeAB = subtract(vertA, vertB);
    Vector3D edgeAC = subtract(vertA, vertC);
    Vector3D rayToA = subtract(vertA, incidentRay.sourcePoint);
    
    double detA = computeDeterminant(edgeAB, edgeAC, incidentRay.pathway);
    if(abs(detA) < 1e-6) return -1;
    
    double detT = computeDeterminant(edgeAB, edgeAC, rayToA);
    double t = detT / detA;
    if(t < 1e-6) return -1;
    
    double detBeta = computeDeterminant(rayToA, edgeAC, incidentRay.pathway);
    double beta = detBeta / detA;
    if(beta < 0 || beta > 1) return -1;
    
    double detGamma = computeDeterminant(edgeAB, rayToA, incidentRay.pathway);
    double gamma = detGamma / detA;
    if(gamma < 0 || gamma > 1 || (beta + gamma) > 1) return -1;
    
    return t;
}

// Triangular surface normal computation
Vector3D computeTriangularNormal(Object* trianglePrimitive, Ray incidentRay, Vector3D vertA, Vector3D vertB, Vector3D vertC) {
    Vector3D edgeBA = subtract(vertB, vertA);
    Vector3D edgeCA = subtract(vertC, vertA);
    Vector3D surfaceNormal = cross(edgeBA, edgeCA);
    normalize(surfaceNormal);
    if(dot(incidentRay.pathway, surfaceNormal) > 0) {
        surfaceNormal = multiply(surfaceNormal, -1);
    }
    return surfaceNormal;
}

// Render triangular primitive
void renderTriangularPrimitive(Object* trianglePrimitive, Vector3D vertA, Vector3D vertB, Vector3D vertC) {
    glPushMatrix();
        glBegin(GL_TRIANGLES);
            glColor3f(trianglePrimitive->color[0], trianglePrimitive->color[1], trianglePrimitive->color[2]);
            glVertex3f(vertA.x, vertA.y, vertA.z);
            glVertex3f(vertB.x, vertB.y, vertB.z);
            glVertex3f(vertC.x, vertC.y, vertC.z);
        glEnd();
    glPopMatrix();
}

// Floor primitive factory function
Object* createPlanarSurface(double planeExtent, double checkerSize) {
    Object* floor = new Object();
    *floor = constructDefaultPrimitive();
    floor->type = FLOOR;
    floor->width = planeExtent;
    floor->length = checkerSize;
    return floor;
}

// Floor intersection computation
double computePlanarIntersection(Object* floorPrimitive, Ray incidentRay, double floorExtent) {
    if(abs(incidentRay.pathway.z) < 1e-6) return -1; // Ray parallel to floor
    
    double t = -incidentRay.sourcePoint.z / incidentRay.pathway.z;
    if(t < 1e-6) return -1;
    
    Vector3D intersectionPoint = add(incidentRay.sourcePoint, multiply(incidentRay.pathway, t));
    if(abs(intersectionPoint.x) > floorExtent || abs(intersectionPoint.y) > floorExtent) {
        return -1;
    }
    
    return t;
}

// Floor surface normal computation - perfectly flat floor
Vector3D computePlanarNormal(Object* floorPrimitive, Ray incidentRay) {
    // Floor normal is always perfectly vertical (0, 0, 1) for flat appearance
    Vector3D floorNormal = makeVector(0, 0, 1);
    
    // Always ensure normal points away from incoming ray for proper lighting
    if(dot(incidentRay.pathway, floorNormal) > 0) {
        floorNormal = makeVector(0, 0, -1);
    }
    
    return floorNormal;
}

// Enhanced floor color computation with proper texture mapping and reflection support
void computeFloorColorPattern(Object* floorPrimitive, Vector3D intersectionPoint, double* outputColor) {
    if (useTexture && textureData) {
        // Proper texture coordinate calculation - one texture per checker square
        double checkerSize = floorPrimitive->length; // This is 20 from scene setup
        if(checkerSize <= 0) checkerSize = 20.0;
        
        // Find which checker square we're in
        int checkX = (int)floor(intersectionPoint.x / checkerSize);
        int checkY = (int)floor(intersectionPoint.y / checkerSize);
        
        // Calculate position within this checker square
        double localX = intersectionPoint.x - checkX * checkerSize;
        double localY = intersectionPoint.y - checkY * checkerSize;
        
        // Handle negative coordinates properly
        if (localX < 0) localX += checkerSize;
        if (localY < 0) localY += checkerSize;
        
        // Convert to texture coordinates (0 to 1) within this checker square
        double u = localX / checkerSize;
        double v = localY / checkerSize;
        
        // Ensure u and v are in [0,1] range
        u = std::max(0.0, std::min(1.0, u));
        v = std::max(0.0, std::min(1.0, v));
        
        // Sample texture color using specification-compliant function
        Color texColor = sampleTexture(u, v);
        outputColor[0] = texColor.r;
        outputColor[1] = texColor.g;
        outputColor[2] = texColor.b;
        
        // Apply natural enhancement for better visibility
        for(int i = 0; i < 3; i++) {
            outputColor[i] = std::min(0.9, outputColor[i] * 1.02);
        }
    } else {
        // Original checkerboard pattern with reflection support
        double checkerSize = floorPrimitive->length;
        if(checkerSize <= 0) checkerSize = 20.0; // Fallback value
        
        int checkX = (int)floor(intersectionPoint.x / checkerSize);
        int checkY = (int)floor(intersectionPoint.y / checkerSize);
        
        // Proper checkerboard pattern that handles negative coordinates correctly
        // Use modulo arithmetic to ensure consistent alternating pattern
        int patternValue = (checkX + checkY) % 2;
        if(patternValue < 0) patternValue += 2; // Handle negative modulo
        
        if(patternValue == 0) {
            // White squares - moderate reflectivity for natural appearance
            outputColor[0] = outputColor[1] = outputColor[2] = 0.7;
        } else {
            // Dark squares - still reflective but less bright
            outputColor[0] = outputColor[1] = outputColor[2] = 0.15;
        }
    }
}

// Render floor primitive
void renderPlanarSurface(Object* floorPrimitive, double floorExtent) {
    double checkerSize = floorPrimitive->length;
    int numCheckers = (int)(2 * floorExtent / checkerSize);
    
    glPushMatrix();
        for(int i = 0; i < numCheckers; i++) {
            for(int j = 0; j < numCheckers; j++) {
                double x1 = -floorExtent + i * checkerSize;
                double y1 = -floorExtent + j * checkerSize;
                double x2 = x1 + checkerSize;
                double y2 = y1 + checkerSize;
                
                if((i + j) % 2 == 0) {
                    glColor3f(1.0, 1.0, 1.0); // White
                } else {
                    glColor3f(0.0, 0.0, 0.0); // Black
                }
                
                glBegin(GL_QUADS);
                    glVertex3f(x1, y1, 0);
                    glVertex3f(x2, y1, 0);
                    glVertex3f(x2, y2, 0);
                    glVertex3f(x1, y2, 0);
                glEnd();
            }
        }
    glPopMatrix();
}

// General quadric surface functions
Object* createQuadricSurface(double coefficients[10], Vector3D referencePoint, double dimX, double dimY, double dimZ) {
    Object* quadric = new Object();
    *quadric = constructDefaultPrimitive();
    quadric->type = GENERAL;
    quadric->reference_point = referencePoint;
    quadric->length = dimX;
    quadric->width = dimY;
    quadric->height = dimZ;
    return quadric;
}

// Quadric intersection computation (simplified ellipsoid)
double computeQuadricIntersection(Object* quadricPrimitive, Ray incidentRay) {
    Vector3D rayOriginRelative = subtract(incidentRay.sourcePoint, quadricPrimitive->reference_point);
    
    // Simplified ellipsoid equation: (x/a)² + (y/b)² + (z/c)² = 1
    double a = (quadricPrimitive->length > 0) ? quadricPrimitive->length : 10.0;
    double b = (quadricPrimitive->width > 0) ? quadricPrimitive->width : 10.0;
    double c = (quadricPrimitive->height > 0) ? quadricPrimitive->height : 10.0;
    
    double A = (incidentRay.pathway.x * incidentRay.pathway.x) / (a * a) + 
               (incidentRay.pathway.y * incidentRay.pathway.y) / (b * b) + 
               (incidentRay.pathway.z * incidentRay.pathway.z) / (c * c);
    
    double B = 2.0 * ((rayOriginRelative.x * incidentRay.pathway.x) / (a * a) + 
                     (rayOriginRelative.y * incidentRay.pathway.y) / (b * b) + 
                     (rayOriginRelative.z * incidentRay.pathway.z) / (c * c));
    
    double C = (rayOriginRelative.x * rayOriginRelative.x) / (a * a) + 
               (rayOriginRelative.y * rayOriginRelative.y) / (b * b) + 
               (rayOriginRelative.z * rayOriginRelative.z) / (c * c) - 1.0;
    
    double discriminant = B * B - 4 * A * C;
    if(discriminant < 0) return -1;
    
    double t1 = (-B - sqrt(discriminant)) / (2 * A);
    double t2 = (-B + sqrt(discriminant)) / (2 * A);
    
    if(t1 > 1e-6) return t1;
    if(t2 > 1e-6) return t2;
    return -1;
}

// Quadric surface normal computation
Vector3D computeQuadricNormal(Object* quadricPrimitive, Vector3D intersectionPoint) {
    Vector3D relativePos = subtract(intersectionPoint, quadricPrimitive->reference_point);
    
    double a = (quadricPrimitive->length > 0) ? quadricPrimitive->length : 10.0;
    double b = (quadricPrimitive->width > 0) ? quadricPrimitive->width : 10.0;
    double c = (quadricPrimitive->height > 0) ? quadricPrimitive->height : 10.0;
    
    Vector3D surfaceNormal = makeVector(
        2 * relativePos.x / (a * a),
        2 * relativePos.y / (b * b),
        2 * relativePos.z / (c * c)
    );
    normalize(surfaceNormal);
    return surfaceNormal;
}

// Render quadric surface (disabled for raytracing-only objects)
void renderQuadricSurface(Object* quadricPrimitive) {
    // According to assignment PDF, general quadric surfaces should NOT be rendered in OpenGL
    // They appear only in raytraced output
}

// IMPROVED REFLECTION AND LIGHTING SYSTEM
// Main raytracing computation with enhanced reflection handling
double intersect(Object* primitive, Ray incomingRay, double* outputColor, int currentLevel, 
                                      Vector3D vertA = makeVector(0,0,0), Vector3D vertB = makeVector(0,0,0), Vector3D vertC = makeVector(0,0,0),
                                      double floorExtent = 1000) {
    
    double intersectionDistance = -1;
    
    // For triangles, use vertices from Object structure if available, otherwise use parameters
    Vector3D actualVertA = vertA, actualVertB = vertB, actualVertC = vertC;
    if(primitive->type == TRIANGLE) {
        // Check if vertices are stored in Object structure (non-zero vertices)
        if(primitive->vertexA.x != 0 || primitive->vertexA.y != 0 || primitive->vertexA.z != 0) {
            actualVertA = primitive->vertexA;
            actualVertB = primitive->vertexB;
            actualVertC = primitive->vertexC;
        }
    }
    
    // Calculate intersection distance based on primitive category
    switch(primitive->type) {
        case SPHERE:
            intersectionDistance = computeSphericalIntersection(primitive, incomingRay);
            break;
        case TRIANGLE:
            intersectionDistance = computeTriangularIntersection(primitive, incomingRay, actualVertA, actualVertB, actualVertC);
            break;
        case GENERAL:
            intersectionDistance = computeQuadricIntersection(primitive, incomingRay);
            break;
        case FLOOR:
            intersectionDistance = computePlanarIntersection(primitive, incomingRay, floorExtent);
            break;
        default:
            return -1;
    }
    
    if(intersectionDistance < 0 || currentLevel == 0) {
        return intersectionDistance;
    }
    
    // Calculate intersection point
    Vector3D intersectionLocation = add(incomingRay.sourcePoint, multiply(incomingRay.pathway, intersectionDistance));
    
    // Get surface color at intersection point
    double surfaceColorAtPoint[3];
    if(primitive->type == FLOOR) {
        computeFloorColorPattern(primitive, intersectionLocation, surfaceColorAtPoint);
    } else {
        for(int i = 0; i < 3; i++) {
            surfaceColorAtPoint[i] = primitive->color[i];
        }
    }
    
    // Initialize with enhanced ambient lighting
    for(int i = 0; i < 3; i++) {
        outputColor[i] = surfaceColorAtPoint[i] * primitive->coEfficients[0] * 1.2;  // Boost ambient
    }
    restrictColorBounds(outputColor);
    
    // Calculate surface normal
    Vector3D surfaceNormal;
    switch(primitive->type) {
        case SPHERE:
            surfaceNormal = computeSphericalNormal(primitive, incomingRay, intersectionLocation);
            break;
        case TRIANGLE:
            surfaceNormal = computeTriangularNormal(primitive, incomingRay, actualVertA, actualVertB, actualVertC);
            break;
        case GENERAL:
            surfaceNormal = computeQuadricNormal(primitive, intersectionLocation);
            break;
        case FLOOR:
            surfaceNormal = computePlanarNormal(primitive, incomingRay);
            break;
    }
    
    // Ensure normal points toward the ray origin (Sadat's improvement)
    if(dot(incomingRay.pathway, surfaceNormal) > 0) {
        surfaceNormal = multiply(surfaceNormal, -1);
    }
    
    // Process point light sources with improved shadow calculation (Sadat's approach)
    for(int i = 0; i < pointLights.size(); i++) {
        Vector3D lightDirection = subtract(pointLights[i].coordinates, intersectionLocation);
        double lightDistance = sqrt(lightDirection.x * lightDirection.x + lightDirection.y * lightDirection.y + lightDirection.z * lightDirection.z);
        
        if(lightDistance < 1e-6) continue; // Light source is at intersection point
        
        normalize(lightDirection);
          // Create shadow ray from intersection point toward light (Sadat's method)
        Ray shadowRay = buildRayFromVectors(intersectionLocation, lightDirection);
        standardizeRayPathway(&shadowRay);
        
        // Move shadow ray origin slightly towards the light to avoid self-intersection
        shadowRay.sourcePoint = add(shadowRay.sourcePoint, multiply(lightDirection, 1e-4));

        bool inShadow = false;
        
        // Check for shadow-casting objects with improved epsilon handling
        for(int k = 0; k < objects.size(); k++) {
            if(objects[k] == primitive) continue; // Skip self-intersection
            
            double tempColorBuffer[3] = {0, 0, 0};
            double shadowDistance = -1;
            
            // Use stored vertices for triangle intersection in shadows
            if(objects[k]->type == TRIANGLE) {
                shadowDistance = intersect(objects[k], shadowRay, tempColorBuffer, 0, 
                                         objects[k]->vertexA, objects[k]->vertexB, objects[k]->vertexC);
            } else {
                shadowDistance = intersect(objects[k], shadowRay, tempColorBuffer, 0);
            }
            
            if(shadowDistance > 1e-4 && shadowDistance + 1e-4 < lightDistance) {
                inShadow = true;
                break;
            }
        }
        
        // Apply lighting if not in shadow
        if(!inShadow) {
            // Diffuse (Lambert) lighting - improved calculation
            double lambertian = dot(surfaceNormal, lightDirection);
            lambertian = max(lambertian, 0.0);
            
            if(lambertian > 1e-6) {
                for(int j = 0; j < 3; j++) {
                    outputColor[j] += pointLights[i].colorChannels[j] * primitive->coEfficients[1] * 
                                    lambertian * surfaceColorAtPoint[j] * 1.0;  // Normal diffuse
                }
                restrictColorBounds(outputColor);
                
                // Specular (Phong) lighting with improved reflection calculation (Sadat's method)
                Vector3D reflectedLight = subtract(lightDirection, multiply(surfaceNormal, 2 * dot(lightDirection, surfaceNormal)));
                normalize(reflectedLight);
                Vector3D viewDirection = multiply(incomingRay.pathway, -1);
                normalize(viewDirection);
                
                double specularFactor = pow(max(dot(reflectedLight, viewDirection), 0.0), primitive->shine);
                for(int j = 0; j < 3; j++) {
                    outputColor[j] += pointLights[i].colorChannels[j] * primitive->coEfficients[2] * 
                                    specularFactor * surfaceColorAtPoint[j] * 1.0;  // Moderate specular
                }
                restrictColorBounds(outputColor);
            }
        }
    }
    
    // Process spot light sources with improved angle calculation (Sadat's approach)
    for(int i = 0; i < spotLights.size(); i++) {
        Vector3D lightDirection = subtract(spotLights[i].coordinates, intersectionLocation);
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
          // Create shadow ray
        Ray shadowRay = buildRayFromVectors(intersectionLocation, lightDirection);
        standardizeRayPathway(&shadowRay);
        
        // Move shadow ray origin slightly towards the light to avoid self-intersection
        shadowRay.sourcePoint = add(shadowRay.sourcePoint, multiply(lightDirection, 1e-4));

        bool inShadow = false;
        
        // Check for shadows with improved epsilon handling
        for(int k = 0; k < objects.size(); k++) {
            if(objects[k] == primitive) continue;
            
            double tempColorBuffer[3] = {0, 0, 0};
            double shadowDistance = intersect(objects[k], shadowRay, tempColorBuffer, 0);
            if(shadowDistance > 1e-6 && shadowDistance + 1e-6 < lightDistance) {
                inShadow = true;
                break;
            }
        }
        
        // Apply spot lighting if not in shadow
        if(!inShadow) {
            double lambertian = dot(surfaceNormal, lightDirection);
            lambertian = max(lambertian, 0.0);
            
            if(lambertian > 1e-6) {
                // Spot light intensity falloff (Sadat's method)
                double beta = angle * 3.1416 / 180.0;
                double spotFalloff = 2.0;
                double spotIntensity = pow(cos(beta), spotFalloff);
                
                for(int j = 0; j < 3; j++) {
                    outputColor[j] += spotLights[i].colorChannels[j] * primitive->coEfficients[1] * 
                                    lambertian * surfaceColorAtPoint[j] * spotIntensity * 1.0;  // Normal spot diffuse
                }
                restrictColorBounds(outputColor);
                
                Vector3D reflectedLight = subtract(lightDirection, multiply(surfaceNormal, 2 * dot(lightDirection, surfaceNormal)));
                normalize(reflectedLight);
                Vector3D viewDirection = multiply(incomingRay.pathway, -1);
                normalize(viewDirection);
                
                double specularFactor = pow(max(dot(reflectedLight, viewDirection), 0.0), primitive->shine);
                for(int j = 0; j < 3; j++) {
                    outputColor[j] += spotLights[i].colorChannels[j] * primitive->coEfficients[2] * 
                                    specularFactor * surfaceColorAtPoint[j] * spotIntensity * 1.0;  // Moderate spot specular
                }
                restrictColorBounds(outputColor);
            }
        }
    }
    
    // Handle recursive reflections with improved reflection calculation (Sadat's method)
    if(currentLevel < recursion_level && primitive->coEfficients[3] > 0.01) {
        Vector3D reflectedDirection = subtract(incomingRay.pathway, multiply(surfaceNormal, 2 * dot(incomingRay.pathway, surfaceNormal)));
        normalize(reflectedDirection);
        
        Ray reflectionRay = buildRayFromVectors(intersectionLocation, reflectedDirection);
        // Better epsilon offset to avoid self-intersection - always use surface normal
        reflectionRay.sourcePoint = add(reflectionRay.sourcePoint, multiply(surfaceNormal, 1e-3));
        standardizeRayPathway(&reflectionRay);
        
        double reflectionColorBuffer[3] = {0, 0, 0};
        double nearestReflectionDistance = 1e9;
        int nearestReflector = -1;
        
        // Find nearest reflecting object with improved distance checking
        for(int k = 0; k < objects.size(); k++) {
            if(objects[k] == primitive) continue; // Skip self-intersection
            
            double tempColorBuffer[3] = {0, 0, 0};
            double reflectionDistance = -1;
            
            // Use stored vertices for triangle intersection in reflections
            if(objects[k]->type == TRIANGLE) {
                reflectionDistance = intersect(objects[k], reflectionRay, tempColorBuffer, 0, 
                                             objects[k]->vertexA, objects[k]->vertexB, objects[k]->vertexC);
            } else {
                reflectionDistance = intersect(objects[k], reflectionRay, tempColorBuffer, 0);
            }
            
            if(reflectionDistance > 1e-4 && reflectionDistance < nearestReflectionDistance) {
                nearestReflectionDistance = reflectionDistance;
                nearestReflector = k;
            }
        }
        
        // Apply reflection contribution with enhanced floor reflection blending
        if(nearestReflector >= 0) {
            // Call intersect with proper vertices for triangles
            if(objects[nearestReflector]->type == TRIANGLE) {
                intersect(objects[nearestReflector], reflectionRay, reflectionColorBuffer, currentLevel + 1,
                         objects[nearestReflector]->vertexA, objects[nearestReflector]->vertexB, objects[nearestReflector]->vertexC);
            } else {
                intersect(objects[nearestReflector], reflectionRay, reflectionColorBuffer, currentLevel + 1);
            }
            
            // Balanced reflection contribution for more natural appearance
            double reflectionCoeff = primitive->coEfficients[3];
            if(primitive->type == FLOOR) {
                // Natural floor reflection - visible but not overwhelming
                reflectionCoeff = std::min(0.5, reflectionCoeff * 0.9); // Keep reflection moderate
                
                // Apply reflection with proper blending for floor
                for(int i = 0; i < 3; i++) {
                    // Balanced reflection blending for natural appearance
                    outputColor[i] = outputColor[i] * (1.0 - reflectionCoeff) + 
                                   reflectionColorBuffer[i] * reflectionCoeff;
                }
            } else {
                // Moderate reflection for other objects
                for(int i = 0; i < 3; i++) {
                    outputColor[i] += reflectionColorBuffer[i] * reflectionCoeff * 0.7;  // Reduced boost
                }
            }
            restrictColorBounds(outputColor);
        }
    }
    
    return intersectionDistance;
}

#endif
