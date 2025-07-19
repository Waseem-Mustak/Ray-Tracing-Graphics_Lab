// Advanced Ray Tracing Implementation
// Compilation: g++ 2005091_main.cpp -o demo.exe -lfreeglut -lglew32 -lopengl32 -lglu32

#include <iostream>
#include <GL/glut.h>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <cmath>

#include "bitmap_image.hpp"
#include "2005091_classes.hpp"

using namespace std;

// === FORWARD DECLARATIONS ===
void processSphereObject(ifstream& file);
void processTriangleObject(ifstream& file);
void processQuadricObject(ifstream& file);
void loadLightingSources(ifstream& file);
void addFloorToScene();

// === GLOBAL CONFIGURATION ===
struct RenderConfig {
    double screenWidth = 500.0;
    double screenHeight = 500.0;
    double fieldOfView = 80.0;
    double imageResolution = 0.0;
    int captureCounter = 1;
};

struct SceneData {
    vector<sphere*> sphericalObjects;
    vector<triangle*> triangularObjects;
    vector<general*> quadricObjects;
    vector<Floor*> floorObjects;
};

// Global instances
camera mainCamera;
RenderConfig renderSettings;
SceneData sceneObjects;

// === INITIALIZATION FUNCTIONS ===
void initializeGraphicsSystem() {
    // Setup initial camera position and orientation
    vector3D cameraPosition = makeVector(115, 115, 50);
    vector3D upDirection = makeVector(0, 0, 1);
    vector3D lookDirection = makeVector(-1/sqrt(2), -1/sqrt(2), 0);

    mainCamera.setCoords(cameraPosition, lookDirection, upDirection);

    // Configure OpenGL settings
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(renderSettings.fieldOfView, 1.0, 1.0, 1000.0);
}

void parseSceneFile() {
    ifstream sceneFile("scene.txt");
    if (!sceneFile.is_open()) {
        cerr << "Error: Cannot open scene.txt file!" << endl;
        return;
    }

    sceneFile >> recursion_level;
    sceneFile >> renderSettings.imageResolution;
    
    int totalObjects;
    sceneFile >> totalObjects;

    // Initialize texture loading
    bool textureLoaded = loadTexture("texture.jpg") || 
                        loadTexture("texture.jpg");
    if (!textureLoaded) {
        printf("Info: Using procedural patterns instead of texture\n");
    }

    // Parse objects from scene file
    for (int objIndex = 0; objIndex < totalObjects; objIndex++) {
        string objectCategory;
        sceneFile >> objectCategory;
        
        if (objectCategory == "sphere") {
            processSphereObject(sceneFile);
        }
        else if (objectCategory == "triangle") {
            processTriangleObject(sceneFile);
        }
        else if (objectCategory == "general") {
            processQuadricObject(sceneFile);
        }
        else {
            cout << "Warning: Unknown object type: " << objectCategory << endl;
        }
    }

    // Load lighting information
    loadLightingSources(sceneFile);
    
    // Add floor to scene
    addFloorToScene();
    
    sceneFile.close();
}

void drawAxes(){
    glBegin(GL_LINES);{ 
        glColor3f(1, 0, 0);
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);
        glColor3f(0, 1, 0);
        glVertex3f(0, 100, 0);
        glVertex3f(0, -100, 0);
        glColor3f(0, 0, 1);
        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    }glEnd();
}

// === RENDERING FUNCTIONS ===
void renderCoordinateAxes() {
    glBegin(GL_LINES);
    {
        // X-axis (Red)
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);
        
        // Y-axis (Green)
        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(0, 100, 0);
        glVertex3f(0, -100, 0);
        
        // Z-axis (Blue)
        glColor3f(0.0, 0.0, 1.0);
        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    }
    glEnd();
}

struct IntersectionResult {
    double distance;
    int objectIndex;
    int objectType; // 0=sphere, 1=triangle, 2=quadric, 3=floor
};

IntersectionResult findNearestIntersection(const Ray& ray) {
    IntersectionResult result = {-1.0, -1, -1};
    double minDistance = INFINITY;
    double tempColor[3] = {0, 0, 0};
    
    // Check spheres
    for (int i = 0; i < sceneObjects.sphericalObjects.size(); i++) {
        double t = sceneObjects.sphericalObjects[i]->intersect(ray, tempColor, 0);
        if (t > 0 && t < minDistance) {
            minDistance = t;
            result = {t, i, 0};
        }
    }
    
    // Check triangles
    for (int i = 0; i < sceneObjects.triangularObjects.size(); i++) {
        double t = sceneObjects.triangularObjects[i]->intersect(ray, tempColor, 0);
        if (t > 0 && t < minDistance) {
            minDistance = t;
            result = {t, i, 1};
        }
    }
    
    // Check quadrics
    for (int i = 0; i < sceneObjects.quadricObjects.size(); i++) {
        double t = sceneObjects.quadricObjects[i]->intersect(ray, tempColor, 0);
        if (t > 0 && t < minDistance) {
            minDistance = t;
            result = {t, i, 2};
        }
    }
    
    // Check floor
    for (int i = 0; i < sceneObjects.floorObjects.size(); i++) {
        double t = sceneObjects.floorObjects[i]->intersect(ray, tempColor, 0);
        if (t > 0 && t < minDistance) {
            minDistance = t;
            result = {t, i, 3};
        }
    }
    
    return result;
}

void calculatePixelColor(const Ray& ray, const IntersectionResult& intersection, double* pixelColor) {
    switch (intersection.objectType) {
        case 0: // Sphere
            sceneObjects.sphericalObjects[intersection.objectIndex]->intersect(ray, pixelColor, 1);
            break;
        case 1: // Triangle
            sceneObjects.triangularObjects[intersection.objectIndex]->intersect(ray, pixelColor, 1);
            break;
        case 2: // Quadric
            sceneObjects.quadricObjects[intersection.objectIndex]->intersect(ray, pixelColor, 1);
            break;
        case 3: // Floor
            sceneObjects.floorObjects[intersection.objectIndex]->intersect(ray, pixelColor, 1);
            break;
    }
    
    // Clamp color values
    for (int i = 0; i < 3; i++) {
        pixelColor[i] = max(0.0, min(1.0, pixelColor[i]));
    }
}

void performRayTracing() {
    cout << "Camera position: " << mainCamera.pos.x << ", " << mainCamera.pos.y << ", " << mainCamera.pos.z << endl;
    
    bitmap_image outputImage(renderSettings.imageResolution, renderSettings.imageResolution);
    outputImage.set_all_channels(0, 0, 0);
    
    // Calculate view plane parameters
    double viewPlaneDistance = (renderSettings.screenHeight / 2.0) / 
                              tan((renderSettings.fieldOfView / 2.0) * (PI / 180.0));
    
    vector3D topLeftCorner = add(
        add(add(mainCamera.pos, multiply(mainCamera.l, viewPlaneDistance)),
            multiply(mainCamera.r, -(renderSettings.screenWidth / 2.0))),
        multiply(mainCamera.u, (renderSettings.screenHeight / 2.0))
    );
    
    double pixelStepX = renderSettings.screenWidth / renderSettings.imageResolution;
    double pixelStepY = renderSettings.screenHeight / renderSettings.imageResolution;
    
    topLeftCorner = subtract(
        add(topLeftCorner, multiply(mainCamera.r, (pixelStepX / 2.0))),
        multiply(mainCamera.u, (pixelStepY / 2.0))
    );
    
    // Ray tracing loop
    for (int row = 0; row < renderSettings.imageResolution; row++) {
        for (int col = 0; col < renderSettings.imageResolution; col++) {
            vector3D pixelPosition = subtract(
                add(topLeftCorner, multiply(mainCamera.r, pixelStepX * col)),
                multiply(mainCamera.u, pixelStepY * row)
            );
            
            vector3D rayDirection = subtract(pixelPosition, mainCamera.pos);
            Ray currentRay = buildRayFromVectors(mainCamera.pos, rayDirection);
            standardizeRayPathway(&currentRay);
            
            IntersectionResult intersection = findNearestIntersection(currentRay);
            
            if (intersection.distance > 0) {
                double pixelColor[3] = {0, 0, 0};
                calculatePixelColor(currentRay, intersection, pixelColor);
                outputImage.set_pixel(col, row, 
                                    pixelColor[0] * 255, 
                                    pixelColor[1] * 255, 
                                    pixelColor[2] * 255);
            }
        }
    }
    
    string filename = "RayTrace_" + to_string(renderSettings.captureCounter) + ".bmp";
    outputImage.save_image(filename);
    outputImage.clear();
    renderSettings.captureCounter++;
    cout << "Image saved as: " << filename << endl;
}

// === DISPLAY AND INTERACTION ===
void renderScene() {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(mainCamera.pos.x, mainCamera.pos.y, mainCamera.pos.z,
              mainCamera.pos.x + mainCamera.l.x, mainCamera.pos.y + mainCamera.l.y, mainCamera.pos.z + mainCamera.l.z,
              mainCamera.u.x, mainCamera.u.y, mainCamera.u.z);

    // Render all objects
    for (auto& sphere : sceneObjects.sphericalObjects) {
        sphere->draw();
    }
    for (auto& triangle : sceneObjects.triangularObjects) {
        triangle->draw();
    }
    for (auto& quadric : sceneObjects.quadricObjects) {
        quadric->draw();
    }
    for (auto& floor : sceneObjects.floorObjects) {
        floor->draw();
    }

    // Render lights
    for (auto& light : pointLights) {
        visualizePointIllumination(light);
    }
    for (auto& light : spotLights) {
        visualizeSpotIllumination(light);
    }

    glutSwapBuffers();
}

void handleArrowKeys(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
            moveCameraForwardFast();
            break;
        case GLUT_KEY_DOWN:
            moveCameraBackwardFast();
            break;
        case GLUT_KEY_LEFT:
            moveCameraLeftFast();
            break;
        case GLUT_KEY_RIGHT:
            moveCameraRightFast();
            break;
        case GLUT_KEY_PAGE_UP:
            moveCameraUpGlobalFast();
            break;
        case GLUT_KEY_PAGE_DOWN:
            moveCameraDownGlobalFast();
            break;
    }
    mainCamera.updateFromGlobalCamera();
}

void handleKeyboard(unsigned char key, int x, int y) {
    switch (key) {
        case '0':
            performRayTracing();
            break;
        case '1':
            rotateCameraLeft();
            break;
        case '2':
            rotateCameraRight();
            break;
        case '3':
            rotateCameraUp();
            break;
        case '4':
            rotateCameraDown();
            break;
        case '5':
            rollCameraRight();
            break;
        case '6':
            rollCameraLeft();
            break;
        case 'w': case 'W':
            moveCameraUp();
            break;
        case 's': case 'S':
            moveCameraDown();
            break;
        case 'a': case 'A':
            moveCameraLeft();
            break;
        case 'd': case 'D':
            moveCameraRight();
            break;
        case 't': case 'T':
            toggleTextureMode();
            break;
        case 27: // ESC key
            exit(0);
            break;
    }
    mainCamera.updateFromGlobalCamera();
}

void refreshDisplay() {
    glutPostRedisplay();
}

// === OBJECT PARSING HELPERS ===
void processSphereObject(ifstream& file) {
    vector3D centerPoint;
    double sphereRadius;
    
    file >> centerPoint.x >> centerPoint.y >> centerPoint.z >> sphereRadius;
    
    sphere* newSphere = new sphere(centerPoint, sphereRadius);
    
    // Read material properties
    double surfaceColor[3], lightingCoeffs[4];
    int reflectivity;
    file >> surfaceColor[0] >> surfaceColor[1] >> surfaceColor[2];
    file >> lightingCoeffs[0] >> lightingCoeffs[1] >> lightingCoeffs[2] >> lightingCoeffs[3];
    file >> reflectivity;
    
    newSphere->setColor(surfaceColor);
    newSphere->setLight(lightingCoeffs);
    newSphere->setShine(reflectivity);
    
    sceneObjects.sphericalObjects.push_back(newSphere);
    objects.push_back(newSphere->spatialForm);
}

void processTriangleObject(ifstream& file) {
    vector3D vertex1, vertex2, vertex3;
    
    file >> vertex1.x >> vertex1.y >> vertex1.z;
    file >> vertex2.x >> vertex2.y >> vertex2.z;
    file >> vertex3.x >> vertex3.y >> vertex3.z;
    
    triangle* newTriangle = new triangle(vertex1, vertex2, vertex3);
    
    // Read material properties
    double surfaceColor[3], lightingCoeffs[4];
    int reflectivity;
    file >> surfaceColor[0] >> surfaceColor[1] >> surfaceColor[2];
    file >> lightingCoeffs[0] >> lightingCoeffs[1] >> lightingCoeffs[2] >> lightingCoeffs[3];
    file >> reflectivity;
    
    newTriangle->setColor(surfaceColor);
    newTriangle->setLight(lightingCoeffs);
    newTriangle->setShine(reflectivity);
    
    sceneObjects.triangularObjects.push_back(newTriangle);
    objects.push_back(newTriangle->spatialForm);
}

void processQuadricObject(ifstream& file) {
    struct QuadricCoeffs {
        double A, B, C, D, E, F, G, H, I, J;
    } coefficients;
    
    file >> coefficients.A >> coefficients.B >> coefficients.C >> coefficients.D >> coefficients.E
         >> coefficients.F >> coefficients.G >> coefficients.H >> coefficients.I >> coefficients.J;
    
    vector3D referencePoint;
    double boundingParams[3];
    file >> referencePoint.x >> referencePoint.y >> referencePoint.z;
    file >> boundingParams[0] >> boundingParams[1] >> boundingParams[2];
    
    double coeffArray[10] = {
        coefficients.A, coefficients.B, coefficients.C, coefficients.D, coefficients.E,
        coefficients.F, coefficients.G, coefficients.H, coefficients.I, coefficients.J
    };
    
    general* newQuadric = new general(coeffArray, referencePoint, 
                                     boundingParams[0], boundingParams[1], boundingParams[2]);
    
    // Read material properties
    double surfaceColor[3], lightingCoeffs[4];
    int reflectivity;
    file >> surfaceColor[0] >> surfaceColor[1] >> surfaceColor[2];
    file >> lightingCoeffs[0] >> lightingCoeffs[1] >> lightingCoeffs[2] >> lightingCoeffs[3];
    file >> reflectivity;
    
    newQuadric->setColor(surfaceColor);
    newQuadric->setLight(lightingCoeffs);
    newQuadric->setShine(reflectivity);
    
    sceneObjects.quadricObjects.push_back(newQuadric);
    objects.push_back(newQuadric->spatialForm);
}

void loadLightingSources(ifstream& file) {
    // Load point lights
    int pointLightCount;
    file >> pointLightCount;
    
    for (int i = 0; i < pointLightCount; i++) {
        vector3D lightPosition;
        double lightColor[3];
        file >> lightPosition.x >> lightPosition.y >> lightPosition.z;
        file >> lightColor[0] >> lightColor[1] >> lightColor[2];
        
        pointLight newPointLight = manufacturePointIllumination(lightPosition, lightColor);
        pointLights.push_back(newPointLight);
    }
    
    // Load spot lights
    int spotLightCount;
    file >> spotLightCount;
    
    for (int i = 0; i < spotLightCount; i++) {
        vector3D lightPosition, lightDirection;
        double lightColor[3], cutoffAngle;
        
        file >> lightPosition.x >> lightPosition.y >> lightPosition.z;
        file >> lightColor[0] >> lightColor[1] >> lightColor[2];
        file >> lightDirection.x >> lightDirection.y >> lightDirection.z;
        file >> cutoffAngle;
        
        spotLight newSpotLight = manufactureSpotIllumination(lightPosition, lightColor, lightDirection, cutoffAngle);
        spotLights.push_back(newSpotLight);
    }
}

void addFloorToScene() {
    Floor* groundPlane = new Floor(1000, 20);
    // Balanced lighting coefficients for realistic reflections
    // [ambient, diffuse, specular, reflection] - moderate reflection for natural appearance
    groundPlane->setLight(new double[4]{0.2, 0.6, 0.2, 0.3}); // Reduced reflection to 0.4
    groundPlane->setShine(50); // Moderate shininess for natural appearance
    
    sceneObjects.floorObjects.push_back(groundPlane);
    objects.push_back(groundPlane->spatialForm);
}

// === MEMORY MANAGEMENT ===
void cleanupResources() {
    // Clean up spherical objects
    for (auto& sphere : sceneObjects.sphericalObjects) {
        delete sphere;
    }
    sceneObjects.sphericalObjects.clear();
    
    // Clean up triangular objects
    for (auto& triangle : sceneObjects.triangularObjects) {
        delete triangle;
    }
    sceneObjects.triangularObjects.clear();
    
    // Clean up quadric objects
    for (auto& quadric : sceneObjects.quadricObjects) {
        delete quadric;
    }
    sceneObjects.quadricObjects.clear();
    
    // Clean up floor objects
    for (auto& floor : sceneObjects.floorObjects) {
        delete floor;
    }
    sceneObjects.floorObjects.clear();
    
    // Clear global collections
    objects.clear();
    pointLights.clear();
    spotLights.clear();
}

// === MAIN PROGRAM ===
int main(int argc, char** argv) {
    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitWindowSize(renderSettings.screenWidth, renderSettings.screenHeight);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Advanced Ray Tracer - 2005091");

    // Setup scene
    initializeGraphicsSystem();
    parseSceneFile();

    // Configure OpenGL
    glEnable(GL_DEPTH_TEST);
    
    // Register callbacks
    glutDisplayFunc(renderScene);
    glutIdleFunc(refreshDisplay);
    glutSpecialFunc(handleArrowKeys);
    glutKeyboardFunc(handleKeyboard);

    // Start main loop
    glutMainLoop();
    
    // Cleanup (though this won't be reached due to glutMainLoop)
    cleanupResources();
    
    return 0;
}