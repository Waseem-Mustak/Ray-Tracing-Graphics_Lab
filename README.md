# Advanced Ray Tracing Implementation - Graphics Lab

A comprehensive ray tracing implementation in C++ using OpenGL/GLUT for real-time 3D visualization and ray tracing image generation. This project demonstrates advanced computer graphics concepts including ray-object intersection, lighting models, shadows, reflections, and texture mapping.

## Features

### 🎯 Core Ray Tracing
- **Multiple Object Types**: Spheres, triangles, and general quadric surfaces
- **Advanced Lighting**: Phong illumination model with ambient, diffuse, and specular components
- **Recursive Reflections**: Configurable recursion depth for realistic reflections
- **Shadow Calculation**: Accurate shadow rendering with multiple light sources
- **Anti-Aliasing**: High-quality image generation with customizable resolution

### 🌟 Object Support
- **Spheres**: Full spherical objects with customizable materials
- **Triangles**: Triangle mesh support for complex geometries
- **Quadric Surfaces**: General quadric equation support (ellipsoids, cylinders, cones, etc.)
- **Textured Floor**: Reflective floor with procedural or texture-based patterns

### 💡 Lighting System
- **Point Lights**: Omnidirectional light sources with color control
- **Spot Lights**: Directional lights with cone angle control
- **Multiple Light Sources**: Support for complex lighting scenarios
- **Shadow Mapping**: Realistic shadow calculation for all light types

### 🎨 Visual Features
- **Texture Mapping**: Image-based textures with toggle functionality
- **Real-time Preview**: OpenGL-based 3D scene visualization
- **Interactive Camera**: Full 6DOF camera control
- **High-Quality Output**: Bitmap image generation at configurable resolutions

## System Requirements

- **Operating System**: Windows (PowerShell support)
- **Compiler**: GCC with C++11 support
- **Graphics Libraries**: 
  - OpenGL
  - GLUT (FreeGLUT recommended)
  - GLEW32
  - GLU32

## Installation & Setup

### 1. Install Dependencies
Ensure you have the required graphics libraries installed:
- FreeGLUT
- GLEW
- OpenGL development libraries

### 2. Compilation
Use the following command to compile the project:

```bash
g++ 2005091_main.cpp -o demo.exe -lfreeglut -lglew32 -lopengl32 -lglu32
```

### 3. Running the Application
After successful compilation, run the executable:

```bash
start demo.exe
```

## Controls & Usage

### 🎮 Camera Movement
| Key | Action | Description |
|-----|--------|-------------|
| `1` | Rotate Left | Rotate camera view left |
| `2` | Rotate Right | Rotate camera view right |
| `3` | Rotate Up | Rotate camera view up |
| `4` | Rotate Down | Rotate camera view down |
| `5` | Roll Right | Roll camera clockwise |
| `6` | Roll Left | Roll camera counter-clockwise |

### ⌨️ Movement Controls
| Key | Action | Description |
|-----|--------|-------------|
| `↑` | Move Forward | Move camera forward (fast) |
| `↓` | Move Backward | Move camera backward (fast) |
| `←` | Move Left | Move camera left (fast) |
| `→` | Move Right | Move camera right (fast) |
| `Page Up` | Move Up | Move camera up globally (fast) |
| `Page Down` | Move Down | Move camera down globally (fast) |
| `W` | Move Up | Move camera up (slow) |
| `S` | Move Down | Move camera down (slow) |
| `A` | Move Left | Move camera left (slow) |
| `D` | Move Right | Move camera right (slow) |

### 🎨 Rendering Controls
| Key | Action | Description |
|-----|--------|-------------|
| `0` | Capture Image | Generate ray traced image |
| `T` | Toggle Texture | Switch between texture and normal mode |
| `ESC` | Exit | Exit the application |

## Scene Configuration

The scene is configured through the `scene.txt` file with the following format:

```
4                    # Recursion level
800                  # Image resolution (800x800)

9                    # Number of objects
sphere               # Object type
0.0 0.0 0.0         # Center coordinates
150.0               # Radius
1.0 1.0 1.0         # Color (RGB)
0.4 0.2 0.2 0.2     # Lighting coefficients [ambient, diffuse, specular, reflection]
10                  # Shininess

# ... more objects ...

4                   # Number of point lights
70.0 70.0 70.0     # Light position
1.0 0.0 0.0        # Light color

1                  # Number of spot lights
100 100 -200       # Spot light position
0 1.0 0.0         # Spot light color
0 0 1             # Spot light direction
12                # Spot light cutoff angle
```

## Project Structure

```
📁 Project Root/
├── 📄 2005091_main.cpp          # Main application file
├── 📄 2005091_classes.hpp       # Object wrapper classes
├── 📄 camera.h                  # Camera control system
├── 📄 vector3D.h               # 3D vector operations
├── 📄 ray.h                    # Ray structure and operations
├── 📄 lights.h                 # Lighting system
├── 📄 objects.h                # Object primitives
├── 📄 bitmap_image.hpp         # Image output handling
├── 📄 stb_image.h              # Image loading library
├── 📄 scene.txt                # Scene configuration file
├── 📄 texture.jpg              # Texture file
├── 📄 demo.exe                 # Compiled executable
├── 📄 RayTrace_*.bmp           # Generated ray traced images
└── 📄 README.md                # This file
```

## Technical Implementation

### Ray Tracing Algorithm
1. **Ray Generation**: Generate rays from camera through each pixel
2. **Intersection Testing**: Find closest object intersection for each ray
3. **Lighting Calculation**: Apply Phong lighting model at intersection points
4. **Reflection Handling**: Recursively trace reflected rays
5. **Shadow Testing**: Cast shadow rays to light sources
6. **Color Composition**: Combine all lighting contributions

### Supported Object Types
- **Spheres**: Using geometric intersection formula
- **Triangles**: Barycentric coordinate intersection
- **Quadric Surfaces**: General quadratic equation solving
- **Floor**: Special case with texture/pattern support

### Lighting Model
- **Ambient**: Base illumination level
- **Diffuse**: Lambert's cosine law
- **Specular**: Phong reflection model
- **Recursive Reflection**: Mirror-like reflections with depth control

## Output Images

The application generates high-quality ray traced images saved as:
- `RayTrace_1.bmp`
- `RayTrace_2.bmp`
- `RayTrace_3.bmp`
- ... (incremental numbering)

Each image is rendered at the resolution specified in `scene.txt`.

## Development Notes

### Performance Considerations
- **Spatial Optimization**: Efficient ray-object intersection testing
- **Shadow Optimization**: Early termination for shadow rays
- **Memory Management**: Proper cleanup of dynamic objects

### Extensibility
- **New Object Types**: Easily add new primitive types
- **Lighting Models**: Modify lighting calculations in object classes
- **Texture Support**: Expandable texture mapping system

## Troubleshooting

### Common Issues
1. **Compilation Errors**: Ensure all graphics libraries are properly linked
2. **Missing textures**: Application will fall back to procedural patterns
3. **Performance**: Reduce image resolution or recursion depth for faster rendering

### Debug Features
- Console output for camera position
- Timing information for ray tracing operations
- Object count and scene statistics

## Credits

**Student ID**: 2005091  
**Course**: Computer Graphics Laboratory  
**Implementation**: Advanced Ray Tracing with Real-time Preview

---

*This project demonstrates advanced computer graphics techniques including ray tracing, lighting models, and real-time 3D visualization using modern C++ and OpenGL.*
