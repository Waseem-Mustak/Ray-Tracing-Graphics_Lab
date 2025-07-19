#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

using namespace std;

struct Vector3D {
    double x;
    double y;
    double z;
    double w;
};

Vector3D makeVector(double x, double y, double z, double w = 1.0) {
    Vector3D v;
    v.x = x;
    v.y = y;
    v.z = z;
    v.w = w;
    return v;
}

void readVector3D(istream &in, Vector3D &vec) {
    in >> vec.x >> vec.y >> vec.z;
    vec.w = 1.0;
}

void printVector3D(ostream &out, const Vector3D &vec) {
    out << fixed << setprecision(7)
         << vec.x << " " << vec.y << " " << vec.z << "\n";
}

Vector3D add(const Vector3D &a, const Vector3D &b) {
    return makeVector(
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    );
}

Vector3D subtract(const Vector3D &a, const Vector3D &b) {
    return makeVector(
        a.x - b.x,
        a.y - b.y,
        a.z - b.z
    );
}

Vector3D multiply(const Vector3D &v, double scalar) {
    return makeVector(
        v.x * scalar,
        v.y * scalar,
        v.z * scalar
    );
}

Vector3D divide(const Vector3D &v, double scalar) {
    if (scalar == 0.0) {
        throw runtime_error("Division by zero");
    }
    return makeVector(
        v.x / scalar,
        v.y / scalar,
        v.z / scalar
    );
}

Vector3D cross(const Vector3D &a, const Vector3D &b) {
    double cx = a.y * b.z - a.z * b.y;
    double cy = a.z * b.x - a.x * b.z;
    double cz = a.x * b.y - a.y * b.x;
    return makeVector(cx, cy, cz);
}

double dot(const Vector3D &a, const Vector3D &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void normalize(Vector3D &v) {
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len != 0.0) {
        v.x /= len;
        v.y /= len;
        v.z /= len;
    }
}

void divideByW(Vector3D &v) {
    if(v.w == 0.0) {
        return;
    }
    v.x /= v.w;
    v.y /= v.w;
    v.z /= v.w;
    v.w = 1.0;
}

void rotateVector(Vector3D &v, Vector3D axis, double angle) {
    normalize(axis);
    Vector3D orig = v;
    Vector3D term1 = multiply(orig, cos(angle));
    double adotv = dot(axis, orig);
    Vector3D term2 = multiply(axis, adotv * (1.0 - cos(angle)));
    Vector3D term3 = multiply(cross(axis, orig), sin(angle));
    Vector3D sum12 = add(term1, term2);
    v = add(sum12, term3);
    // v.w = 1.0;
}

#endif