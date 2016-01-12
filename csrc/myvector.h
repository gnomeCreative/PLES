
/*
 * File:   vector.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:11 PM
 */

#ifndef VECTOR_H
#define	VECTOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

enum ProblemName {NONE, SHEARCELL, SPLASH, AVALANCHE, DRUM, BOX, NET, DIFF};

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// BASIC TYPE DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// lists of  basic variables
typedef std::vector<unsigned int> unsIntList;
typedef std::vector<int> intList;
typedef std::vector<double> doubleList;

// vector variable type
// defined with 3 double numbers
typedef class tinyVector tVect;
// quaternion variable type
// defined with 4 double numbers (should be normalized, ||q||=1)
typedef class quaternion tQuat;
// matrix variable type
// defined with 9 double numbers
typedef class tinyMatrix tMat;

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// CLASS  DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// VECTOR CLASS
class tinyVector {
private:
    double x,y,z;
public:
    tinyVector(){
        x=0.0;
        y=0.0;
        z=0.0;
    }
    tinyVector(double a, double b, double c){
        x=a;
        y=b;
        z=c;
    }
    void get(std::ifstream& inputFile);
    void show() const;
    void printLine(std::ofstream& outputFile) const;
    void print(std::ofstream& outputFile) const;
    double max() const;
    void reset();
    // overloaded operators
    tVect operator+(const tVect& vec) const;
    tVect& operator+=(const tVect& vec);
    tVect operator-(const tVect& vec) const;
    tVect& operator-=(const tVect& vec);
    tVect operator*(const double& scalar) const;
    tVect& operator*=(const double& scalar);
    tVect operator/(const double& scalar) const;
    tVect& operator/=(const double& scalar);
    friend tVect operator*(const double& scalar, const tVect& vec);
    // mathematical operations
    tVect transport() const;
    double dot(const tVect& vec) const;
    double dot2(const tVect& vec) const;
    tVect cross(const tVect& vec) const;
    tVect compProd(const tVect& vec) const;
    tMat outer(const tVect& vec) const;
    double norm() const;
    double norm2() const;
    int linearizePosition(double cellWidth[], unsigned int nCells[]) const;
    friend tVect newtonAcc(tVect moment, tVect I, tVect wBf);
    friend tQuat vec2quat(tVect vec);
    friend tVect project(tVect vec, tQuat quat);
    // geometric position functions
    bool insideSphere(tVect center, double radius) const;
    bool insideCylinder(tVect p1, tVect naxes, double R) const;
    bool insidePlane(tVect p, tVect n) const;
    double isoparameterSphere(tVect center, double radius) const;
};


// direction vectors
#define X tVect(1.0,0.0,0.0)
#define Y tVect(0.0,1.0,0.0)
#define Z tVect(0.0,0.0,1.0)
#define Xm tVect(-1.0,0.0,0.0)
#define Ym tVect(0.0,-1.0,0.0)
#define Zm tVect(0.0,0.0,-1.0)
const tVect Zero=tVect(0.0,0.0,0.0);

// QUATERNION CLASS
class quaternion{
private:
    double q0,q1,q2,q3;
public:
    quaternion(){
        q0=1.0;
        q1=q2=q3=0.0;
    }
    quaternion(double a, double b, double c, double d){
        q0=a;
        q1=b;
        q2=c;
        q3=d;
    }
    void show() const;
    // overloaded operators
    tQuat operator+(const tQuat& quat) const;
    tQuat operator-(const tQuat& quat) const;
    tQuat operator*(const double& scalar) const;
    tQuat operator/(const double& scalar) const;
    friend tQuat operator*(const double& scalar, const tQuat& quat);
    // mathematical operations
    void normalize();
    tQuat adjoint() const;
    double norm() const;
    double norm2() const;
    tQuat inverse() const;
    tQuat multiply(tQuat r) const;
    friend tVect quat2vec(tQuat quat);
    friend tVect project(tVect vec, tQuat quat);
    friend tQuat quatAcc(tVect waBf, tQuat Q1);
};

// MATRIX CLASS
class tinyMatrix {
private:
    double m00,m01,m02,m10,m11,m12,m20,m21,m22;
public:
    tinyMatrix(){
        m00=m01=m02=0.0;
        m10=m11=m12=0.0;
        m20=m21=m22=0.0;
    }
    tinyMatrix(double a00, double a01, double a02,
                         double a10, double a11, double a12,
                         double a20, double a21, double a22){
        m00=a00;
        m01=a01;
        m02=a02;
        m10=a10;
        m11=a11;
        m12=a12;
        m20=a20;
        m21=a21;
        m22=a22;
    }
    tinyMatrix(tVect v1, tVect v2){
        tMat matrix=v1.outer(v2);
        m00=matrix.m00;
        m01=matrix.m01;
        m02=matrix.m02;
        m10=matrix.m10;
        m11=matrix.m11;
        m12=matrix.m12;
        m20=matrix.m20;
        m21=matrix.m21;
        m22=matrix.m22;
    }
    void show() const;
    // overloaded operators
    tMat operator+(const tMat& mat) const;
    tMat& operator+=(const tMat& mat);
    tMat operator-(const tMat& mat) const;
    tMat& operator-=(const tMat& mat);
    tMat operator*(const double& scalar) const;
    tMat& operator*=(const double& scalar);
    tMat operator/(const double& scalar) const;
    tMat& operator/=(const double& scalar);
    friend tMat operator*(const double& scalar, const tMat& mat);
    // mathematical operations
    double magnitude() const;
};

// INTER-CLASS FUNCTIONS
// overloaded operators
tVect operator*(const double& scalar, const tVect& vec);
tQuat operator*(const double& scalar, const tQuat& quat);
tMat operator*(const double& scalar, const tMat& mat);
// mathematical operations
tVect quat2vec(tQuat quat);
tQuat vec2quat(tVect vec);
tVect project(tVect vec, tQuat quat);
tVect newtonAcc(tVect moment, tVect I, tVect wBf);
tQuat quatAcc(tVect waBf, tQuat Q1);

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// COMPOSITE TYPE DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// lists of vectors and list of lists of vectors
typedef std::vector<tVect> vecList;
typedef std::vector< std::vector<tVect> > setList;

// list of matrices
typedef std::vector<tMat> matList;

// list of particles, elements and ghosts
class particle;
typedef std::vector<particle> particleList;
class elmt;
typedef std::vector<elmt> elmtList;
class ghost;
typedef std::vector<ghost> ghostList;
class object;
typedef std::vector<object> objectList;

// list of geometric entities
class wall;
typedef std::vector<wall> wallList;
class cylinder;
typedef std::vector<cylinder> cylinderList;
class pbc;
typedef std::vector<pbc> pbcList;

#endif	/* VECTOR_H */

