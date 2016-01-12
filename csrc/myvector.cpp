

#include "myvector.h"

using namespace std;

/*/////////////////////////////////////////////////////////////////////////////////////////
// VECTOR
/////////////////////////////////////////////////////////////////////////////////////////*/

// basic functions

void tinyVector::show() const {
    cout<<"("<<x<<", "<<y<<", "<<z<<")";
}

void tinyVector::printLine(std::ofstream& outputFile) const {
    outputFile<<x<<" "<<y<<" "<<z<<"\n";
}

void tinyVector::print(std::ofstream& outputFile) const {
    outputFile<<x<<" "<<y<<" "<<z<<" ";
}


void tinyVector::get(std::ifstream& inputFile) {
        inputFile>>x;
        inputFile>>y;
        inputFile>>z;
    }

double tinyVector::max() const {
    return std::max(std::abs(x),std::max(std::abs(y),std::abs(z)));
}

void tinyVector::reset() {
    x=0.0;
    y=0.0;
    z=0.0;
}

// overloaded operators

tVect tinyVector::operator+(const tVect& vec) const {
    return tVect(
        x+vec.x,
        y+vec.y,
        z+vec.z);
}

tVect& tinyVector::operator+=(const tVect& vec) {
    x+=vec.x;
    y+=vec.y;
    z+=vec.z;
    return *this;
}

tVect tinyVector::operator-(const tVect& vec) const {
    return tinyVector(
        x-vec.x,
        y-vec.y,
        z-vec.z);
}

tVect& tinyVector::operator-=(const tVect& vec) {
    x-=vec.x;
    y-=vec.y;
    z-=vec.z;
    return *this;
}

tVect tinyVector::operator*(const double& scalar) const {
    return tVect(
        x*scalar,
        y*scalar,
        z*scalar);
}

tVect& tinyVector::operator*=(const double& scalar) {
    x*=scalar;
    y*=scalar;
    z*=scalar;
    return *this;
}

tVect tinyVector::operator/(const double& scalar) const {
    return tVect(
        x/scalar,
        y/scalar,
        z/scalar);
}

tVect& tinyVector::operator/=(const double& scalar) {
    x/=scalar;
    y/=scalar;
    z/=scalar;
    return *this;
}

tVect operator *(const double& scalar, const tVect& vec) {
    return tVect(
        vec.x*scalar,
        vec.y*scalar,
        vec.z*scalar);
}

// mathematical functions

tVect tinyVector::transport() const {
    return tVect(y*y+z*z,z*z+x*x,x*x+y*y);
}

double tinyVector::dot(const tVect& vec) const {
    return x*vec.x+y*vec.y+z*vec.z;
}

double tinyVector::dot2(const tVect& vec) const {
    return pow(this->dot(vec),2.0);
}

tVect tinyVector::cross(const tVect& vec) const {
    return tVect(
            y*vec.z-z*vec.y,
            z*vec.x-x*vec.z,
            x*vec.y-y*vec.x);
}

tVect tinyVector::compProd(const tVect& vec) const {
    return tVect(
            x*vec.x,
            y*vec.y,
            z*vec.z);
}

tMat tinyVector::outer(const tVect& vec) const {
    return tMat(x*vec.x,x*vec.y,x*vec.z,y*vec.x,y*vec.y,y*vec.z,z*vec.x,z*vec.y,z*vec.z);
}

double tinyVector::norm() const {
    return sqrt(x*x+y*y+z*z);
}

double tinyVector::norm2() const {
    return x*x+y*y+z*z;
}

int tinyVector::linearizePosition(double cellWidth[], unsigned int nCells[]) const  {
    return floor(x/cellWidth[0])+1+nCells[0]*( floor(y/cellWidth[1])+1+( nCells[1]*( floor(z/cellWidth[2]))+1 ) );
}

// geometric position functions

bool tinyVector::insideSphere(tVect center, double radius) const {
    //distance between center of the sphere and point
    tVect dist=*this -center;
//    return ((z-center.z)*(z-center.z)+(y-center.y)*(y-center.y)+(x-center.x)*(x-center.x)<radius*radius);
    return (dist.norm2()<radius*radius);
}

bool tinyVector::insideCylinder(tVect p1, tVect naxes, double R) const {
    // distance to point 1 of axes
    tVect p1dist;
    p1dist=*this-p1;
    // same but projected on the axes
    tVect p1distax=(p1dist.dot(naxes))*naxes;
    // distance center to cylinder
    tVect p1distcylinder=p1dist-p1distax;
    // condition for being inside
    return (p1distcylinder.norm2()<R*R);
}

bool tinyVector::insidePlane(tVect p, tVect n) const {
    tVect pdist=*this-p;
    return n.dot(pdist)<0.0;
}

double tinyVector::isoparameterSphere(tVect center, double radius) const {
    //distance between center of the sphere and point
    tVect dist=*this -center;
//    return ((z-center.z)*(z-center.z)+(y-center.y)*(y-center.y)+(x-center.x)*(x-center.x)<radius*radius);
    return (radius*radius)/dist.norm2();
}

/*/////////////////////////////////////////////////////////////////////////////////////////
// QUATERNION
/////////////////////////////////////////////////////////////////////////////////////////*/

// basic functions

void quaternion::show() const {
    std::cout<<" ("<<q0<<", "<<q1<<", "<<q2<<", "<<q3<<")";
}

// overloaded operators

tQuat quaternion::operator+(const tQuat& quat) const {
    return tQuat(
        q0+quat.q0,
        q1+quat.q1,
        q2+quat.q2,
        q3+quat.q3);
}

tQuat quaternion::operator-(const tQuat& quat) const {
    return tQuat(
        q0-quat.q0,
        q1-quat.q1,
        q2-quat.q2,
        q3-quat.q3);
}

tQuat quaternion::operator*(const double& scalar) const {
    return tQuat(
        q0*scalar,
        q1*scalar,
        q2*scalar,
        q3*scalar);
}

tQuat quaternion::operator/(const double& scalar) const {
    return tQuat(
        q0/scalar,
        q1/scalar,
        q2/scalar,
        q3/scalar);
}

tQuat operator *(const double& scalar, const tQuat& quat){
    return tQuat(
        quat.q0*scalar,
        quat.q1*scalar,
        quat.q2*scalar,
        quat.q3*scalar);
}

// mathematical functions

void quaternion::normalize(){
    double norm=this->norm();
    q0/=norm;
    q1/=norm;
    q2/=norm;
    q3/=norm;
}

double quaternion::norm() const {
    return sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
}

double quaternion::norm2() const {
    return q0*q0+q1*q1+q2*q2+q3*q3;
}

tQuat quaternion::adjoint() const {
    return tQuat(
        q0,
        -q1,
        -q2,
        -q3);
}

tQuat quaternion::inverse() const {
    // hopefully this is not needed, since for a unit quaternion inverse and adjoint are equal
    tQuat dummyQ;
    dummyQ=dummyQ.adjoint();
    dummyQ.normalize();
    return dummyQ;
}

tQuat quaternion::multiply(quaternion quat) const {
    // note that quaternion multiplication is non-commutative
    // here the product q*r is solved, which in general is different from r*q
    // refer to the internet (mathworks page on quaternion multiplication) for documentation
    return tQuat(
        quat.q0*q0-quat.q1*q1-quat.q2*q2-quat.q3*q3,
        quat.q0*q1+quat.q1*q0-quat.q2*q3+quat.q3*q2,
        quat.q0*q2+quat.q1*q3+quat.q2*q0-quat.q3*q1,
        quat.q0*q3-quat.q1*q2+quat.q2*q1+quat.q3*q0);
}

/*/////////////////////////////////////////////////////////////////////////////////////////
// MATRIX
/////////////////////////////////////////////////////////////////////////////////////////*/

// basic functions

void tinyMatrix::show() const {
    std::cout<<" (";
    std::cout<<m00<<", "<<m01<<", "<<m02<<"; ";
    std::cout<<m10<<", "<<m11<<", "<<m12<<"; ";
    std::cout<<m20<<", "<<m21<<", "<<m22<<") ";
}

// overloaded operators

tMat tinyMatrix::operator+(const tMat& mat) const {
    return tMat(
        m00+mat.m00,
        m01+mat.m01,
        m02+mat.m02,
        m10+mat.m10,
        m11+mat.m11,
        m12+mat.m12,
        m20+mat.m20,
        m21+mat.m21,
        m22+mat.m22);
}

tMat& tinyMatrix::operator+=(const tMat& mat){
    m00+=mat.m00;
    m01+=mat.m01;
    m02+=mat.m02;
    m10+=mat.m10;
    m11+=mat.m11;
    m12+=mat.m12;
    m20+=mat.m20;
    m21+=mat.m21;
    m22+=mat.m22;
    return *this;
}

tMat tinyMatrix::operator-(const tMat& mat) const {
    return tMat(
        m00-mat.m00,
        m01-mat.m01,
        m02-mat.m02,
        m10-mat.m10,
        m11-mat.m11,
        m12-mat.m12,
        m20-mat.m20,
        m21-mat.m21,
        m22-mat.m22);
}

tMat& tinyMatrix::operator-=(const tMat& mat){
    m00-=mat.m00;
    m01-=mat.m01;
    m02-=mat.m02;
    m10-=mat.m10;
    m11-=mat.m11;
    m12-=mat.m12;
    m20-=mat.m20;
    m21-=mat.m21;
    m22-=mat.m22;
    return *this;
}

tMat tinyMatrix::operator*(const double& scalar) const {
    return tMat(
        m00*scalar,
        m01*scalar,
        m02*scalar,
        m10*scalar,
        m11*scalar,
        m12*scalar,
        m20*scalar,
        m21*scalar,
        m22*scalar);
}

tMat& tinyMatrix::operator*=(const double& scalar){
    m00*=scalar;
    m01*=scalar;
    m02*=scalar;
    m10*=scalar;
    m11*=scalar;
    m12*=scalar;
    m20*=scalar;
    m21*=scalar;
    m22*=scalar;
    return *this;
}

tMat tinyMatrix::operator/(const double& scalar) const {
    return tMat(
        m00/scalar,
        m01/scalar,
        m02/scalar,
        m10/scalar,
        m11/scalar,
        m12/scalar,
        m20/scalar,
        m21/scalar,
        m22/scalar);
}

tMat& tinyMatrix::operator/=(const double& scalar){
    m00/=scalar;
    m01/=scalar;
    m02/=scalar;
    m10/=scalar;
    m11/=scalar;
    m12/=scalar;
    m20/=scalar;
    m21/=scalar;
    m22/=scalar;
    return *this;
}

tMat operator *(const double& scalar, const tMat& mat){
    return tMat(
        mat.m00*scalar,
        mat.m01*scalar,
        mat.m02*scalar,
        mat.m10*scalar,
        mat.m11*scalar,
        mat.m12*scalar,
        mat.m20*scalar,
        mat.m21*scalar,
        mat.m22*scalar);
}

// mathematical functions

double tinyMatrix::magnitude() const {
//    return 2.0*sqrt((m01*m10+m20*m02+m12*m21-(m00*m11+m11*m22+m22*m00))); // this gives nan
    return sqrt(0.5*(m00*m00+m11*m11+m22*m22+2.0*(m01*m10+m20*m02+m12*m21))); // this WORKS
//    return 2.0*sqrt(m00*m00+m11*m11+m22*m22+2.0*(m01*m10+m20*m02+m12*m21)); // this does not work
}

/*/////////////////////////////////////////////////////////////////////////////////////////
// INTER-CLASS
/////////////////////////////////////////////////////////////////////////////////////////*/

tVect quat2vec(tQuat quat) {
    if (std::abs(quat.q0)>0.001*std::abs(quat.q1)){
//    cout<<"transform error="<<quat.q0<<"\n";
    }
    return tVect(
            quat.q1,
            quat.q2,
            quat.q3);
}

tQuat vec2quat(tVect vec){
    return tQuat(
            0.0,
            vec.x,
            vec.y,
            vec.z);
}

tVect project(tVect vec, tQuat quat){
    // projection of a vector  v from a reference frame to another.
    // the reference frame is identified by a quaternion q
    // v'=qvq*
    // v=q*v'q
    tQuat qAdjoint=quat.adjoint();
    tQuat vQuat=vec2quat(vec);
    tQuat rotQuat;
    rotQuat=quat.multiply(vQuat);
    rotQuat= rotQuat.multiply(qAdjoint);
    return quat2vec(rotQuat);
}

tVect newtonAcc(tVect moment, tVect I, tVect wBf){
    double xA=(moment.x+(I.y-I.z)*wBf.y*wBf.z)/I.x;
    double yA=(moment.y+(I.z-I.x)*wBf.z*wBf.x)/I.y;
    double zA=(moment.z+(I.x-I.y)*wBf.x*wBf.y)/I.z;
    return tVect(xA,yA,zA);
}

tQuat quatAcc(tVect waBf, tQuat Q1){
    tQuat waQuat=vec2quat(waBf);
    waQuat.q0=-2.0*Q1.norm2();
    return waQuat;
}
