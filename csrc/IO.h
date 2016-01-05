/*
 * File:   IO.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:06 PM
 */

#ifndef IO_H
#define	IO_H


#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <iomanip>
// only to create a directory
#include <sys/stat.h>
#include <sys/types.h>
//
#include "DEM.h"
#include "vector.h"
#include "elmt.h"

using namespace std;
extern ProblemName problemName;

class IO{
public:
    // config file
    string configFile;
    // time integration
    unsigned int currentTimeStep, maximumTimeSteps;
    // maximum real time for simulation
    double maximumRealTime;
    // maximum time variable
    double maxTime;
    // intervals for output
    double screenExpTime, partExpTime, objectExpTime, recycleExpTime, outputExpTime;
    // support variables for output
    unsigned int lastScreenExp, lastPartExp, lastObjectExp, lastOutputExp, lastRecycleExp;
    // formats for output file names
    string partFileFormat, objectFileFormat, recycleFileFormat;
    // support variable to store file names
    char filePathBuffer [1024];
    // boolean for exit control
    bool problem;
    // simulation time in real units
    double realTime;
    // locations
    string workDirectory, partDirectory;
    // output file names
    string exportFileName, maxSpeedFileName, plasticityFileName;
    string forceFileName, viscFileName, diffFileName;
    // generic file output streams
    ofstream exportFile, maxSpeedFile, plasticityFile;
    ofstream forceFile, viscFile, diffFile;
    // avalanche files
    string obstacleFileName;
    ofstream obstacleFile;

    // to be fixed
    unsigned int samplingStart;
    // time-averaged quantities
    double msdSum, msdSum2;
    unsigned int msdSumN;
    tVect msdVecSum, msdVecSum2;
    //
    int timeStepCount;
public:
    void initialize();
    void outputStep(const double& demTime, const DEM& dem);
    void outputFinal();
    void exportRecycleParticles(const elmtList& elmts);
private:
    // paraview
    void createParaviewFiles(const DEM& dem);
    void exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile);
    void exportParaviewObjects(const objectList& particles, const string& objectFile);
    void exportParaviewBox(const double& time);
    // print export stuff
    void exportMaxSpeedParticles(const DEM& dem);
    void exportForces(const DEM& dem);
    void exportForceObstacle(const DEM& dem);   // obstacle not used but good to look to learn
    void exportDiffusion(const DEM& dem);
    // data elaboration
    double totParticleMass(const elmtList& elmts) const;
    tVect centerOfMassParticles(const elmtList& elmts) const;
    tVect totForceObject(const objectList& objects) const;
    void  diffusionCoefficient(const elmtList& elmts, double& msd, tVect& msdComp) const;
    double collisionForceTot(const elmtList& elmts) const;
};


#endif	/* IO_H */

