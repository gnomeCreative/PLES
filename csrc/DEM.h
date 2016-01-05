/*
 * File:   DEM.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:07 PM
 */

#ifndef DEM_H
#define	DEM_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
//
#include "vector.h"
#include "elmt.h"
#include "utils.h"
#include "getpot.h"
#include "macros.h"

using namespace std;
extern ProblemName problemName;

class DEM{
private:
    // domain size: DEM grid is orthogonal
    doubleList demSize;
    // boundaries
    unsIntList boundary;
    // neighbor list (here COMMENTS are needed)
    // neighbor parameters
    double maxDisp;
    double nebrRange;
    double nebrShell;
    double cellWidth[3];
    unsigned int nCells[3];
    // neighbor tables
    intList cellTable;
    unsIntList neighborTable;
    unsIntList nearWallTable;
    unsIntList nearObjectTable;
    unsIntList nearCylinderTable;
    // hydrodynamic forces stabilization: number of time step for average
    // numerical viscosity to avoid energy problems
    double numVisc;

    // prototype shape collection
    std::vector <vecList> prototypes;
public:
    // time step
    double deltat;
    // force field (if not constant need to be calculated inside cycle) // should not be here!!!!!
    tVect demF;
    // material for elements // should not be here!!!
    material sphereMat;
    // time
    double demTime;
    // solid walls
    wallList walls;
    // solid walls
    cylinderList cylinders;
    // periodic boundary conditions
    pbcList pbcs;
    // list with elements
    elmtList elmts;
    // list with particles
    particleList particles;
    // list with objects
    objectList objects;
    // particles contains both standard particles and ghost particles. The total number of standard particles needs to be saved
    unsigned int stdParticles;
    // list with ghost elements
    ghostList ghosts;
    // true if indices for the LB need to be regenerated
    bool newNeighborList;
    // stuff for drum
    double drumSpeed;
public:
    DEM(){
        newNeighborList=false;
        demSize.resize(3);
        demSize[0]=demSize[1]=demSize[2]=1.0;
        boundary.resize(6);
        boundary[0]=boundary[1]=boundary[2]=boundary[3]=boundary[4]=boundary[5]=1;
        demTime=0.0;
        deltat=1.0;
        demF.reset();
        walls.clear();
        pbcs.clear();
        elmts.clear();
        particles.clear();
        objects.clear();
        stdParticles=0;
        ghosts.clear();
//        ghosts.clear();
        cellTable.clear();
        neighborTable.clear();
        nearWallTable.clear();
        nearObjectTable.clear();
        nearCylinderTable.clear();
        // neighbor list variables
        maxDisp=0.0;
        nebrRange=0.0;
        cellWidth[0]=cellWidth[0]=cellWidth[0]=0.0;
        nCells[0]=nCells[0]=nCells[0]=0;
        prototypes.clear();
        // problem-related variables
        drumSpeed=0.0;
}
    void discreteElementStep();
    void discreteElementGet(GetPot& lbmCfgFile, GetPot& command_line);
    void discreteElementInit();
    void evolveBoundaries();
    void reset();
private:
    // initialization functions
    void compositeProperties();
    void initializeWalls();
    void initializeCylinders();
    void initializePbcs();
    // integration functions
    void predictor();
    void corrector();
    void evalForces();
    void updateParticlesPredicted();
    void updateParticlesCorrected();
    double criticalTimeStep() const;
    // neighbor list functions
    void evalNeighborTable(); // should be private: correct!
    void initNeighborParameters();
    void evalMaxDisp();
    void evalCellTable();
    void evalNearWallTable();
    void evalNearObjectTable();
    void evalNearCylinderTable();
    // periodicity functions
    void createGhosts();
    void pbcShift();
    // force computation functions
    void particleParticleContacts();
    void wallParticleContacts();
    void cylinderParticelContacts();
    void objectParticleContacts();
    inline void particleParticleCollision(const particle *partI, const particle *partJ, const tVect& vectorDistance);
    inline void wallParticleCollision(wall *walli, const particle *partJ, const double& overlap);
    inline void cylinderParticleCollision(cylinder *cylinderI, const particle *partJ, const double& overlap);
    inline void objectParticleCollision(object *iObject, const particle *partJ, const tVect& vectorDistance);
    double normalContact(const double& overlap, const double& vrelnnorm, const double& rEff, const double& massEff) const;
    double tangentialContact(const double& vreltNorm, const double& fn, const double& effRad, const double& effMass, const double& friction) const;
};

#endif	/* DEM_H */

