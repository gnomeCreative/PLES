
#include "DEM.h"

#define USE_MATH_DEFINES

using namespace std;

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

void DEM::discreteElementGet(GetPot& lbmCfgFile, GetPot& command_line){


    PARSE_CLASS_MEMBER(lbmCfgFile, demSize[0], "demSizeX",1.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, demSize[1], "demSizeY",1.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, demSize[2], "demSizeZ",1.0);

    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[0], "boundary1",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[1], "boundary2",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[2], "boundary3",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[3], "boundary4",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[4], "boundary5",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[5], "boundary6",1);

    double readForceX(0.0),readForceY(0.0), readForceZ(0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, readForceX, "demFX",1.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, readForceY, "demFY",1.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, readForceZ, "demFZ",1.0);
    demF=tVect(readForceX,readForceY,readForceZ);

    // getting time step
    PARSE_CLASS_MEMBER(lbmCfgFile, deltat, "deltat",1.0);
    ASSERT(deltat>=0.0);

    // getting material properties
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.density, "density",1.0);
    ASSERT(sphereMat.density>0.0);

    string contactModelString;
    PARSE_CLASS_MEMBER(lbmCfgFile, contactModelString, "contactModel","none");
    if (contactModelString=="HERTZIAN") sphereMat.contactModel=HERTZIAN;
    else if (contactModelString=="LINEAR") sphereMat.contactModel=LINEAR;


    // Hertzian contact model /////////////
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.youngMod, "youngMod",1.0);
    ASSERT(sphereMat.youngMod>0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.poisson, "poisson",1.0);
    ASSERT(sphereMat.poisson>=0.0);
    // normal stiffness (constant part, calculated here to avoid repetition)
    sphereMat.knConst=2.0/3.0*sphereMat.youngMod/(1-sphereMat.poisson*sphereMat.poisson);
    sphereMat.ksConst=2.0*sphereMat.youngMod/(2.0-sphereMat.poisson)*(1.0+sphereMat.poisson);

    // linear contact model /////////////
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.linearStiff, "linearStiff",1.0);
    ASSERT(sphereMat.linearStiff>=0.0);

    // normal damping ///////////////////////
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.restitution, "restitution",1.0);
    ASSERT(sphereMat.restitution>0.0);
    ASSERT(sphereMat.restitution<=1.0);
    // calculating coefficient for normal damping
    switch (sphereMat.contactModel) {
        case HERTZIAN: {
            // see "On the Determination of the Damping Coefficient of Non-linear Spring-dashpot System to Model Hertz Contact for Simulation by Discrete Element Method"
            // Hu, Hu, Jian, Liu, Wan, Journal of Computers, 6 (2011) OR BETTER Antypov & Elliott
            sphereMat.dampCoeff=-1.0*sqrt(5)*log(sphereMat.restitution)/sqrt((log(sphereMat.restitution)*log(sphereMat.restitution)+M_PI*M_PI));
        }
        case LINEAR: {
            sphereMat.dampCoeff=-1.0*log(sphereMat.restitution)/sqrt((log(sphereMat.restitution)*log(sphereMat.restitution)+M_PI*M_PI));
        }
    }

    // tangential model //////////////////////
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.viscTang, "viscTang",1.0);
    ASSERT(sphereMat.viscTang>=0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.frictionCoefPart, "frictionCoefPart",1.0);
    ASSERT(sphereMat.frictionCoefPart>=0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.frictionCoefWall, "frictionCoefWall",1.0);
    ASSERT(sphereMat.frictionCoefWall>=0.0);



    // particle initial state //////////////////////
    string particleFile;
    PARSE_CLASS_MEMBER(lbmCfgFile, particleFile, "particleFile","particles.dat");
    double translateX(0.0), translateY(0.0), translateZ(0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, translateX, "translateX",0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, translateY, "translateY",0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, translateZ, "translateZ",0.0);
    tVect translate(translateX,translateY,translateZ);
    double scaleX(1.0), scaleY(1.0), scaleZ(1.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, scaleX, "scaleX",1.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, scaleY, "scaleY",1.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, scaleZ, "scaleZ",1.0);
    tVect scale(scaleX,scaleY,scaleZ);
    double scaleMin=std::min(std::min(scaleX,scaleY),scaleZ); // FIX THIS!
    cout<<"scale factor="<<scaleMin<<endl;

    ifstream particleFileID;
    particleFileID.open(particleFile.c_str(),ios::in);
    ASSERT(particleFileID.is_open());

    unsigned int totElmt;
    particleFileID>>totElmt;

    for (int n=0; n<totElmt; ++n) {
        elmt dummyElmt;
        double x,y,z;

        // import variables
        particleFileID>>dummyElmt.index;
        particleFileID>>dummyElmt.size;
        particleFileID>>dummyElmt.radius;
        dummyElmt.radius=dummyElmt.radius*scaleMin;
        particleFileID>>x;
        particleFileID>>y;
        particleFileID>>z;
        dummyElmt.x0=tVect(x*scaleX,y*scaleY,z*scaleY)+translate;
        particleFileID>>x;
        particleFileID>>y;
        particleFileID>>z;
        dummyElmt.x1=tVect(x,y,z);

        elmts.push_back(dummyElmt);
    }

    // objects initial state //////////////////////
    string objectFile;
    PARSE_CLASS_MEMBER(lbmCfgFile, objectFile, "objectFile","objects.dat");
    ifstream objectFileID;
    objectFileID.open(objectFile.c_str(),ios::in);
    ASSERT(objectFileID.is_open());

    unsigned int totObjects;
    objectFileID>>totObjects;

    for (int n=0; n<totObjects; ++n) {
        object dummyObject;
        double x,y,z;
        double size; //dummy, just to use same style

        // import variables
        objectFileID>>dummyObject.index;
        objectFileID>>size; // must be one
        objectFileID>>dummyObject.radius;
        objectFileID>>x;
        objectFileID>>y;
        objectFileID>>z;
        dummyObject.x0=tVect(x,y,z);
        objectFileID>>x;
        objectFileID>>y;
        objectFileID>>z;
        dummyObject.x1=tVect(x,y,z);

        objects.push_back(dummyObject);
    }

    // numerical viscosity for stability
    PARSE_CLASS_MEMBER(lbmCfgFile, numVisc, "numVisc",0.0);

}

void DEM::discreteElementInit(){

    // initializing DEM parameters from lattice parameters
    cout<<endl;
    cout<<"DEM INITIALIZATION"<<endl;

    // DEM time step
    // if multistep is 0, it should be calculated by the program here
    if (elmts.size()) {
        if (deltat==0) {
            // find critical deltaT
            const double crit=0.005*criticalTimeStep();
            deltat=crit;
        }
        // multistep can also be imposed by the user
    }

    // initializing particles
    const double partDensity=sphereMat.density;

    // initializing composite particle properties
    compositeProperties();
    // clear particle list
    particles.clear();

    unsigned int globalIndex=0;
    for (int n=0; n<elmts.size(); ++n){
        // initialize element
        elmts[n].initialize(partDensity, prototypes, demF);
        // generate particles
        elmts[n].generateParticles(globalIndex, particles, prototypes);
    }
    // the number of standard particles (=not ghosts) is now fixed, and WILL NOT BE CHANGED
    stdParticles=particles.size();

    // initializing wall for DEM
    initializeWalls();
    // initializing cylinders for DEM
    initializeCylinders();
    // initializing periodic boundary conditions for DEM
    initializePbcs();

    // initialize neighbor list parameters (also sets ghosts)
    if (elmts.size()) {
        initNeighborParameters();
        evalNeighborTable();
    }

    double totMass=0.0;
    for (int n=0; n<elmts.size(); ++n){
        // calculate mass
        totMass+=elmts[n].m;
    }

    cout<<"DEM parameters\n";
    cout<<"Tot elements: "<<elmts.size()<<";\t";
    cout<<"Tot standard particles: "<<stdParticles<<endl;
    cout<<"Deltat ="<<deltat<<endl;
    cout<<"Total mass ="<<totMass<<endl;
    switch (sphereMat.contactModel) {
        case LINEAR: {
            cout<<"Contact model: linear dashpot"<<endl;
            cout<<"Normal stiffness = "<<sphereMat.linearStiff<<endl;
            break;
        }
        case HERTZIAN: {
            cout<<"Contact model: damped Hertzian contact"<<endl;
            break;
        }
    }
    cout<<"Damping ratio = "<<sphereMat.dampCoeff<<", equivalent to a coefficient of restitution of c="<<sphereMat.restitution<<endl;
    cout<<"Tangential viscosity = "<<sphereMat.viscTang<<endl;
    cout<<"Particle-particle friction = "<<sphereMat.frictionCoefPart<<", wall-particle friction = "<<sphereMat.frictionCoefWall<<endl;
}

void DEM::discreteElementStep(){

    // set trigger for new neighbor list
    const double neighListTrigger=0.25*nebrRange;//0.25*nebrRange;   // 0.25 -> update particle when particle is inside 0.25 other particle domain

    // neighbor management
    evalMaxDisp();
    if (maxDisp>neighListTrigger) { // maxDisp>0.25*(nebrRange-2.0*cutOff)
        maxDisp=0.0;
        evalNeighborTable();
    }

    // predictor step
    predictor();

    // particles generation
    updateParticlesPredicted();

    // force evaluation
    evalForces();

    // corrector step
    corrector();

    // particles re-generation
    updateParticlesCorrected();

}

void DEM::reset(){
    // resets velocity so that the next configuration can be used as initial configuration
    for (int n=0; n<elmts.size(); ++n) {
        // initialize element
        elmts[n].resetVelocity();
    }
    updateParticlesCorrected();
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// initialization functions

void DEM::compositeProperties() {

    vecList prototype1,prototype2, prototype3, prototype4;


     // prototypes for shapes
    // every vector defines the position of a particle in the object reference frame
    // unit is radius

    prototypes.resize(5);
    prototype1.resize(1);
    prototype1[0].reset();
    prototypes[1]=prototype1;
    prototype2.resize(2);
    prototype2[0]=tVect(0.5,0.0,0.0);
    prototype2[1]=tVect(-0.5,0.0,0.0);
    prototypes[2]=prototype2;
    prototype3.resize(3);
    prototype3[0].reset();
    prototype3[1]=tVect(1.0,0.0,0.0);
    prototype3[2]=tVect(-1.0,0.0,0.0);
    prototypes[3]=prototype3;
    prototype4.resize(4);
    prototype4[0]=tVect(0.0,0.0,1.0);
    prototype4[1]=tVect(0.0,2.0*sqrt(2)/3.0,-1.0/3.0);
    prototype4[2]=tVect(2.0*sqrt(6)/6.0,-2.0*sqrt(2)/6.0,-1.0/3.0);
    prototype4[3]=tVect(-2.0*sqrt(6)/6.0,-2.0*sqrt(2)/6.0,-1.0/3.0);
    prototypes[4]=prototype4;

}

void DEM::initializeWalls() {
    walls.clear();
    wall dummyWall;
    unsigned int index=0;
    double boxVel=0.15/1.5;  // to delete

    // wall #0
    if ((boundary[0]==5)||(boundary[0]==6)||(boundary[0]==7)||(boundary[0]==8)) {
        dummyWall.p=tVect(0.0,0.0,0.0);
        dummyWall.n=tVect(1.0,0.0,0.0);
        dummyWall.index=index;
        dummyWall.flag=false;
        dummyWall.translating=false;
        dummyWall.trans.reset();
        if (boundary[0]==5) {
            dummyWall.moving=false;
            dummyWall.slip=true;
        }
        else if (boundary[0]==6) {
            dummyWall.moving=true;
            dummyWall.slip=true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        }
        else if (boundary[0]==7) {
            dummyWall.moving=false;
            dummyWall.slip=false;
        }
        else if (boundary[0]==8) {
            dummyWall.moving=true;
            dummyWall.slip=false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel=tVect(0.0,0.0,0.0); // HERE
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #1
    if ((boundary[1]==5)||(boundary[1]==6)||(boundary[1]==7)||(boundary[1]==8)) {
        dummyWall.p=tVect(demSize[0],0.0,0.0);
        dummyWall.n=tVect(-1.0,0.0,0.0);
        dummyWall.index=index;
        dummyWall.flag=false;
        dummyWall.translating=false;
        dummyWall.trans.reset();
        if (boundary[1]==5) {
            dummyWall.moving=false;
            dummyWall.slip=true;
        }
        else if (boundary[1]==6) {
            dummyWall.moving=true;
            dummyWall.slip=true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        }
        else if (boundary[1]==7) {
            dummyWall.moving=false;
            dummyWall.slip=false;
        }
        else if (boundary[1]==8) {
            dummyWall.moving=true;
            dummyWall.slip=false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel=tVect(0.0,0.0,0.0); // HERE
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #2
    if ((boundary[2]==5)||(boundary[2]==6)||(boundary[2]==7)||(boundary[2]==8)) {
        dummyWall.p=tVect(0.0,0.0,0.0);
        dummyWall.n=tVect(0.0,1.0,0.0);
        dummyWall.index=index;
        dummyWall.flag=false;
        dummyWall.translating=false;
        dummyWall.trans.reset();
        if (boundary[2]==5) {
            dummyWall.moving=false;
            dummyWall.slip=true;
        }
        else if (boundary[2]==6) {
            dummyWall.moving=true;
            dummyWall.slip=true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        }
        else if (boundary[2]==7) {
            dummyWall.moving=false;
            dummyWall.slip=false;
        }
        else if (boundary[2]==8) {
            dummyWall.moving=true;
            dummyWall.slip=false;
            dummyWall.rotCenter=tVect(0.0,0.4,1.260);
            dummyWall.omega=tVect(0.0,0.0,0.0); //0.942
            dummyWall.vel=tVect(-0.0,0.0,0.0);
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #3
    if ((boundary[3]==5)||(boundary[3]==6)||(boundary[3]==7)||(boundary[3]==8)) {
        dummyWall.p=tVect(0.0,demSize[1],0.0);
        dummyWall.n=tVect(0.0,-1.0,0.0);
        dummyWall.index=index;
        dummyWall.flag=false;
        dummyWall.translating=false;
        dummyWall.trans.reset();
        if (boundary[3]==5) {
            dummyWall.moving=false;
            dummyWall.slip=true;
        }
        else if (boundary[3]==6) {
            dummyWall.moving=true;
            dummyWall.slip=true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        }
        else if (boundary[3]==7) {
            dummyWall.moving=false;
            dummyWall.slip=false;
        }
        else if (boundary[3]==8) {
            dummyWall.moving=true;
            dummyWall.slip=false;
            dummyWall.rotCenter=tVect(0.0,demSize[1]+0.4,1.260);
            dummyWall.omega=tVect(0.0,0.0,0.0); //0.942
            dummyWall.vel=tVect(-0.0,0.0,0.0);
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #4
    if ((boundary[4]==5)||(boundary[4]==6)||(boundary[4]==7)||(boundary[4]==8)) {
        dummyWall.p=tVect(0.0,0.0,0.0);
        dummyWall.n=tVect(0.0,0.0,1.0);
        dummyWall.index=index;
        dummyWall.flag=false;
        dummyWall.translating=false;
        dummyWall.trans.reset();
        if (boundary[4]==5) {
            dummyWall.moving=false;
            dummyWall.slip=true;
        }
        else if (boundary[4]==6) {
            dummyWall.moving=true;
            dummyWall.slip=true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        }
        else if (boundary[4]==7) {
            dummyWall.moving=false;
            dummyWall.slip=false;
        }
        else if (boundary[4]==8) {
            dummyWall.moving=true;
            dummyWall.slip=false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel=tVect(0.0,0.0,0.0);
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #5
    if ((boundary[5]==5)||(boundary[5]==6)||(boundary[5]==7)||(boundary[5]==8)) {
        dummyWall.p=tVect(0.0,0.0,demSize[2]);
        dummyWall.n=tVect(0.0,0.0,-1.0);
        dummyWall.index=index;
        dummyWall.flag=false;
        dummyWall.translating=true;
        dummyWall.trans=tVect(0.0,0.0,0.0);
        if (boundary[5]==5) {
            dummyWall.moving=false;
            dummyWall.slip=true;
        }
        else if (boundary[5]==6) {
            dummyWall.moving=true;
            dummyWall.slip=true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        }
        else if (boundary[5]==7) {
            dummyWall.moving=false;
            dummyWall.slip=false;
        }
        else if (boundary[5]==8) {
            dummyWall.moving=true;
            dummyWall.slip=false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel=tVect(0.0,0.0,0.0);
        }
        ++index;
        walls.push_back(dummyWall);
    }

    // additional walls
    switch (problemName) {
        case DRUM: {
            // left wall
            if ((boundary[2]==5)||(boundary[2]==6)||(boundary[2]==7)||(boundary[2]==8)) {
                dummyWall.p=tVect(0.0+0.6+0.35,0.0,0.0);
                dummyWall.n=tVect(0.0,1.0,0.0);
                dummyWall.index=index;
                dummyWall.flag=false;
                dummyWall.translating=false;
                dummyWall.trans.reset();
                if (boundary[2]==5) {
                    dummyWall.moving=false;
                    dummyWall.slip=true;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                    dummyWall.vel.reset();
                }
                else if (boundary[2]==6) {
                    dummyWall.moving=true;
                    dummyWall.slip=true;
                    dummyWall.rotCenter=tVect(0.0+0.6+0.35,0.0,1.260);
                    dummyWall.omega=tVect(0.0,drumSpeed,0.0);
                    dummyWall.vel=tVect(-0.0,0.0,0.0);
                }
                else if (boundary[2]==7) {
                    dummyWall.moving=false;
                    dummyWall.slip=false;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                    dummyWall.vel.reset();
                }
                else if (boundary[2]==8) {
                    dummyWall.moving=true;
                    dummyWall.slip=false;
                    dummyWall.rotCenter=tVect(0.0+0.6+0.35,0.0,1.260);
                    dummyWall.omega=tVect(0.0,drumSpeed,0.0);
                    dummyWall.vel=tVect(-0.0,0.0,0.0);
                }
                ++index;
                walls.push_back(dummyWall);
            }
            // right wall
            if ((boundary[3]==5)||(boundary[3]==6)||(boundary[3]==7)||(boundary[3]==8)) {
                dummyWall.p=tVect(0.0+0.6+0.35,demSize[1],0.0);
                dummyWall.n=tVect(0.0,-1.0,0.0);
                dummyWall.index=index;
                dummyWall.flag=false;
                dummyWall.translating=false;
                dummyWall.trans.reset();
                if (boundary[3]==5) {
                    dummyWall.moving=false;
                    dummyWall.slip=true;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                    dummyWall.vel.reset();
                }
                else if (boundary[3]==6) {
                    dummyWall.moving=true;
                    dummyWall.slip=true;
                    dummyWall.rotCenter=tVect(0.0+0.6+0.35,demSize[1],1.260);
                    dummyWall.omega=tVect(0.0,drumSpeed,0.0);
                    dummyWall.vel=tVect(-0.0,0.0,0.0);
                }
                else if (boundary[3]==7) {
                    dummyWall.moving=false;
                    dummyWall.slip=false;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                    dummyWall.vel.reset();
                }
                else if (boundary[3]==8) {
                    dummyWall.moving=true;
                    dummyWall.slip=false;
                    dummyWall.rotCenter=tVect(0.0+0.6+0.35,demSize[1],1.260);
                    dummyWall.omega=tVect(0.0,drumSpeed,0.0);
                    dummyWall.vel=tVect(-0.0,0.0,0.0);
                }
                ++index;
                walls.push_back(dummyWall);
            }
            break;
        }
    }

    for (int n=0; n<walls.size(); ++n) {
        walls[n].wallShow();
    }
}

void DEM::initializeCylinders() {
    cylinders.clear();
    unsigned int index=0;

    switch (problemName) {
        case DRUM: {
            cylinder dummyCylinder;
            dummyCylinder.index=index;
            dummyCylinder.p1=tVect(0.0+0.6+0.35,0.0,1.26);
            dummyCylinder.p2=tVect(0.0+0.6+0.35,1.0,1.26);
            dummyCylinder.R=1.243;
            dummyCylinder.omega=tVect(0.0,drumSpeed,0.0);
            dummyCylinder.initAxes();
            dummyCylinder.moving=true;
            dummyCylinder.slip=false;
            cylinders.push_back(dummyCylinder);
            ++index;
            break;
        }
    }

    for (int n=0; n<cylinders.size(); ++n) {
        cylinders[n].cylinderShow();
    }
}

void DEM::initializePbcs() {

    pbcs.clear();
    pbc dummyPbc;
    unsigned int index=0;

    if (boundary[0]==4) {
        if (boundary[1]!=4) {
            cout<<"Periodic boundaries error!"<<endl;
            exit(0);
        }
        else {
            dummyPbc.p=tVect(0.0,0.0,0.0);
            dummyPbc.v=tVect(demSize[0],0.0,0.0);
            dummyPbc.setPlanes();
            dummyPbc.index=index;
            ++index;
            pbcs.push_back(dummyPbc);
        }
    }
    if (boundary[2]==4) {
        if (boundary[3]!=4) {
            cout<<"Periodic boundaries error!"<<endl;
            exit(0);
        }
        else {
            dummyPbc.p=tVect(0.0,0.0,0.0);
            dummyPbc.v=tVect(0.0,demSize[1],0.0);
            dummyPbc.setPlanes();
            dummyPbc.index=index;
            ++index;
            pbcs.push_back(dummyPbc);
        }
    }
    if (boundary[4]==4) {
        if (boundary[5]!=4) {
            cout<<"Periodic boundaries error!"<<endl;
            exit(0);
        }
        else {
            dummyPbc.p=tVect(0.0,0.0,0.0);
            dummyPbc.v=tVect(0.0,0.0,demSize[2]);
            dummyPbc.setPlanes();
            dummyPbc.index=index;
            ++index;
            pbcs.push_back(dummyPbc);
        }
    }
    for (int n=0; n<pbcs.size(); ++n) {
        pbcs[n].pbcShow();
    }
}

// integration functions

void DEM::predictor() {
    static const double c[5]={deltat, deltat*deltat/2.0, deltat*deltat*deltat/6.0, deltat*deltat*deltat*deltat/24.0, deltat*deltat*deltat*deltat*deltat/120.0};

    #pragma omp parallel for
    for (int n=0; n<elmts.size(); ++n) {
        elmts[n].predict(c);
    }
 }

void DEM::corrector() {
    static double gear[6] = {3.0/16.0, 251.0/360.0, 1.0, 11.0/18.0, 1.0/6.0, 1.0/60.0};
    static const double c[5]={deltat, deltat*deltat/2.0, deltat*deltat*deltat/6.0, deltat*deltat*deltat*deltat/24.0, deltat*deltat*deltat*deltat*deltat/120.0};
    static const double coeff[6]={gear[0]*c[1], gear[1]*c[1]/c[0], gear[2]*c[1]/c[1], gear[3]*c[1]/c[2], gear[4]*c[1]/c[3], gear[5]*c[1]/c[4]};
//    doubleList coeff;
//    coeff.resize(6);
//    coeff[0]=gear[0]*c[1];
//    coeff[1]=gear[1]*c[1]/c[0];
//    coeff[2]=gear[2]*c[1]/c[1];
//    coeff[3]=gear[3]*c[1]/c[2];
//    coeff[4]=gear[4]*c[1]/c[3];
//    coeff[5]=gear[5]*c[1]/c[4];

    #pragma omp parallel for
    for (int n=0; n<elmts.size(); n++) {
        elmts[n].correct(coeff);
    }
}

void DEM::evalForces() {

    for (int n=0; n<elmts.size(); ++n) {
        elmts[n].FParticle.reset();
        elmts[n].FWall.reset();
        elmts[n].MParticle.reset();
        elmts[n].MWall.reset();
    }

    for (int w=0; w<walls.size(); ++w) {
        walls[w].FParticle.reset();
    }

    for (int o=0; o<objects.size(); ++o) {
        objects[o].force.reset();
    }

    // forces due to particle overlap and lubrication
    particleParticleContacts();
    // forces due to contact with plane walls and lubrication
    wallParticleContacts();
    // forces due to contact with stationary spheres (objects)
    objectParticleContacts();
    // forces due to contact with cylinders
    cylinderParticelContacts();

    //  Newton equations solution
    for (int n=0; n<elmts.size(); ++n) {



        // numerical viscosity for stability
        // see "Viscous torque on a sphere under arbitrary rotation" by Lei,  Yang, and Wu, Applied Physics Letters 89, 181908 (2006)
        const tVect FVisc=-6.0*M_PI*numVisc*elmts[n].radius*elmts[n].x1;
        const tVect MVisc=-8.0*M_PI*numVisc*elmts[n].radius*elmts[n].radius*elmts[n].radius*elmts[n].wpGlobal;


        // translational motion
        // double massCoeff=(elmts[n].m-elmts[n].fluidMass)/elmts[n].m;
        // acceleration
        elmts[n].x2 = (FVisc+elmts[n].FParticle+elmts[n].FWall)/elmts[n].m + demF; //(sphereMat.density-1000.0)/sphereMat.density*

        // rotational motion
        // adjoint of orientation quaternion
        const tQuat q0adj=elmts[n].q0.adjoint();
        // rotational velocity (body-fixed reference frame)
        const tVect wBf=2.0*quat2vec( q0adj.multiply( elmts[n].q1 ) );
        // moment in global reference frame
        tVect moment=MVisc+elmts[n].MParticle+elmts[n].MWall;
        // moment in body-fixed reference frame
        moment=project(moment,elmts[n].q0.adjoint());
        // rotational acceleration (body-fixed reference frame) (Newton equation for principal system)
        const tVect waBf=newtonAcc(moment, elmts[n].I, wBf);
        // rotational acceleration (quaternion)
        const tQuat waQuat=quatAcc(waBf, elmts[n].q1);
        elmts[n].q2=0.5*elmts[n].q0.multiply(waQuat);
    }


}

void DEM::updateParticlesPredicted() {
    // this subroutine assumes that there is no change in the structure of the particles
    // any creation, deletion or reorganization of particles requires a different subroutine

//    #pragma omp parallel for
    for (int p=0; p<stdParticles; ++p) {
        //getting belonging element index
        const unsigned int clusterIndex=particles[p].clusterIndex;
        particles[p].updatePredicted(elmts[clusterIndex],prototypes);
    }

    if (ghosts.size()!=0) {
        // updating ghost particles
        for (int g=0; g<ghosts.size(); ++g) {
            //getting original particle index
            const unsigned int originParticleIndex=ghosts[g].ghostIndex;
            const unsigned int ghostParticleIndex=stdParticles+g;
            // reconstructing ghost particle
            particles[ghostParticleIndex].ghostUpdate(particles[originParticleIndex], ghosts[g].pbcVector);
        }
    }
}

void DEM::updateParticlesCorrected() {
    // this subroutine assumes that there is no change in the structure of the particles
    // any creation, deletion or reorganization of particles requires a different subroutine

//    #pragma omp parallel for
    for (int p=0; p<stdParticles; ++p) {
        //getting belonging element index
        const unsigned int clusterIndex=particles[p].clusterIndex;
        particles[p].updateCorrected(elmts[clusterIndex],prototypes);
    }

    if (ghosts.size()!=0) {
        // updating ghost particles
        for (int g=0; g<ghosts.size(); ++g) {
            //getting original particle index
            const unsigned int originParticleIndex=ghosts[g].ghostIndex;
            const unsigned int ghostParticleIndex=stdParticles+g;
            // reconstructing ghost particle
            particles[ghostParticleIndex].ghostUpdate(particles[originParticleIndex], ghosts[g].pbcVector);
        }
    }
}

double DEM::criticalTimeStep() const {
    // determines the critical time step based on the stiffness and mass of the elements
    // we use YADE documentation, see https://yade-dem.org/doc/formulation.html

    double minRad=elmts[0].radius;
    double minMass=elmts[0].m;
    for (int n=0; n<elmts.size(); ++n) {
        minRad=std::min(minRad,elmts[n].radius);
        minMass=std::min(minMass,elmts[n].m);
    }
    for (int n=0; n<elmts.size(); ++n) {
        minRad=std::min(minRad,elmts[n].radius);
        minMass=std::min(minMass,elmts[n].m);
    }


    // double const k=8.0/15.0*sphereMat.youngMod/(1-sphereMat.poisson*sphereMat.poisson)*sqrt(minRad);

    double deltaTCrit=0.0;
    switch (sphereMat.contactModel) {
        case HERTZIAN: {
            // maximum length scale
            const double maxDist=std::max(demSize[0],std::max(demSize[1],demSize[2]));
            // modulus of acceleration
            const double maxAccel=demF.norm();
            // estimate maximum velocity (assuming governed by acceleration field))
            const double maxVel=std::sqrt(2.0*maxAccel*maxDist);
            // see Landau & Lifshitz, or better Antypov & Elliott
            deltaTCrit=2.214*(2.0*minRad)*std::pow(sphereMat.density/sphereMat.youngMod,2.0/5.0)*std::pow(maxVel,1.0/5.0);
            // deltaTCrit=minRad*sqrt(sphereMat.density/sphereMat.youngMod);
            break;
        }
        case LINEAR: {
            deltaTCrit=M_PI/sqrt(sphereMat.knConst/minMass-sphereMat.dampCoeff);
            break;
        }
    }
    return deltaTCrit;
}

// neighbor list functions

void DEM::initNeighborParameters() {
    // initializes all parameters useful for neighbor list algorithm

    cout<<"Initialize neighbor list\n";
    // maximum radius
    double maxRad=0.0;
    for (int i=0; i<elmts.size(); ++i){
        if (maxRad<elmts[i].radius) {
            maxRad=elmts[i].radius;
        }
    }

    // if there are no particles, just to avoid nonsense numbers, we use an only cell
    // therefore the cell width is the same as the simulation box (which comes from the LB settings)
    if (elmts.size()==0) {
        cellWidth[0]=demSize[0];
        cellWidth[1]=demSize[1];
        cellWidth[2]=demSize[2];
    }
    // otherwise the cell size is determined by the size of the particles
    else {
        cellWidth[0]=std::min(maxRad*5.0,demSize[0]);
        cellWidth[1]=std::min(maxRad*5.0,demSize[1]);
        cellWidth[2]=std::min(maxRad*5.0,demSize[2]);
    }

    for (int k=0; k<3; ++k) {
        // ncells = numbers of cells for the linked cell algorithm
        nCells[k]=int( ceil(demSize[k]/cellWidth[k]) );
        // width of the cells (actual)
        cellWidth[k]=demSize[k]/double(nCells[k]);
        // increase by two to give a container for ghost cells, in case of periodicity
        nCells[k]+=2;
    }

    // may need a revision
    nebrRange=std::max(maxRad*3.0,0.5*std::min(cellWidth[0],std::min(cellWidth[1],cellWidth[2])));
    maxDisp=0.5*nebrRange;
    cout<<"Neighbor list parameters\n";
    cout<<"Number of Cells "<<nCells[0]<<" "<<nCells[1]<<" "<<nCells[2]<<"\n";
    cout<<"Cell width "<<cellWidth[0]<<" "<<cellWidth[1]<<" "<<cellWidth[2]<<"\n";
    cout<<"Range "<<nebrRange<<"\n";
}

void DEM::evalMaxDisp() {

    double maxVel2=0.0;
    for (int n=0; n<elmts.size(); ++n){
        const double thisVel2=elmts[n].x1.norm2();
        if (thisVel2>maxVel2){
            maxVel2=elmts[n].x1.norm2();
        }
    }
    const double maxVel=sqrt(maxVel2);
    maxDisp+=maxVel*deltat;
}

void DEM::evalCellTable() {
    // updates cellTable, a table which associates every particle to the belonging cell

    // delete precedent cell table
    cellTable.clear();
    // create a new cell table of size: number of particles + number of cells
    cellTable.resize(particles.size()+nCells[0]*nCells[1]*nCells[2]);

    // assigns value -1 to all elements of cellTable (the reason becomes clear when looking at the last two lines of this function)
    for (int n=0; n<cellTable.size(); ++n) {
        // -1 is the default value for unassigned node
        cellTable[n]=-1;
    }

    // cycle through n = cellTable size (number of particles + number of cells)
    for (int n=0; n<particles.size(); ++n) {

        // c is a identifier for the belonging cell of a particle. If particle A belongs to cell (2,5,7) in a cell
        // system with 4*6*9 cells, then it will be c = particleSize + 2 + (5*4 + 7)*6 = particleSize + 164 <- bullshit :D
        // basically is a system to store the cell coordinate of a particle in a line vector
        // Here the code: floor(x/cellWidth[0]) + nCells[0]*( floor(y/cellWidth[1]) + (nCells[1]*floor(z/cellWidth[2])) );
        const int c=particles[n].x0.linearizePosition(cellWidth, nCells)+particles.size();
        particles[n].tableCell=c;
        // if the coordinate exceeds the borders of the box, a message is displayed
        if (c>cellTable.size() || c<0) { // just put control over ghost, here !!!!!!!!!!!!!!!!!!!!!!
            cout<<"#neighborList, "<<c<<" initCellTable: particle "<<n<<" outside box, ignoring for force calculation."<<endl;
//            exit(0);
//            continue;
        }

        // cellTable is a structure to contain data in an efficient ways
        // every element of the array contains the pointer to another particle of the cell, until -1 is found
        // this explains why it needs to have particle.size() + cell.size() elements:
        // every particle points to another particle, and we need cell.size() elements with value -1
        cellTable[n] = cellTable[c];
        cellTable[c] = n;
    }

}

void DEM::evalNeighborTable() {
//    cout<<"NEW NEIGHBOR TABLE"<<endl;
// cout<<"Negihbor table evaluation"<<endl;
    // prototype neighbors
    const static int shifCell[14][3] = {
        {0, 0, 0},
        {1, 0, 0},
        {1, 1, 0},
        {0, 1, 0},
        {-1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {1, 1, 1},
        {0, 1, 1},
        {-1, 1, 1},
        {-1, 0, 1},
        {-1, -1, 1},
        {0, -1, 1},
        {1, -1, 1}
    };
    //updateParticles();
    // if periodicity is active: shift particles, identify ghosts and set the update signal for the LBM
    if (pbcs.size()!=0) {
        // cleaning up
        ghosts.clear();
        // resizing containers
        particles.resize(stdParticles);
        //resizing components of elements
        for (int n=0; n<elmts.size(); ++n) {
            elmts[n].components.resize(elmts[n].size);
        }
        // periodicity shift
        pbcShift();
        // regenerate particles
        updateParticlesCorrected();
        // identify ghosts
        createGhosts();
        // update particle listing to get ready for neighborList
        updateParticlesCorrected();

        // update signal for LBM
        newNeighborList=true;

    }

    // every element in its cell
    evalCellTable();
    // delete precedent neighbor table
    neighborTable.clear();
    // elements size
    unsigned int particleSize=particles.size();

    // cycle through all cells ->i0
    // variables for cell coordinate
    unsigned int i0[3];
    for (i0[0]=0; i0[0]<nCells[0]; ++i0[0]) {
        for (i0[1]=0; i0[1]<nCells[1]; ++i0[1]) {
            for (i0[2]=0; i0[2]<nCells[2]; ++i0[2]) {
                // linearized coordinate of cell0 (first cell of couples)
                const unsigned int cell0=(i0[2]*nCells[1]+i0[1])*nCells[0]+i0[0]+particleSize;
                // cycle through neighbors vectors -> [s][k]
                for (int s=0; s<14; s++) {
                        // variables for cell coordinate
                    unsigned int i1[3];
                    for (int k=0; k<3; k++) {
                        // determines the neighbor cells starting from the considered cell and the
                        // prototype for neighbors (shifCell)
                        i1[k]=i0[k]+shifCell[s][k];
                        // if neighbor cell is in upper border
                        if (i1[k]==nCells[k]) {
                            i1[k]=0;
                        }
                        // if neighbor cell is in lower border
                        if (i1[k]<0) {
                            i1[k]=nCells[k]-1;
                        }
                    }

                    // linearized coordinate of cell0 (first cell of couples)
                    const unsigned int cell1=(i1[2]*nCells[1]+i1[1])*nCells[0]+i1[0]+particleSize;

                    // this cycles through the elements of a cell, checking all neighbors
                    // the storage system is efficient but not trivial
                    int n0=cellTable[cell0];
                    while (n0>-1) {
                        int n1=cellTable[cell1];
                        while (n1>-1) {
                            // we save only the pair n1<n0 (otherwise twice the number of pairs), the pair n0<n1 (which is just the same) is excluded
                            if (cell0!=cell1 || n1<n0) {
                                // we exclude particles belonging to the same cluster
                                if (particles[n0].clusterIndex!=particles[n1].clusterIndex) {
                                    // we also exclude too far away pairs
                                    const tVect r10=particles[n1].x0-particles[n0].x0;
                                    if (r10.norm2()<nebrRange*nebrRange) {
                                        neighborTable.push_back(n0);
                                        neighborTable.push_back(n1);
                                    }
                                }
                            }
                            // to the next element (cell1)
                            n1=cellTable[n1];
                        }
                        // to the next element (cell0)
                        n0=cellTable[n0];
                    }
                }
            }
        }
    }

    if (walls.size()!=0) {
        evalNearWallTable();
    }
    if (cylinders.size()!=0) {
        evalNearCylinderTable();
    }
    if (objects.size()!=0) {
        evalNearObjectTable();
    }

}

void DEM::evalNearWallTable() {
    // evaluate the distance of all particles to the walls. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearWallTable.clear();

    for (int n=0; n<stdParticles; ++n) {
        for (wallList::iterator ip=walls.begin(); ip!=walls.end(); ++ip) {
            if (ip->dist(particles[n].x0)<nebrRange) {
                nearWallTable.push_back(n);
                break;
            }
        }
    }
}

void DEM::evalNearObjectTable() {
    // evaluate the distance of all particles to the walls. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearObjectTable.clear();

    const double nebrRange2=nebrRange*nebrRange;

    for (int n=0; n<stdParticles; ++n) {
        for (int o=0; o<objects.size(); ++o) {
            tVect x0ij = particles[n].x0 - objects[o].x0;
            if (x0ij.norm2()<nebrRange2) {
                nearObjectTable.push_back(n);
                break;
            }
        }
    }
}

void DEM::evalNearCylinderTable() {
    // evaluate the distance of all particles to the cylinders. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearCylinderTable.clear();

    for (int n=0; n<stdParticles; ++n) {
        for (cylinderList::iterator ip=cylinders.begin(); ip!=cylinders.end(); ++ip) {
            if (ip->dist(particles[n].x0) < nebrRange) {
                nearCylinderTable.push_back(n);
                break;
            }
        }
    }
}

// periodicity functions

void DEM::pbcShift() {
    //shifts elements and particle according to boundary conditions

    // this subroutine works on elements: particles need be generated/updated elsewhere
    for (int b=0; b<pbcs.size(); ++b) {
        const tVect pbcVector=pbcs[b].v;
        // cycle through nodes
        for (int n=0; n<elmts.size(); ++n) {
            const double leftDist=pbcs[b].pl1.dist(elmts[n].x0);
            const double rightDist=pbcs[b].pl2.dist(elmts[n].x0);
            // first plane of couple (left)
            if (leftDist<0.0) {
                // cout<<"Shift element number"<<elmts[n].index<<"\n";
                elmts[n].translate(pbcVector);
                // shift particles (not the ghosts)
//                for (int j=0; j<elmts[n].size; ++j) {
//                    particles[elmts[n].components[j]].x0+=pbcVector;
//                }
            }
            // second plane of couple (right)
            if (rightDist<0.0) {
//                cout<<"Shift element number"<<elmts[n].index<<"\n";
                elmts[n].translate(-1.0*pbcVector);
                // shift particles (not the ghosts)
//                for (int j=0; j<elmts[n].size; ++j) {
//                    particles[elmts[n].components[j]].x0-=pbcVector;
//                }
            }
        }
    }
}

void DEM::createGhosts() {

    ghosts.clear();

    // this subroutine implies that elements (and their non-ghost particles) have already been shifted
    for (int b=0; b<pbcs.size(); ++b) {
        const tVect pbcVector=pbcs[b].v;
        // cycle through elements
        for (int n=0; n<elmts.size(); ++n) {
            // cycle through standard particles
            for (int j=0; j<elmts[n].size; ++j) {
                const unsigned int p=elmts[n].components[j];
                // distances from the periodic walls
                const double leftDist=pbcs[b].pl1.dist(particles[p].x0);
                const double rightDist=pbcs[b].pl2.dist(particles[p].x0);
    //            ASSERT(leftDist>0);
    //            ASSERT(rightDist>0)
                // first plane of couple (left)
                if (leftDist<nebrRange) {
                    ghost dummyGhost;
                    dummyGhost.ghostIndex=particles[p].particleIndex;
                    dummyGhost.pbcVector=pbcVector;
                    ghosts.push_back(dummyGhost);
                }
                // second plane of couple (right), we copy only one time
                else if (rightDist<nebrRange) {
                    ghost dummyGhost;
                    dummyGhost.ghostIndex=particles[p].particleIndex;
                    dummyGhost.pbcVector=-1.0*pbcVector;
                    ghosts.push_back(dummyGhost);
                }
            }
        }
//        for (int p=0; p<stdParticles; ++p) {
//            // distances from the periodic walls
//            const double leftDist=pbcs[b].pl1.dist(particles[p].x0);
//            const double rightDist=pbcs[b].pl2.dist(particles[p].x0);
////            ASSERT(leftDist>0);
////            ASSERT(rightDist>0)
//            // first plane of couple (left)
//            if (leftDist<nebrRange) {
//                ghost dummyGhost;
//                dummyGhost.ghostIndex=particles[p].particleIndex;
//                dummyGhost.pbcVector=pbcVector;
//                ghosts.push_back(dummyGhost);
//            }
//            // second plane of couple (right), we copy only one time
//            else if (rightDist<nebrRange) {
//                ghost dummyGhost;
//                dummyGhost.ghostIndex=particles[p].particleIndex;
//                dummyGhost.pbcVector=-1.0*pbcVector;
//                ghosts.push_back(dummyGhost);
//            }
//        }
    }
    // now we need to check for particles at corners
    // these are particles have already been shifted by two periodic walls IMPORTANT: (the case of three periodic walls is missing)
    ghostList cornerGhosts;
    cornerGhosts.clear();
    for (int g1=0; g1<ghosts.size(); ++g1) {
        for (int g2=g1+1; g2<ghosts.size(); ++g2) {
            if (ghosts[g1].ghostIndex==ghosts[g2].ghostIndex) {
                // we need to create a corner ghost
                ghost dummyGhost;
                dummyGhost.ghostIndex=ghosts[g1].ghostIndex;
                dummyGhost.pbcVector=ghosts[g1].pbcVector+ghosts[g2].pbcVector;
                cornerGhosts.push_back(dummyGhost);
            }
        }
    }
//    cout<<cornerGhosts.size()<<endl;
    ghosts.insert(ghosts.end(), cornerGhosts.begin(), cornerGhosts.end());

    // resizing containers
    particles.resize(stdParticles);
    //resizing components of elements
    for (int t=0; t<elmts.size(); ++t) {
        elmts[t].components.resize(elmts[t].size);
    }

    // creating particles from ghosts
    // index for ghost particles
    unsigned int ind=stdParticles;
    // adding ghost particles to particle list
    for (int g=0; g<ghosts.size(); ++g) {
        //getting original particle index
        const unsigned int originParticleIndex=ghosts[g].ghostIndex;
        // reconstructing ghost particle
        particle dummyParticle=particles[originParticleIndex];
//        dummyParticle.x0+=ghosts[g].pbcVector;
        dummyParticle.particleIndex=ind;
        // adding ghost particle to components of the mother element
        elmts[dummyParticle.clusterIndex].components.push_back(dummyParticle.particleIndex);
        // updating particle list
        particles.push_back(dummyParticle);
        // updating particle index
        ind++;
    }
}

// force computation functions

void DEM::particleParticleContacts() {

    for (unsIntList::iterator iti = neighborTable.begin(); iti!=neighborTable.end(); iti = iti + 2){
        // couple of contact candidates
        unsIntList::iterator itj = iti + 1;
         // pointers to particles
        const particle *parti=&particles[*iti];
        const particle *partj=&particles[*itj];
        // checking for overlap
        const double ri=parti->r;
        const double rj=partj->r;
        const double sigij=ri+rj;
        const double sigij2=sigij*sigij;
        // distance between centers
        const tVect vectorDistance=partj->x0-parti->x0;
        const double distance2=vectorDistance.norm2();
        // check for contact
        if (distance2<sigij2) {
            particleParticleCollision(parti, partj, vectorDistance);
        }
    }
}

void DEM::wallParticleContacts() {

    // to keep conventions, the index i refers to the wall, and j the particle

    // cycling through walls
    for (wallList::iterator ip = walls.begin(); ip!=walls.end(); ip++) {
        wall *wallI=&*ip;
        // cycling through particle in wall neighbor list
        for (unsIntList::iterator it = nearWallTable.begin(); it!=nearWallTable.end(); it++) {
            // particle
            const particle *partJ=&particles[*it];
            // radius
            const double rj=partJ->r;
            // distance from wall (norm)
            const double distance=wallI->dist(partJ->x0);
            bool wallIsOk=true;
            if (wallIsOk) {
                // distance before contact
                const double overlap=rj-distance;
                if (overlap>0.0){
                    wallParticleCollision(wallI, partJ, overlap);
                }
            }
        }
    }
}

void DEM::cylinderParticelContacts() {

    // to keep conventions, the index i refers to the cylinder, and j the particle

    // cycling through the cylinders
    for (cylinderList::iterator ip = cylinders.begin(); ip!=cylinders.end(); ip++) {
        cylinder *cylinderI=&*ip;
        // cycling through the cylinder neighbor particles
        for (unsIntList::iterator it = nearCylinderTable.begin(); it!=nearCylinderTable.end(); it++) {
            // particle
            const particle *partJ=&particles[*it];
            // radius
            const double rj=partJ->r;
            // distance from wall (norm)
            const double distance=cylinderI->dist(partJ->x0);
            // distance before contact
            const double overlap =rj-distance;
            // check for contact
            if (overlap>0.0) {
                cylinderParticleCollision(cylinderI, partJ, overlap);
            }
        }
    }
}

void DEM::objectParticleContacts() {

    // to keep conventions, the index i refers to the object, and j the particle

    for (objectList::iterator io=objects.begin(); io!=objects.end(); io++) {
        object *objectI=&*io;
        for (unsIntList::iterator it=nearObjectTable.begin(); it!=nearObjectTable.end(); it++) {
            // particle
            const particle *partJ=&particles[*it];
            // radius
            const double rj=partJ->r;
            // distance from object (vector)
            const tVect vectorDistance=partJ->x0-objectI->x0;
            // distance from object (norm)
            const double distance=vectorDistance.norm();
            // distance before contact
            double overlap=rj+objectI->radius-distance;
            if (overlap>0.0) {
//                cout<<demTime<<"contact!\n";
//                cin.get();
                objectParticleCollision(objectI, partJ, vectorDistance);
            }
        }
    }
}

inline void DEM::particleParticleCollision(const particle *partI, const particle *partJ, const tVect& vectorDistance) {

    // pointers to elements
    elmt *elmtI=&elmts[partI->clusterIndex];
    elmt *elmtJ=&elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry /////////////////////////////
    const double radI=partI->r;
    const double radJ=partJ->r;
    // distance norm
    const double distance=vectorDistance.norm();
    // overlap
    const double overlap=radI+radJ-distance;
    // relative velocity
    const tVect relVel=partJ->x1-partI->x1;
    // first local unit vector (normal)
    const tVect en=vectorDistance/distance;
    // relative normal velocity (modulus)
    const double normRelVel=relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel=en*normRelVel;
    // effective mass
    const double effMass=elmtI->m*elmtJ->m/(elmtI->m+elmtJ->m);
    // effective radius
    const double effRad=radI*radJ/(radI+radJ);

    // force computation /////////////////////////////////
    const double normNormalForce=normalContact(overlap, normRelVel, effRad, effMass);

    //    // Overlap elastic potential energy
//    energy.elastic+=0.5*fn*xi*xi;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce=en*normNormalForce;

    // cluster geometry
    // vectorized radii
    const tVect vecRadI=radI*en;
    const tVect vecRadj=-radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistI=vecRadI+(partI->radiusVec);
    const tVect centerDistJ=vecRadj+(partJ->radiusVec);

    if (partI->particleIndex<stdParticles) {
        elmtI->FParticle=elmtI->FParticle-normalForce;
        //  moment generated in non-spherical particles
        elmtI->MParticle=elmtI->MParticle-centerDistI.cross(normalForce);
    }
    if (partJ->particleIndex<stdParticles) {
        elmtJ->FParticle=elmtJ->FParticle+normalForce;
        //  moment generated in non-spherical particles
        elmtJ->MParticle=elmtJ->MParticle+centerDistJ.cross(normalForce);
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wI=elmtI->wpGlobal; //2.0*quat2vec( elmtI->qp1.multiply( elmtI->qp0.adjoint() ) );
    const tVect wJ=elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact=relVel-wI.cross(vecRadI)+wJ.cross(vecRadj);
    // tangential component of relative velocity
    const tVect tangRelVelContact=relVelContact-normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact=tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact!=0.0) {
        // tangential force
        double normTangForce=tangentialContact(normTangRelVelContact, normNormalForce, effRad, effMass, sphereMat.frictionCoefPart);
        // second local unit vector (tangential)
        const tVect et=tangRelVelContact/normTangRelVelContact;
        // vectorial tangential force
        const tVect tangForce=normTangForce*et;
        // torque updating
        if (partI->particleIndex<stdParticles) {
            elmtI->MParticle=elmtI->MParticle+centerDistI.cross(tangForce);
            elmtI->FParticle=elmtI->FParticle+tangForce;
        }
        if (partJ->particleIndex<stdParticles) {
            elmtJ->MParticle=elmtJ->MParticle-centerDistJ.cross(tangForce);
            elmtJ->FParticle=elmtJ->FParticle-tangForce;
        }
    }
}

inline void DEM::wallParticleCollision(wall *wallI, const particle *partJ, const double& overlap) {

   // pointers to element
    elmt *elmtJ=&elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ=partJ->r;
    // first local unit vector (normal)
    const tVect en=wallI->n;
    // speed of the wall at contact point
    const tVect contactPointVelocity=wallI->getSpeed(partJ->x0); // fix this, contact point not defined
    // relative velocity
    const tVect relVel=partJ->x1-contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel=relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel=en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce=normalContact(2.0*overlap, normRelVel, radJ, elmtJ->m);


    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce=en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ=-radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ=vecRadJ+(partJ->x0-elmtJ->xp0);

    // force updating
    elmtJ->FWall=elmtJ->FWall+normalForce;
    wallI->FParticle=wallI->FParticle-normalForce;
    // torque updating
    elmtJ->MWall=elmtJ->MWall+centerDistJ.cross(normalForce);

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ=elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact=relVel+wJ.cross(vecRadJ);
    // tangential component of relative velocity
    const tVect tangRelVelContact=relVelContact-normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact=tangRelVelContact.norm();

    // checking if there is any tangential motion
    if (normTangRelVelContact!=0.0) {
        // tangential force
        double normTangForce=tangentialContact(normTangRelVelContact, normNormalForce, radJ, elmtJ->m, sphereMat.frictionCoefWall);
        // second local unit vector (tangential)
        const tVect et=tangRelVelContact/tangRelVelContact.norm();
        // vectorial tangential force
        const tVect tangForce=normTangForce*et;

        // torque updating
        elmtJ->MWall=elmtJ->MWall-centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FWall=elmtJ->FWall-tangForce;
        wallI->FParticle = wallI->FParticle+tangForce;
    }
}

inline void DEM::cylinderParticleCollision(cylinder *cylinderI, const particle *partJ, const double& overlap) {

   // pointers to element
    elmt *elmtJ=&elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ=partJ->r;
    // vectorial distance
    const tVect vecDistance=cylinderI->vecDist(partJ->x0);
    // contact point
    const tVect contactPoint=partJ->x0-vecDistance;
    // first local unit vector (normal)
    const tVect en=vecDistance/(radJ-overlap);
    // speed of the cylinder at contact point
    const tVect contactPointVelocity=cylinderI->getSpeed(contactPoint);
    // relative velocity
    const tVect relVel = partJ->x1-contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel=relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel=en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce=normalContact(2.0*overlap, normRelVel, radJ, elmtJ->m);

    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce=en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ=-radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ=vecRadJ+(partJ->x0-elmtJ->xp0);

    // force updating
    elmtJ->FWall=elmtJ->FWall+normalForce;
    // wallI->FParticle=wallI->FParticle-fnv;
    // torque updating
    elmtJ->MWall=elmtJ->MWall+centerDistJ.cross(normalForce);

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ=elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact=relVel+wJ.cross(vecRadJ); // couldn't we just use elmtJ.w?
    // tangential component of relative velocity
    const tVect tangRelVelContact=relVelContact-normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact=tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact!=0.0) {
        // tangential force
        double normTangForce=tangentialContact(normTangRelVelContact, normNormalForce, radJ, elmtJ->m, sphereMat.frictionCoefWall);
        // second local unit vector (tangential)
        const tVect et=tangRelVelContact/tangRelVelContact.norm();
        // vectorial tangential force
        const tVect tangForce=normTangForce*et;

        // torque updating
        elmtJ->MWall=elmtJ->MWall-centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FWall=elmtJ->FWall-tangForce;
        //wallI->FParticle = wallI->FParticle+ftv;
    }
}

inline void DEM::objectParticleCollision(object *objectI, const particle *partJ, const tVect& vectorDistance) {

   // pointers to element
    elmt *elmtJ=&elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ=partJ->r;
    // distance from object (norm)
    const double distance=vectorDistance.norm();
    // distance before contact
    double overlap=radJ+objectI->radius-distance;
    // first local unit vector (normal)
    const tVect en=(1.0/distance)*vectorDistance;
    // speed of the wall at contact point
    const tVect contactPointVelocity=objectI->x1; // fix this, contact point not defined
    // relative velocity
    const tVect relVel=partJ->x1-contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel=relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel=en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce=normalContact(overlap, normRelVel, radJ, elmtJ->m);

    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce=en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ=-radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ=vecRadJ+(partJ->x0-elmtJ->xp0);

    // force updating
    elmtJ->FWall=elmtJ->FWall+normalForce;
    objectI->force=objectI->force-normalForce;
    // torque updating
    elmtJ->MWall=elmtJ->MWall+centerDistJ.cross(normalForce);

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ=elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact=relVel+wJ.cross(vecRadJ);
    // tangential component of relative velocity
    const tVect tangRelVelContact=relVelContact-normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact=tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact!=0.0) {
        // tangential force
        double normTangForce=tangentialContact(normTangRelVelContact, normNormalForce, radJ, elmtJ->m, sphereMat.frictionCoefPart);
        // second local unit vector (tangential)
        const tVect et=tangRelVelContact/tangRelVelContact.norm();
        // vectorial tangential force
        const tVect tangForce=normTangForce*et;

        // torque updating
        elmtJ->MWall=elmtJ->MWall-centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FWall=elmtJ->FWall-tangForce;
        objectI->force = objectI->force+tangForce;
    }

}

double DEM::normalContact(const double& overlap, const double& vrelnnorm, const double& effRad, const double& effMass) const {

    // total normal force
    double fn=0.0;

    switch (sphereMat.contactModel) {
        case HERTZIAN: {
            // square root of effective radius
            const double sqrtEffRad=sqrt(effRad);
            // square root of overlap
            const double sqrtOverlap=sqrt(overlap);
            // normal stiffness = 2/3*(Reff^1/2)*Y/(1-nu²)
            const double kn=sphereMat.knConst*sqrtEffRad*sqrtOverlap;
            // damping
            const double gamman=2.0*sphereMat.dampCoeff*sqrt(kn*effMass);
            const double dumpfn=-gamman*vrelnnorm;
            // elastic normal force (Hertzian contact)
            const double elasticfn=kn*overlap;
            // total normal force (limit to repulsive force)
            fn=std::max(elasticfn+dumpfn,0.0);
            break;
        }
        case LINEAR: {
            // elastic normal force (linear contact)
            const double elasticfn=sphereMat.linearStiff*overlap;
            // damping
            const double gamman=2.0*sphereMat.dampCoeff*sqrt(sphereMat.linearStiff*effMass);
            const double dumpfn=-gamman*vrelnnorm;
            // total normal force (limit to repulsive force)
            fn=std::max(elasticfn+dumpfn,0.0);
            break;
        }
    }

    return fn;

}

double DEM::tangentialContact(const double& vreltNorm, const double& fn, const double& effRad, const double& effMass, const double& friction) const {

    // tangential force
    double fs=0.0;

    switch (sphereMat.contactModel) {
        case HERTZIAN: {
            static const double power=1.0/3.0; // plus or minus???????
            // square root of effective radius
            const double sqrtEffRad=sqrt(effRad);
            // tangent stiffness -> not physically sound, static friction is actually missing
            const double ks=sphereMat.ksConst*sqrtEffRad*pow(fn,power);
            // maximum force due to dynamic friction
            const double fsMax=friction *fn;
            // viscous force
            const double viscousForce=2.0*sphereMat.viscTang*sqrt(effMass*ks)*vreltNorm;
            // capping with dynamic friction
            fs=std::min(viscousForce, fsMax);
            break;
        }
        case LINEAR: {
            // tangent stiffness -> not physically sound, static friction is actually missing
            const double ks=sphereMat.linearStiff;
            // maximum force due to dynamic friction
            const double fsMax=friction *fn;
            // viscous force
            const double viscousForce=2.0*sphereMat.viscTang*sqrt(effMass*ks)*vreltNorm;
            // capping with dynamic friction
            fs=std::min(viscousForce, fsMax);
            break;
        }
    }

    return fs;

}
