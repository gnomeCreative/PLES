#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <csignal>

#include "getpot.h"

#include "IO.h"
#include "DEM.h"

/*
 *      HYBIRD_shelly
 *      DEM-LBM-FEM coupling code
 *
 *      authors: Alessandro Leonardi & Roman Vetter
 *
 */


ProblemName problemName=NONE;

void goCycle(IO& io, DEM& dem) {

    // advance one step in time
    io.outputStep(dem.demTime, dem);

    ++io.currentTimeStep;
    dem.demTime+=dem.deltat;

    dem.discreteElementStep();

}

void parseCommandLine(IO& io, GetPot& command_line) {

    // create the data directory in case it doesn't exist
    io.workDirectory = "data";
    if (command_line.search("-d")) {
        io.workDirectory = command_line.next(io.workDirectory);
    }

    // if the name of the simulation is specified create subfolder
    std::string programName = "name";
    if (command_line.search("-n")) {
        // creating work directory
        programName=command_line.next(programName);
        // if the specified name is "time" then create a folder named with the initial time
        if (programName=="time") {
            string simulationTime;
            time_t t;
            time(&t);
            simulationTime=ctime(&t);
            io.workDirectory=io.workDirectory+"/"+simulationTime;
            io.workDirectory=io.workDirectory.substr(0,io.workDirectory.size()-1);
        }
        else {
            io.workDirectory=io.workDirectory+"/"+programName;
        }

//        io.workDirectory = command_line.next(io.workDirectory);
    }

    // lbm configuration file
    io.demConfigFile = "config.cfg";
    if (command_line.search("-c")) {
        io.demConfigFile = command_line.next(io.demConfigFile.c_str());
        cout<<"Using "<<io.demConfigFile<<" for DEM parameters\n";
        // copy the lbm configuration file into the data folder for later reproducibility
        std::system((std::string("mkdir -p '") + io.workDirectory + "'; cp '" + io.demConfigFile + "' '" + io.workDirectory + "'").c_str());
    }

    // make sure the config files can be read
    std::ifstream lbmCfgFile(io.demConfigFile.c_str());
    if (!lbmCfgFile) {
        cout << "ERROR: Can't open config file \"" << io.demConfigFile << "\" for reading!\n";
//        return ERROR;
    }
    lbmCfgFile.close();
}

void parseDemConfigFile(IO& io, DEM& dem, GetPot& lbmCfgFile, GetPot& command_line) {

    // PROBLEM NAME //////////////
    // necessary for hard coded sections of the code
    string problemNameString;
    PARSE_CLASS_MEMBER(lbmCfgFile, problemNameString, "problemName","none");
    if (problemNameString=="DRUM") problemName=DRUM;
    else if (problemNameString=="SHEARCELL") problemName=SHEARCELL;
    else if (problemNameString=="NONE") problemName=NONE;
    else if (problemNameString=="AVALANCHE") problemName=AVALANCHE;
    else if (problemNameString=="SPLASH") problemName=SPLASH;
    else if (problemNameString=="BOX") problemName=BOX;
    else if (problemNameString=="NET") problemName=NET;
    else if (problemNameString=="DIFF") problemName=DIFF;
    // else if (problemNameString=="SETT") problemName=SETT;


    // GETTING SIMULATION PARAMETERS  /////////
    // maximum time variable value
    PARSE_CLASS_MEMBER(lbmCfgFile, io.maxTime, "maxTime", 0.0);
    ASSERT(io.maxTime >= 0);
    // output parameters
    PARSE_CLASS_MEMBER(lbmCfgFile, io.samplingStart, "samplingStart",0);
    ASSERT(io.samplingStart >= 0);
    // stuff taken from shell code
    PARSE_CLASS_MEMBER(lbmCfgFile, io.maximumTimeSteps, "time_steps", 0);
    ASSERT(io.maximumTimeSteps >= 0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.maximumRealTime, "walltime_limit", 0);
    ASSERT(io.maximumRealTime >= 0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.screenExpTime, "screenExpTime", 0.0);
    ASSERT(io.screenExpTime  >= 0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.partExpTime, "partExpTime", 0.0);
    ASSERT(io.partExpTime   >= 0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.objectExpTime, "objectExpTime", 0.0);
    ASSERT(io.objectExpTime   >= 0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.recycleExpTime, "recycleExpTime", 0.0);
    ASSERT(io.recycleExpTime   >= 0);

    // getting DEM parameters and particle initial status
    dem.discreteElementGet(lbmCfgFile, command_line);

    // problem specific parameters
    switch (problemName) {
        case DRUM: {
            // drum rotational speed
            PARSE_CLASS_MEMBER(lbmCfgFile, dem.drumSpeed, "drumSpeed", 0.0);
            cout<<"DRUM SPEED ="<<dem.drumSpeed<<"\n";
            break;
        }
    }
}

void printUfo(GetPot& command_line, GetPot& lbmCfgFile){
    // warn about unused parameters

    std::vector<std::string> ufos=lbmCfgFile.unidentified_variables();
    if (ufos.size() > 0) {
        cout << INFO << "WARNING: Unrecognized lbm config file parameter(s):";
        for (unsigned int i = 0; i < ufos.size(); ++i) {
            cout << " " << ufos[i];
        }
        cout << endl;
    }

    ufos = command_line.unidentified_arguments();
    if (ufos.size() > 0) {
        cout << INFO << "WARNING: Unrecognized command line argument(s):";
        for (unsigned int i = 0; i < ufos.size(); ++i) {
            cout << " " << ufos[i];
        }
        cout << endl;
    }
}

int demMain(int argc, char** argv) {

    // Checking number of processes involved
    cout<<"Program starts with threads:";
    #pragma omp parallel  // create parallel memory
    {
        #pragma omp critical  // intruduce command to be carried out sequentially (e.g. for a for loop)
        cout<<" X ";
    }
    cout<<"\n";

    // DECLARATION OF VARIABLES - Input-Output ///////////////
    IO io;

    // DECLARATION OF VARIABLES - DEM ///////////////
    DEM dem;

    time_t start_time(time(NULL));

    // parsing command line and configuration file
    GetPot command_line(argc, argv);
    parseCommandLine(io, command_line);

    // parsing LBM input file
    GetPot lbmCfgFile(io.demConfigFile);
    parseDemConfigFile(io, dem, lbmCfgFile, command_line);

    io.currentTimeStep = 0;

    /* /////////////////////////////////////////////////////////////////////////
    // PROGRAM CORE  ///////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////*/
    cout<<"PROBLEM ID: "<<problemName<<"\n";

    //  GLOBAL INITIALIZATION ///////////////////////////////////////////////

    // initializing DEM parameters
    dem.discreteElementInit();

    printUfo(command_line, lbmCfgFile);

    // setting time
    dem.demTime=0.0;

    // initializing input-output
    io.initialize();

    // CYCLE /////////////////////////////
    // integrate in time
    bool run=true;
    while (run) {
        if (io.realTime!=0.0 && io.realTime>io.maxTime) {
            run=false;
        }
        // exit normally if the maximum simulation time has been reached
        else if (io.maximumTimeSteps && io.currentTimeStep>=io.maximumTimeSteps) {
            run=false;
        }
        else {
            // core of the code, performs time steps
            goCycle(io, dem);
        }
    }
  return 0;
}
