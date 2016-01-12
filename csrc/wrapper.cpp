/*
 *      LES by UNITS + DEM particles by Alessandro Leonardi
 *
 *      authors: Alessandro Leonardi & more
 *
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <csignal>
#include <stdio.h>
#include <string>
#include <cstring>
#include <mpi.h>


#include "getpot.h"

#include "IO.h"
#include "DEM.h"

using namespace std;

ProblemName problemName=NONE;

// stuff from Fortran
extern "C" {
void les_initialize();
void les_core();
void les_finalize();
}

extern int nproc,myid;

void goCycle(IO& io, DEM& dem) {

	// advance one step in time
	io.outputStep(dem.demTime, dem);

	++io.currentTimeStep;
	dem.demTime+=dem.deltat;

	dem.discreteElementStep();

}

void parseCommandLine(IO& io, GetPot& command_line) {

	extern int myid;

	if (myid==0) {
		cout<<"Parsing command line"<<endl;
	}

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
	}

	// lbm configuration file
	io.configFile = "config.cfg";
	if (command_line.search("-c")) {
		io.configFile = command_line.next(io.configFile.c_str());
		cout<<"Using "<<io.configFile<<" for DEM parameters\n";
		// copy the lbm configuration file into the data folder for later reproducibility
		std::system((std::string("mkdir -p '") + io.workDirectory + "'; cp '" + io.configFile + "' '" + io.workDirectory + "'").c_str());
	}

	// make sure the config files can be read
	std::ifstream lbmCfgFile(io.configFile.c_str());
	if (!lbmCfgFile) {
		cout << "ERROR: Can't open config file \"" << io.configFile << "\" for reading!\n";
		//        return ERROR;
	}
	lbmCfgFile.close();
}

void printUfo(GetPot& command_line, GetPot& lbmCfgFile){

	// warn about unused parameters

	std::vector<std::string> ufos=lbmCfgFile.unidentified_variables();
	if (ufos.size() > 0) {
		cout << INFO << "WARNING: Unrecognized config file parameter(s):";
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

void parseConfigFile(IO& io, DEM& dem, GetPot& lbmCfgFile, GetPot& command_line) {

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

void parseLesConfigFile(GetPot& command_line, GetPot& lesCfgFile) {

	if (myid==0) {
		cout<<"Parsing LES configuration file"<<endl;
	}

	// Grid size
	PARSE_EXTERNAL(lesCfgFile, int, n1, "n1", -1);
	ASSERT(n1 > 0);
	PARSE_EXTERNAL(lesCfgFile, int, n2, "n2", -1);
	ASSERT(n2 > 0);
	PARSE_EXTERNAL(lesCfgFile, int, n3, "n3", -1);
	ASSERT(n3 > 0);

	// Number of scalar equations
	PARSE_EXTERNAL(lesCfgFile, int, nscal, "nscal", -1);
	ASSERT(nscal > 0);

	// AGENERALE ---------------------------------------------------------
	// Iteration, restart, print
	PARSE_EXTERNAL(lesCfgFile, int, niter, "niter", -1);
	ASSERT(niter >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, i_rest, "i_rest", -1);
	ASSERT(i_rest >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, print_iter_or_time, "print_iter_or_time", -1);
	ASSERT(print_iter_or_time >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, i_print, "i_print", -1);
	ASSERT(i_print >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, i_time, "i_time", -1);
	ASSERT(i_time >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, i_printfile, "i_printfile", -1);
	ASSERT(i_printfile >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, i_paraview, "i_paraview", -1);
	ASSERT(i_paraview >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, i_cumulative, "i_cumulative", -1);
	ASSERT(i_cumulative >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, iformat_newres, "iformat_newres", -1);
	ASSERT(iformat_newres >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, iformat_grid, "iformat_grid", -1);
	ASSERT(iformat_grid >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ifolder, "ifolder", -1);
	ASSERT(ifolder >= 0);
	//string string_newres_format_C,string_grid_format_C;
	//PARSE_EXTERNAL(lesCfgFile, string, string_newres_format_C, "string_newres_format", "aa");
	//PARSE_EXTERNAL(lesCfgFile, string, string_grid_format_C, "string_grid_format", "aa");

	//char char_newres_format_C[12]=c_str(string_newres_format_C);
	//extern int string_newres_format( char string[12] );
	//convertToFortran(string_newres_format, sizeof string_newres_format, char_newres_format_C);
	//ASSERT(string_newres_format >= 0);

	//extern int niter;
	//PARSE_CLASS_MEMBER(lesCfgFile, string_grid_format, "string_grid_format", 0);
	//ASSERT(string_grid_format >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, inf, "inf", -1);
	ASSERT(inf >= 0);

	// Characteristic numbers
	PARSE_EXTERNAL(lesCfgFile, double, re, "re", -1.0);
	ASSERT(re >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, rich, "rich", -1.0);
	ASSERT(rich >= 0.0);

	// Characteristic numbers
	PARSE_EXTERNAL(lesCfgFile, int, ind_cou, "ind_cou", -1);
	ASSERT(ind_cou >= 0);
	PARSE_EXTERNAL(lesCfgFile, double, cou, "cou", -1.0);
	ASSERT(cou >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, dt, "dt", -1.0);
	ASSERT(dt >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, int, espl, "espl", -1);
	ASSERT(espl >= 0);

	// Numerical schemes, equations to solve
	PARSE_EXTERNAL(lesCfgFile, int, attiva_scal, "attiva_scal", -1);
	ASSERT(attiva_scal >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, potenziale, "potenziale", -1);
	ASSERT(potenziale >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, insc, "insc", -1);
	ASSERT(insc >= 0);

	// Pressure condition, pressure solver
	PARSE_EXTERNAL(lesCfgFile, double, omega, "om_pre", -1.0);
	ASSERT(omega >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, int, jpos, "jpos", -1);
	ASSERT(jpos >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, islor, "islor", -1);
	ASSERT(islor >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, nlevmultimax, "nlevmultimax", -1);
	ASSERT(nlevmultimax >= 0);
	cout<<"nlevmultimax = "<<nlevmultimax<<endl;
	PARSE_EXTERNAL(lesCfgFile, double, bbx, "bbx", -1.0);
	ASSERT(bbx >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, bby, "bby", -1.0);
	ASSERT(bby >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, bbz, "bbz", -1.0);
	ASSERT(bbz >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, eps, "eps", -1.0);
	ASSERT(eps >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, int, ficycle, "ficycle", -1);
	ASSERT(ficycle >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, bodypressure, "bodypressure", -1);
	ASSERT(bodypressure >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ipress_cart, "ipress_cart", -1);
	ASSERT(ipress_cart >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, freesurface, "freesurface", -1);
	ASSERT(freesurface >= 0);

	// Buondary condition
	PARSE_EXTERNAL(lesCfgFile, int, lett, "lett", -1);
	ASSERT(lett >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ibb, "ibb", -1);
	ASSERT(ibb >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, bodyforce, "bodyforce", -1);
	ASSERT(bodyforce >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, num_iter, "num_iter", -1);
	ASSERT(num_iter >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, coef_wall, "coef_wall", -1);
	ASSERT(coef_wall >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, integrale, "integrale", -1);
	ASSERT(integrale >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, rough, "rough", -1);
	ASSERT(rough >= 0);
	PARSE_EXTERNAL(lesCfgFile, double, z0, "z0", -1.0);
	ASSERT(z0 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, int, att_wm_sgs, "att_wm_sgs", -1);
	ASSERT(att_wm_sgs >= 0);

	// Turbulence model
	PARSE_EXTERNAL(lesCfgFile, int, nsgs, "nsgs", -1);
	ASSERT(nsgs >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, inmod, "inmod", -1);
	ASSERT(inmod >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, inmodrho, "inmodrho", -1);
	ASSERT(inmodrho >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, isotropo, "isotropo", -1);
	ASSERT(isotropo >= 0);
	PARSE_EXTERNAL(lesCfgFile, double, cost, "cost", -1.0);
	ASSERT(cost >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, costh, "costH", -1.0);
	ASSERT(costh >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, costv, "costV", -1.0);
	ASSERT(costv >= 0.0);

	// Forcing
	PARSE_EXTERNAL(lesCfgFile, int, indm, "indm", -1);
	ASSERT(indm >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, windyes, "windyes", -1);
	ASSERT(windyes >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, wavebk, "wavebk", -1);
	ASSERT(wavebk >= 0);
	PARSE_EXTERNAL(lesCfgFile, double, alpha, "alpha", -1.0);
	ASSERT(alpha >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, c10, "c10", -1.0);
	ASSERT(c10 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, l_0, "l_0", -1.0);
	ASSERT(l_0 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, int, langyes, "langyes", -1);
	ASSERT(langyes >= 0);
	PARSE_EXTERNAL(lesCfgFile, double, lamb, "lamb", -1.0);
	ASSERT(lamb >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, h_0, "h_0", -1.0);
	ASSERT(h_0 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, a1, "A1", -1.0);
	ASSERT(a1 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, a2, "A2", -1.0);
	ASSERT(a2 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, a3, "A3", -1.0);
	ASSERT(a3 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, u0, "U0", -1.0);
	ASSERT(u0 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, v0, "V0", -1.0);
	ASSERT(v0 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, w0, "W0", -1.0);
	ASSERT(w0 >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, latitude, "LATITUDE", -1.0);
	ASSERT(latitude >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, omegam2, "omegaM2", -1.0);
	ASSERT(omegam2 >= 0.0);

	// Others
	PARSE_EXTERNAL(lesCfgFile, int, visualizzo, "visualizzo", -1);
	ASSERT(visualizzo >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, lagr, "lagr", -1);
	ASSERT(lagr >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, i_sta, "i_sta", -1);
	ASSERT(i_sta >= 0);

	// Movies
	PARSE_EXTERNAL(lesCfgFile, int, imovieprint, "imovie", -1);
	ASSERT(imovieprint >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ktime_movie, "ktime_movie", -1);
	ASSERT(ktime_movie >= 0);
	PARSE_EXTERNAL(lesCfgFile, double, dt_movie, "dt_movie", -1.0);
	ASSERT(dt_movie >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, double, dt_delay, "dt_delay", -1.0);
	ASSERT(dt_delay >= 0.0);
	PARSE_EXTERNAL(lesCfgFile, int, i_movie, "i_movie", -1);
	ASSERT(i_movie >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, j_movie, "j_movie", -1);
	ASSERT(j_movie >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, k_movie, "k_movie", -1);
	ASSERT(k_movie >= 0);

	// ABOUNDARY ---------------------------------------------------------
	// Periodicity
	PARSE_EXTERNAL(lesCfgFile, int, ip, "ip", -1);
	ASSERT(ip >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, jp, "jp", -1);
	ASSERT(jp >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, kp, "kp", -1);
	ASSERT(kp >= 0);

	// For open boundary, if periodicity on, these are not used!
	PARSE_EXTERNAL(lesCfgFile, int, infout1, "infout1", -1);
	ASSERT(infout1 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, infout2, "infout2", -1);
	ASSERT(infout2 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, infout3, "infout3", -1);
	ASSERT(infout3 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, infout4, "infout4", -1);
	ASSERT(infout4 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, infout5, "infout5", -1);
	ASSERT(infout5 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, infout6, "infout6", -1);
	ASSERT(infout6 >= 0);

	// Settings for solid wall, only if infout=2
	PARSE_EXTERNAL(lesCfgFile, int, iboun1, "iboun1", -1);
	ASSERT(iboun1 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, iboun2, "iboun2", -1);
	ASSERT(iboun2 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, iboun3, "iboun3", -1);
	ASSERT(iboun3 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, iboun4, "iboun4", -1);
	ASSERT(iboun4 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, iboun5, "iboun5", -1);
	ASSERT(iboun5 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, iboun6, "iboun6", -1);
	ASSERT(iboun6 >= 0);

	// wall function at solid wall, only wfp3 and wfp4 are currently used
	PARSE_EXTERNAL(lesCfgFile, int, wfp1, "wfp1", -1);
	ASSERT(wfp1 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, wfp2, "wfp2", -1);
	ASSERT(wfp2 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, wfp3, "wfp3", -1);
	ASSERT(wfp3 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, wfp4, "wfp4", -1);
	ASSERT(wfp4 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, wfp5, "wfp5", -1);
	ASSERT(wfp5 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, wfp6, "wfp6", -1);
	ASSERT(wfp6 >= 0);

	// bodyforce in buffer region to enhance turbulence, if ibb=1 in Agenerale.in
	PARSE_EXTERNAL(lesCfgFile, int, ibodybuffer1, "ibodybuffer1", -1);
	ASSERT(ibodybuffer1 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ibodybuffer2, "ibodybuffer2", -1);
	ASSERT(ibodybuffer2 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ibodybuffer3, "ibodybuffer3", -1);
	ASSERT(ibodybuffer3 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ibodybuffer4, "ibodybuffer4", -1);
	ASSERT(ibodybuffer4 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ibodybuffer5, "ibodybuffer5", -1);
	ASSERT(ibodybuffer5 >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, ibodybuffer6, "ibodybuffer6", -1);
	ASSERT(ibodybuffer6 >= 0);

	if(ibb == 0) {
		ASSERT(ibodybuffer1 == 0);
		ASSERT(ibodybuffer2 == 0);
		ASSERT(ibodybuffer3 == 0);
		ASSERT(ibodybuffer4 == 0);
		ASSERT(ibodybuffer5 == 0);
		ASSERT(ibodybuffer6 == 0);
	}

	PARSE_EXTERNAL(lesCfgFile, int, ispon, "ispon", -1);
	ASSERT(ispon >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, kspon, "kspon", -1);
	ASSERT(kspon >= 0);
	PARSE_EXTERNAL(lesCfgFile, int, corr_factor, "corr_factor", -1);
	ASSERT(corr_factor >= 0);

	// ASCALARE  ---------------------------------------------------------
	PARSE_EXTERNAL(lesCfgFile, int, re_analogy, "re_analogy", -1);
	ASSERT(re_analogy >= 0);
	PARSE_EXTERNAL_ARRAY(lesCfgFile, double, nscal, c_pran, "pran", -1.0);
	PARSE_EXTERNAL_ARRAY(lesCfgFile, double, nscal, c_prsc, "prsc", -1.0);

	PARSE_EXTERNAL(lesCfgFile, int, imoist, "imoist", -1);
	ASSERT(imoist >= 0);
	PARSE_EXTERNAL(lesCfgFile, double, tref, "Tref", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, qref, "Qref", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, betat, "betaT", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, betaq, "betaQ", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, ma, "Ma", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, mv, "Mv", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, lv, "Lv", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, gdry, "Gdry", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, rd, "Rd", -1.0);
	PARSE_EXTERNAL(lesCfgFile, double, cpd, "cpd", -1.0);

	// APIANISONDE  ---------------------------------------------------------
	PARSE_EXTERNAL(lesCfgFile, int, npiani, "npiani", -1);
	ASSERT(npiani >= 0);
	PARSE_EXTERNAL_ARRAY(lesCfgFile, int, npiani, c_piani, "piani", -1);

	PARSE_EXTERNAL(lesCfgFile, int, nsonde, "nsonde", -1);
	ASSERT(nsonde >= 0);
	PARSE_EXTERNAL_ARRAY(lesCfgFile, int, nsonde, c_sondeindexi, "indexi", -1);
	PARSE_EXTERNAL_ARRAY(lesCfgFile, int, nsonde, c_sondeindexj, "indexj", -1);
	PARSE_EXTERNAL_ARRAY(lesCfgFile, int, nsonde, c_sondeindexk, "indexk", -1);

	// AFILTRAGGIO  ---------------------------------------------------------

	PARSE_EXTERNAL(lesCfgFile, int, ifiltro, "ifiltro", -1);
	PARSE_EXTERNAL(lesCfgFile, int, nfiltro, "nfiltro", -1);

	PARSE_EXTERNAL(lesCfgFile, int, filtrou, "filtrou", -1);
	PARSE_EXTERNAL(lesCfgFile, int, filtrov, "filtrov", -1);
	PARSE_EXTERNAL(lesCfgFile, int, filtrow, "filtrow", -1);
	PARSE_EXTERNAL(lesCfgFile, int, filtrorho, "filtrorho", -1);
	PARSE_EXTERNAL(lesCfgFile, int, filtrofi, "filtrofi", -1);

	PARSE_EXTERNAL(lesCfgFile, int, xstart, "xstart", -1);
	PARSE_EXTERNAL(lesCfgFile, int, xend, "xend", -1);
	PARSE_EXTERNAL(lesCfgFile, int, ystart, "ystart", -1);
	PARSE_EXTERNAL(lesCfgFile, int, yend, "yend", -1);
	PARSE_EXTERNAL(lesCfgFile, int, zstart, "zstart", -1);
	PARSE_EXTERNAL(lesCfgFile, int, zend, "zend", -1);


	MPI::COMM_WORLD.Barrier();

	//	int dim=3;
	//	extern int** c_sonde;
	//	c_sonde = new int* [nsonde];
	//	for (unsigned int n=0; n<nsonde; ++n) {
	//		c_sonde[n]=new int [3];
	//	}


	//	for (unsigned int sondeIndex=0; sondeIndex<nsonde; ++sondeIndex) {
	//		const unsigned int correctSondeIndex=sondeIndex+1;
	//		ostringstream convert;
	//		convert<<correctSondeIndex;
	//		const string stringIndex = convert.str();
	//		string stringName[3];
	//		stringName[0] = stringIndex + "/indexi";
	//		stringName[1] = stringIndex + "/indexj";
	//		stringName[2] = stringIndex + "/indexk";
	//		for (unsigned int dimIndex=0; dimIndex<3; ++dimIndex) {
	//			// names to point the right place in configuration file and command line
	//			const string stringNameHere=stringName[dimIndex];
	//			const string optionName = "-" + stringNameHere;
	//
	//			if (!lesCfgFile.have_variable(stringNameHere) && !command_line.search(optionName)) {
	//				cout << INFO << "WARNING: Parameter " << stringNameHere << " is missing! Using default: " << -1 << std::endl;
	//			}
	//			c_sonde[sondeIndex][dimIndex]=lesCfgFile(stringNameHere, -1);;
	//			//<<c_sonde[sondeIndex,dimIndex]<<"/("<<sondeIndex<<" "<<dimIndex<<") = "<<lesCfgFile(stringNameHere, -1)<<endl;
	//			if (command_line.search(optionName)) {
	//				c_sonde[sondeIndex][dimIndex] = command_line.next(c_sonde[sondeIndex][dimIndex]);
	//			}
	//		}
	//	}

	if (myid==0) {
		cout<<"Number of piani :"<<npiani<<endl;
		for (unsigned int n=0; n<npiani; ++n) {
			cout<<"Piano number "<<n<<", position: "<<c_piani[n]<<endl;
		}
	}
	extern int myid;
	if (myid==0) {
		cout<<"Number of sonde :"<<nsonde<<endl;
		for (unsigned int n=0; n<nsonde; ++n) {
			cout<<"Sonda number "<<n<<", position: ("<<c_sondeindexi[n]<<" "<<c_sondeindexj[n]<<" "<<c_sondeindexk[n]<<")"<<endl;
		}
	}

	if (myid==0) {
		printUfo(command_line, lesCfgFile);
		cout<<"LES parsing finished"<<endl;
	}

}

int main(int argc, char* argv[]){

	extern int niter,ktime;
	extern bool termina;

	MPI::Status status;

	MPI::Init();

	// DECLARATION OF VARIABLES - Input-Output ///////////////
	IO io;

	// DECLARATION OF VARIABLES - DEM ///////////////
	DEM dem;

	// parsing command line and configuration file
	GetPot command_line(argc, argv);
	parseCommandLine(io, command_line);

	// parsing LBM input file
	GetPot lbmCfgFile(io.configFile);
	parseConfigFile(io, dem, lbmCfgFile, command_line);

	myid = MPI::COMM_WORLD.Get_rank();
	nproc = MPI::COMM_WORLD.Get_size();

	string lesConfigFile="lesCgfFile.cfg";
	GetPot lesCfgFile(lesConfigFile);
	parseLesConfigFile(command_line,lesCfgFile);

	MPI::COMM_WORLD.Barrier();

	les_initialize();

	for (ktime=1; ktime<=niter; ktime++) {
		// nesting
		if (termina) {
			break;
		}
		les_core();
	}

	les_finalize();

	MPI::Finalize();

	if (false) {



		io.currentTimeStep = 0;

    	/* /////////////////////////////////////////////////////////
        // PROGRAM CORE  ///////////////////////////////////////////
        //////////////////////////////////////////////////////////*/
		cout<<"PROBLEM ID: "<<problemName<<"\n";

		//  GLOBAL INITIALIZATION //////////////////////////////////

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
	}
	return 0;


}
