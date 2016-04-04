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
#include <time.h>


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
void get_cpp_strings(const char*,const char*,const char*);
void les_domain_size();
}

extern int nproc,myid;

void parseCommandLine(IO& io,GetPot& command_line) {

	// std::system((std::string("mkdir -p '") + io.workDirectory + "'; cp '" + io.lbmCfgName + "' '" + io.workDirectory + "'").c_str());
	extern int myid;

	// PART 1 - SIMULATION FOLDER ---------------------------------------------------

	// create the data directory in case it doesn't exist
	// standard is ./results
	io.workDirectory = "./results";
	if (command_line.search("-d")) {
		io.workDirectory = command_line.next(io.workDirectory);
		if (myid==0) {
			std::system((std::string("mkdir -p '") + io.workDirectory + "'").c_str());
		}
	}

	// if the name of the simulation is specified create subfolder
	// standard is just a string with the current time
	std::string programName = "time";
	if (command_line.search("-n")) {
		// creating work directory
		programName=command_line.next(programName);
	}

	// if there is no specified name,create a folder named with the initial time
	if (programName=="time") {
		char* buff;
		if (myid==0) {
			time_t now = time(NULL);
			const struct tm* timeptr=localtime(&now);
			strftime(buff,20,"%Y%m%d_%H%M%S",timeptr);
			//time_t rawtime;
			//  time (&rawtime);
		}
		MPI_Bcast(buff,20,MPI_CHAR,0,MPI::COMM_WORLD);
		MPI::COMM_WORLD.Barrier();
		io.workDirectory=io.workDirectory+"/"+std::string(buff);
	}
	// otherwise use the name as input
	else {
		io.workDirectory=io.workDirectory+"/"+programName;
	}

	// make simulation folder
	if (myid==0) {
		std::system((std::string("mkdir -p '") + io.workDirectory + "'").c_str());
	}

	if (myid==0) {
		cout<<" Simulation folder "<<io.workDirectory<<"\n";
	}

	// PART 2 - CONFIGURATION FILES ---------------------------------------------------

	// DEM configuration file
	io.demConfigFile = "demConfigFile.cfg";
	if (command_line.search("-dem")) {
		io.demConfigFile = command_line.next(io.demConfigFile.c_str());
	}
	if (myid==0) {
		cout<<" Using '"<<io.demConfigFile<<"' for DEM parameters\n";
	}
	// LES configuration file
	io.lesConfigFile = "lesConfigFile.cfg";
	if (command_line.search("-les")) {
		io.lesConfigFile = command_line.next(io.lesConfigFile.c_str());
	}
	if (myid==0) {
		cout<<" Using '"<<io.lesConfigFile<<"' for LES parameters\n";
	}
	// copy the configuration files into the data folder for later reproducibility
	if (myid==0) {
		std::system(("cp '" + io.demConfigFile + "' '" + io.workDirectory + "'").c_str());
		std::system(("cp '" + io.lesConfigFile + "' '" + io.workDirectory + "'").c_str());
	}
	// LES grid file
	io.lesGridFile = "gri3dp_in.dat";
	if (command_line.search("-grid")) {
		io.lesGridFile = command_line.next(io.lesGridFile.c_str());
	}
	if (myid==0) {
		cout<<" Using grid file '"<<io.lesGridFile<<"' for LES\n";
	}

	// make sure the config files can be read
	std::ifstream demConfigStream(io.demConfigFile.c_str());
	std::ifstream lesConfigStream(io.lesConfigFile.c_str());
	if (!demConfigStream) {
		if (myid==0) {
			cout << " ERROR: Can't open DEM config file \"" << io.demConfigFile << "\" for reading!\n";
		}
		MPI::Finalize();
		exit(0);
	}
	if (!lesConfigStream) {
		if (myid==0) {
			cout << " ERROR: Can't open LES config file \"" << io.lesConfigFile << "\" for reading!\n";
		}
		MPI::Finalize();
		exit(0);
	}
	demConfigStream.close();
	lesConfigStream.close();

	// LES restart file
	io.lesRestartFile = "old_res_form";
	if (command_line.search("-restart")) {
		io.lesRestartFile = command_line.next(io.lesRestartFile.c_str());
		if (myid==0) {
			cout<<" Using restart file '"<<io.lesRestartFile<<"' for LES\n";
		}
	}

	// PART 3 - SEND TO FORTAN ---------------------------------------------------

	get_cpp_strings(io.workDirectory.c_str(),io.lesGridFile.c_str(),io.lesRestartFile.c_str());

}

void printUfo(GetPot& command_line,GetPot& lbmCfgFile){

	// warn about unused parameters

	std::vector<std::string> ufos=lbmCfgFile.unidentified_variables();
	if (ufos.size() > 0) {
		cout << " WARNING: Unrecognized config file parameter(s):"; //<< INFO
		for (unsigned int i = 0; i < ufos.size(); ++i) {
			cout << " " << ufos[i];
		}
		cout << endl;
	}

	ufos = command_line.unidentified_arguments();
	if (ufos.size() > 0) {
		cout << " WARNING: Unrecognized command line argument(s):"; //<< INFO
		for (unsigned int i = 0; i < ufos.size(); ++i) {
			cout << " " << ufos[i];
		}
		cout << endl;
	}
}

void parseDemConfigFile(IO& io,DEM& dem,GetPot& demCfgFile,GetPot& command_line) {

	// PROBLEM NAME //////////////
	// necessary for hard coded sections of the code
	string problemNameString;
	PARSE_CLASS_MEMBER(demCfgFile,problemNameString,"problemName","none");
	if (problemNameString=="DRUM") problemName=DRUM;
	else if (problemNameString=="SHEARCELL") problemName=SHEARCELL;
	else if (problemNameString=="NONE") problemName=NONE;
	else if (problemNameString=="AVALANCHE") problemName=AVALANCHE;
	else if (problemNameString=="SPLASH") problemName=SPLASH;
	else if (problemNameString=="BOX") problemName=BOX;
	else if (problemNameString=="NET") problemName=NET;
	else if (problemNameString=="DIFF") problemName=DIFF;

	// GETTING SIMULATION PARAMETERS  /////////
	// stuff taken from shell code
	PARSE_CLASS_MEMBER(demCfgFile,io.screenExpTime,"screenExpTime",0.0);
	ASSERT(io.screenExpTime>= 0);
	PARSE_CLASS_MEMBER(demCfgFile,io.partExpTime,"partExpTime",0.0);
	ASSERT(io.partExpTime>= 0);
	PARSE_CLASS_MEMBER(demCfgFile,io.objectExpTime,"objectExpTime",0.0);
	ASSERT(io.objectExpTime>= 0);
	PARSE_CLASS_MEMBER(demCfgFile,io.recycleExpTime,"recycleExpTime",0.0);
	ASSERT(io.recycleExpTime>= 0);

	// getting DEM parameters and particle initial status
	dem.discreteElementGet(demCfgFile,command_line);

	// problem specific parameters
	switch (problemName) {
	case DRUM: {
		// drum rotational speed
		PARSE_CLASS_MEMBER(demCfgFile,dem.drumSpeed,"drumSpeed",0.0);
		break;
	}
	}
}

void parseLesConfigFile(GetPot& command_line,GetPot& lesCfgFile) {

	if (myid==0) {
		cout<<" Parsing LES configuration file"<<endl;
	}

	// AGENERALE ---------------------------------------------------------
	// Iteration,restart,print
	PARSE_EXTERNAL(lesCfgFile,int,niter,"niter",-1);
	ASSERT(niter>=0);
	PARSE_EXTERNAL(lesCfgFile,int,i_rest,"i_rest",-1);
	ASSERT(i_rest>=0);
	PARSE_EXTERNAL(lesCfgFile,int,print_iter_or_time,"print_iter_or_time",-1);
	ASSERT(print_iter_or_time>=0);
	ASSERT(print_iter_or_time<=1);
	PARSE_EXTERNAL(lesCfgFile,int,i_print,"i_print",-1);
	ASSERT(i_print>=0);
	PARSE_EXTERNAL(lesCfgFile,int,i_time,"i_time",-1);
	ASSERT(i_time>=0);
	PARSE_EXTERNAL(lesCfgFile,int,i_printfile,"i_printfile",-1);
	ASSERT(i_printfile>=0);
	ASSERT(i_printfile<=2);
	PARSE_EXTERNAL(lesCfgFile,bool,i_paraview,"i_paraview",false);
	PARSE_EXTERNAL(lesCfgFile,bool,i_newres,"i_newres",false);
	PARSE_EXTERNAL(lesCfgFile,bool,i_medietempo,"i_medietempo",false);
	PARSE_EXTERNAL(lesCfgFile,bool,i_cumulative,"i_cumulative",false);
	PARSE_EXTERNAL(lesCfgFile,int,iformat_newres,"iformat_newres",-1);
	ASSERT(iformat_newres>=0);
	ASSERT(iformat_newres<=2);
	PARSE_EXTERNAL(lesCfgFile,bool,iformat_grid,"iformat_grid",false);

	PARSE_EXTERNAL(lesCfgFile,int,inf,"inf",-1);
	ASSERT(inf>=0);

	// Characteristic numbers
	PARSE_EXTERNAL(lesCfgFile,double,re,"re",-1.0);
	ASSERT(re>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,rich,"rich",-1.0);
	ASSERT(rich>=0.0);

	// Characteristic numbers
	PARSE_EXTERNAL(lesCfgFile,int,ind_cou,"ind_cou",-1);
	ASSERT(ind_cou>=0);
	PARSE_EXTERNAL(lesCfgFile,double,cou,"cou",-1.0);
	ASSERT(cou>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,dt,"dt",-1.0);
	ASSERT(dt>=0.0);
	PARSE_EXTERNAL(lesCfgFile,int,espl,"espl",-1);
	ASSERT(espl>=0);

	// Numerical schemes,equations to solve
	PARSE_EXTERNAL(lesCfgFile,int,attiva_scal,"attiva_scal",-1);
	ASSERT(attiva_scal>=0);
	PARSE_EXTERNAL(lesCfgFile,int,potenziale,"potenziale",-1);
	ASSERT(potenziale>=0);
	PARSE_EXTERNAL(lesCfgFile,int,insc,"insc",-1);
	ASSERT(insc>=0);

	// Pressure condition,pressure solver
	PARSE_EXTERNAL(lesCfgFile,double,omega,"om_pre",-1.0);
	ASSERT(omega>=0.0);
	PARSE_EXTERNAL(lesCfgFile,int,jpos,"jpos",-1);
	ASSERT(jpos>=0);
	PARSE_EXTERNAL(lesCfgFile,int,islor,"islor",-1);
	ASSERT(islor>=0);
	PARSE_EXTERNAL(lesCfgFile,int,nlevmultimax,"nlevmultimax",-1);
	ASSERT(nlevmultimax>=0);
	PARSE_EXTERNAL(lesCfgFile,double,bbx,"bbx",-1.0);
	ASSERT(bbx>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,bby,"bby",-1.0);
	ASSERT(bby>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,bbz,"bbz",-1.0);
	ASSERT(bbz>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,eps,"eps",-1.0);
	ASSERT(eps>=0.0);
	PARSE_EXTERNAL(lesCfgFile,int,ficycle,"ficycle",-1);
	ASSERT(ficycle>=0);
	PARSE_EXTERNAL(lesCfgFile,int,bodypressure,"bodypressure",-1);
	ASSERT(bodypressure>=0);
	PARSE_EXTERNAL(lesCfgFile,int,ipress_cart,"ipress_cart",-1);
	ASSERT(ipress_cart>=0);
	PARSE_EXTERNAL(lesCfgFile,int,freesurface,"freesurface",-1);
	ASSERT(freesurface>=0);

	// Buondary condition
	PARSE_EXTERNAL(lesCfgFile,int,lett,"lett",-1);
	ASSERT(lett>=0);
	PARSE_EXTERNAL(lesCfgFile,int,ibb,"ibb",-1);
	ASSERT(ibb>=0);
	PARSE_EXTERNAL(lesCfgFile,bool,bodyforce,"bodyforce",-1);
	PARSE_EXTERNAL(lesCfgFile,bool,particles,"particles",false);
	PARSE_EXTERNAL(lesCfgFile,int,num_iter,"num_iter",-1);
	ASSERT(num_iter>=0);
	PARSE_EXTERNAL(lesCfgFile,int,coef_wall,"coef_wall",-1);
	ASSERT(coef_wall>=0);
	PARSE_EXTERNAL(lesCfgFile,int,integrale,"integrale",-1);
	ASSERT(integrale>=0);
	PARSE_EXTERNAL(lesCfgFile,int,rough,"rough",-1);
	ASSERT(rough>=0);
	PARSE_EXTERNAL(lesCfgFile,double,z0,"z0",-1.0);
	ASSERT(z0>=0.0);
	PARSE_EXTERNAL(lesCfgFile,int,att_wm_sgs,"att_wm_sgs",-1);
	ASSERT(att_wm_sgs>=0);

	// Turbulence model
	PARSE_EXTERNAL(lesCfgFile,int,nsgs,"nsgs",-1);
	ASSERT(nsgs>=0);
	PARSE_EXTERNAL(lesCfgFile,int,inmod,"inmod",-1);
	ASSERT(inmod>=0);
	PARSE_EXTERNAL(lesCfgFile,int,inmodrho,"inmodrho",-1);
	ASSERT(inmodrho>=0);
	PARSE_EXTERNAL(lesCfgFile,int,isotropo,"isotropo",-1);
	ASSERT(isotropo>=0);
	PARSE_EXTERNAL(lesCfgFile,double,cost,"cost",-1.0);
	ASSERT(cost>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,costh,"costH",-1.0);
	ASSERT(costh>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,costv,"costV",-1.0);
	ASSERT(costv>=0.0);

	// Forcing
	PARSE_EXTERNAL(lesCfgFile,int,windyes,"windyes",-1);
	ASSERT(windyes>=0);
	PARSE_EXTERNAL(lesCfgFile,int,wavebk,"wavebk",-1);
	ASSERT(wavebk>=0);
	PARSE_EXTERNAL(lesCfgFile,double,alpha,"alpha",-1.0);
	ASSERT(alpha>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,c10,"c10",-1.0);
	ASSERT(c10>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,l_0,"l_0",-1.0);
	ASSERT(l_0>=0.0);
	PARSE_EXTERNAL(lesCfgFile,int,langyes,"langyes",-1);
	ASSERT(langyes>=0);
	PARSE_EXTERNAL(lesCfgFile,double,lamb,"lamb",-1.0);
	ASSERT(lamb>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,h_0,"h_0",-1.0);
	ASSERT(h_0>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,latitude,"LATITUDE",-1.0);
	ASSERT(latitude>=0.0);

	// Others
	PARSE_EXTERNAL(lesCfgFile,int,visualizzo,"visualizzo",-1);
	ASSERT(visualizzo>=0);
	PARSE_EXTERNAL(lesCfgFile,int,lagr,"lagr",-1);
	ASSERT(lagr>=0);
	PARSE_EXTERNAL(lesCfgFile,int,i_sta,"i_sta",-1);
	ASSERT(i_sta>=0);

	// Movies
	PARSE_EXTERNAL(lesCfgFile,int,imovieprint,"imovie",-1);
	ASSERT(imovieprint>=0);
	PARSE_EXTERNAL(lesCfgFile,int,ktime_movie,"ktime_movie",-1);
	ASSERT(ktime_movie>=0);
	PARSE_EXTERNAL(lesCfgFile,double,dt_movie,"dt_movie",-1.0);
	ASSERT(dt_movie>=0.0);
	PARSE_EXTERNAL(lesCfgFile,double,dt_delay,"dt_delay",-1.0);
	ASSERT(dt_delay>=0.0);
	PARSE_EXTERNAL(lesCfgFile,int,i_movie,"i_movie",-1);
	ASSERT(i_movie>=0);
	PARSE_EXTERNAL(lesCfgFile,int,j_movie,"j_movie",-1);
	ASSERT(j_movie>=0);
	PARSE_EXTERNAL(lesCfgFile,int,k_movie,"k_movie",-1);
	ASSERT(k_movie>=0);

	// ABOUNDARY ---------------------------------------------------------
	// Periodicity
	PARSE_EXTERNAL(lesCfgFile,int,ip,"ip",-1);
	ASSERT(ip>=0);
	PARSE_EXTERNAL(lesCfgFile,int,jp,"jp",-1);
	ASSERT(jp>=0);
	PARSE_EXTERNAL(lesCfgFile,int,kp,"kp",-1);
	ASSERT(kp>=0);

	// For open boundary,if periodicity on,these are not used!
	PARSE_EXTERNAL(lesCfgFile,int,infout1,"infout1",-1);
	ASSERT(infout1>=0);
	PARSE_EXTERNAL(lesCfgFile,int,infout2,"infout2",-1);
	ASSERT(infout2>=0);
	PARSE_EXTERNAL(lesCfgFile,int,infout3,"infout3",-1);
	ASSERT(infout3>=0);
	PARSE_EXTERNAL(lesCfgFile,int,infout4,"infout4",-1);
	ASSERT(infout4>=0);
	PARSE_EXTERNAL(lesCfgFile,int,infout5,"infout5",-1);
	ASSERT(infout5>=0);
	PARSE_EXTERNAL(lesCfgFile,int,infout6,"infout6",-1);
	ASSERT(infout6>=0);

	// Settings for solid wall,only if infout=2
	PARSE_EXTERNAL(lesCfgFile,int,iboun1,"iboun1",-1);
	ASSERT(iboun1>=0);
	PARSE_EXTERNAL(lesCfgFile,int,iboun2,"iboun2",-1);
	ASSERT(iboun2>=0);
	PARSE_EXTERNAL(lesCfgFile,int,iboun3,"iboun3",-1);
	ASSERT(iboun3>=0);
	PARSE_EXTERNAL(lesCfgFile,int,iboun4,"iboun4",-1);
	ASSERT(iboun4>=0);
	PARSE_EXTERNAL(lesCfgFile,int,iboun5,"iboun5",-1);
	ASSERT(iboun5>=0);
	PARSE_EXTERNAL(lesCfgFile,int,iboun6,"iboun6",-1);
	ASSERT(iboun6>=0);

	// wall function at solid wall,only wfp3 and wfp4 are currently used
	PARSE_EXTERNAL(lesCfgFile,int,wfp1,"wfp1",-1);
	ASSERT(wfp1>=0);
	PARSE_EXTERNAL(lesCfgFile,int,wfp2,"wfp2",-1);
	ASSERT(wfp2>=0);
	PARSE_EXTERNAL(lesCfgFile,int,wfp3,"wfp3",-1);
	ASSERT(wfp3>=0);
	PARSE_EXTERNAL(lesCfgFile,int,wfp4,"wfp4",-1);
	ASSERT(wfp4>=0);
	PARSE_EXTERNAL(lesCfgFile,int,wfp5,"wfp5",-1);
	ASSERT(wfp5>=0);
	PARSE_EXTERNAL(lesCfgFile,int,wfp6,"wfp6",-1);
	ASSERT(wfp6>=0);

	// bodyforce in buffer region to enhance turbulence,if ibb=1 in Agenerale.in
	PARSE_EXTERNAL(lesCfgFile,int,ibodybuffer1,"ibodybuffer1",-1);
	ASSERT(ibodybuffer1>=0);
	PARSE_EXTERNAL(lesCfgFile,int,ibodybuffer2,"ibodybuffer2",-1);
	ASSERT(ibodybuffer2>=0);
	PARSE_EXTERNAL(lesCfgFile,int,ibodybuffer3,"ibodybuffer3",-1);
	ASSERT(ibodybuffer3>=0);
	PARSE_EXTERNAL(lesCfgFile,int,ibodybuffer4,"ibodybuffer4",-1);
	ASSERT(ibodybuffer4>=0);
	PARSE_EXTERNAL(lesCfgFile,int,ibodybuffer5,"ibodybuffer5",-1);
	ASSERT(ibodybuffer5>=0);
	PARSE_EXTERNAL(lesCfgFile,int,ibodybuffer6,"ibodybuffer6",-1);
	ASSERT(ibodybuffer6>=0);

	if(ibb == 0) {
		ASSERT(ibodybuffer1 == 0);
		ASSERT(ibodybuffer2 == 0);
		ASSERT(ibodybuffer3 == 0);
		ASSERT(ibodybuffer4 == 0);
		ASSERT(ibodybuffer5 == 0);
		ASSERT(ibodybuffer6 == 0);
	}

	PARSE_EXTERNAL(lesCfgFile,int,ispon,"ispon",-1);
	ASSERT(ispon>=0);
	PARSE_EXTERNAL(lesCfgFile,int,kspon,"kspon",-1);
	ASSERT(kspon>=0);
	PARSE_EXTERNAL(lesCfgFile,int,corr_factor,"corr_factor",-1);
	ASSERT(corr_factor>=0);

	// ASCALARE  ---------------------------------------------------------
	// Number of scalar equations
	PARSE_EXTERNAL(lesCfgFile,int,nscal,"nscal",-1);
	ASSERT(nscal>=1);
	PARSE_EXTERNAL(lesCfgFile,int,re_analogy,"re_analogy",-1);
	ASSERT(re_analogy>=0);
	PARSE_EXTERNAL_ARRAY(lesCfgFile,double,nscal,c_pran,"pran",-1.0);
	PARSE_EXTERNAL_ARRAY(lesCfgFile,double,nscal,c_prsc,"prsc",-1.0);

	// APIANISONDE  ---------------------------------------------------------
	PARSE_EXTERNAL(lesCfgFile,int,npiani,"npiani",-1);
	ASSERT(npiani>=0);
	PARSE_EXTERNAL_ARRAY(lesCfgFile,int,npiani,c_piani,"piani",-1);

	PARSE_EXTERNAL(lesCfgFile,int,nsonde,"nsonde",-1);
	ASSERT(nsonde>=0);
	PARSE_EXTERNAL_ARRAY(lesCfgFile,int,nsonde,c_sondeindexi,"indexi",-1);
	PARSE_EXTERNAL_ARRAY(lesCfgFile,int,nsonde,c_sondeindexj,"indexj",-1);
	PARSE_EXTERNAL_ARRAY(lesCfgFile,int,nsonde,c_sondeindexk,"indexk",-1);

	// AFILTRAGGIO  ---------------------------------------------------------

	PARSE_EXTERNAL(lesCfgFile,int,ifiltro,"ifiltro",-1);
	PARSE_EXTERNAL(lesCfgFile,int,nfiltro,"nfiltro",-1);

	PARSE_EXTERNAL(lesCfgFile,int,filtrou,"filtrou",-1);
	PARSE_EXTERNAL(lesCfgFile,int,filtrov,"filtrov",-1);
	PARSE_EXTERNAL(lesCfgFile,int,filtrow,"filtrow",-1);
	PARSE_EXTERNAL(lesCfgFile,int,filtrorho,"filtrorho",-1);
	PARSE_EXTERNAL(lesCfgFile,int,filtrofi,"filtrofi",-1);

	PARSE_EXTERNAL(lesCfgFile,int,xstart,"xstart",-1);
	PARSE_EXTERNAL(lesCfgFile,int,xend,"xend",-1);
	PARSE_EXTERNAL(lesCfgFile,int,ystart,"ystart",-1);
	PARSE_EXTERNAL(lesCfgFile,int,yend,"yend",-1);
	PARSE_EXTERNAL(lesCfgFile,int,zstart,"zstart",-1);
	PARSE_EXTERNAL(lesCfgFile,int,zend,"zend",-1);


	MPI::COMM_WORLD.Barrier();

	if (myid==0) {
		cout<<" Number of piani: "<<npiani<<endl;
		for (unsigned int n=0; n<npiani; ++n) {
			cout<<" Piano number "<<n<<",position: "<<c_piani[n]<<endl;
		}
	}
	extern int myid;
	if (myid==0) {
		cout<<" Number of sonde: "<<nsonde<<endl;
		for (unsigned int n=0; n<nsonde; ++n) {
			cout<<" Sonda number "<<n<<",position: ("<<c_sondeindexi[n]<<" "<<c_sondeindexj[n]<<" "<<c_sondeindexk[n]<<")"<<endl;
		}
	}

	if (myid==0) {
		printUfo(command_line,lesCfgFile);
	}

}

int main(int argc,char* argv[]){

	extern int niter,ktime;
	extern bool termina;

	MPI::Status status;

	MPI::Init();

	myid = MPI::COMM_WORLD.Get_rank();
	nproc = MPI::COMM_WORLD.Get_size();

	if (myid==0) {
		cout<<" -------------->PROGRAM BEGINS<---------------\n";
	}

	// DECLARATION OF VARIABLES - Input-Output ///////////////
	IO io;

	// DECLARATION OF VARIABLES - DEM ///////////////
	DEM dem;

	// parsing command line and configuration file
	GetPot command_line(argc,argv);
	parseCommandLine(io,command_line);

	// parsing DEM input file
	GetPot demCfgFile(io.demConfigFile);
	parseDemConfigFile(io,dem,demCfgFile,command_line);

	// parsing LES input file
	GetPot lesCfgFile(io.lesConfigFile);
	parseLesConfigFile(command_line,lesCfgFile);

	MPI::COMM_WORLD.Barrier();

	//  INITIALIZATION //////////////////////////////////
	// first step: read domain size from grid
	les_domain_size();

	ktime=0;

	// initializing input-output
	io.initialize();
	MPI::COMM_WORLD.Barrier();

	// initializing DEM parameters
	dem.discreteElementInit();
	MPI::COMM_WORLD.Barrier();
	if (myid==0) {
		printUfo(command_line,demCfgFile);
	}

	les_initialize();
	dem.discreteElementStep1();

	// CYCLE /////////////////////////////
	for (ktime=1; ktime<=niter; ktime++) {

		if (termina) {
			break;
		}

		MPI::COMM_WORLD.Barrier();

		// corrector
		dem.discreteElementStep2();
		// advance one step in time
		extern double dt;
		dem.demTime+=dt;
		if (myid==0) {
			cout<<" ~~~~~~~ Enter DEM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
			cout<<" DEM time = "<<dem.demTime<<endl;
			io.outputStep(dem.demTime,dem);
			++io.currentTimeStep;
		}
		// predictor
		dem.discreteElementStep1();

		MPI::COMM_WORLD.Barrier();
		if (myid==0) {
			cout<<" ~~~~~~~ Enter LES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		}
		les_core();

	}

	les_finalize();

	MPI::Finalize();

	return 0;

}
