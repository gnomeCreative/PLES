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
#include <limits>

#include "getpot.h"

#include "IO.h"
#include "DEM.h"

using namespace std;

ProblemName problemName=NONE;
bool lesSolve=true;
bool demSolve=true;

// stuff from Fortran
extern "C" {
void les_initialize_part1();
void les_initialize_part2();
void les_core();
void les_finalize();
void get_cpp_strings(const char*,const char*,const char*,const char*);
void les_domain_size();
}

extern int nproc,myid;

void parseCommandLine(IO& io, GetPot& command_line) {

	if (myid==0) {
		cout<<" Parsing comand line"<<endl;
	}

	// PART 1 - SIMULATION FOLDER ---------------------------------------------------
	// create the data directory in case it doesn't exist
	// standard is ./results
	io.workDirectory = "./results";
	if (command_line.search("-d")) {
		io.workDirectory = command_line.next(io.workDirectory);
		if (myid==0) {
			string createFolder=std::string("mkdir -p '") + io.workDirectory + "'";
			cout<<"Creating folder with "<<createFolder<<endl;
			std::system(createFolder.c_str());
		}
	}
	MPI::COMM_WORLD.Barrier();

	if (command_line.search("-noles")) {
		lesSolve=false;
		if (myid==0) {
			cout<<"<<<<<<<<<<<<  LES deactivated  >>>>>>>>>>>>"<<endl;
		}
	}
	else {
		lesSolve=true;
	}
	if (command_line.search("-nodem")) {
		demSolve=false;
		if (myid==0) {
			cout<<"<<<<<<<<<<<<  DEM deactivated  >>>>>>>>>>>>"<<endl;
		}
	}
	else {
		demSolve=true;
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
        ostringstream convertTime;
        convertTime<<std::setfill('0') << std::setw(4) <<(io.now->tm_year + 1900)
                <<std::setfill('0') << std::setw(2) <<(io.now->tm_mon+1)
                <<std::setfill('0') << std::setw(2) <<io.now->tm_mday
                <<"_"
                <<std::setfill('0') << std::setw(2) <<io.now->tm_hour
                <<std::setfill('0') << std::setw(2) <<io.now->tm_min
                <<std::setfill('0') << std::setw(2) <<io.now->tm_sec;
        string simulationTime=convertTime.str();
        io.workDirectory=io.workDirectory+"/"+simulationTime;
        io.workDirectory=io.workDirectory.substr(0,io.workDirectory.size());
    }
	// otherwise use the name as input
	else {
		io.workDirectory=io.workDirectory+"/"+programName;
	}

    // make simulation folder
	if (myid==0) {
		std::system((std::string("mkdir -p '") + io.workDirectory + "'").c_str());
		cout<<" Simulation folder "<<io.workDirectory<<"\n";
	}

	// PART 2 - CONFIGURATION FILES ---------------------------------------------------

	// DEM configuration file
	if (demSolve) {
		io.demConfigFileName = "demConfigFileName.cfg";
		if (command_line.search("-dem")) {
			io.demConfigFileName = command_line.next(io.demConfigFileName.c_str());
		}
		if (myid==0) {
			cout<<" Using '"<<io.demConfigFileName<<"' for DEM parameters\n";
			// copy the configuration files into the data folder for later reproducibility
			std::system(("cp '" + io.demConfigFileName + "' '" + io.workDirectory + "'").c_str());
		}
		// make sure the config files can be read
		std::ifstream demConfigStream(io.demConfigFileName.c_str());
		if (!demConfigStream) {
			if (myid==0) {
				cout << " ERROR: Can't open DEM config file \"" << io.demConfigFileName << "\" for reading!\n";
			}
			MPI::Finalize();
			exit(0);
		}
		demConfigStream.close();
	}

	// LES configuration file
	if (lesSolve) {
		io.lesConfigFileName = "lesConfigFileName.cfg";
		if (command_line.search("-les")) {
			io.lesConfigFileName = command_line.next(io.lesConfigFileName.c_str());
		}
		if (myid==0) {
			cout<<" Using '"<<io.lesConfigFileName<<"' for LES parameters\n";
			// copy the configuration files into the data folder for later reproducibility
			std::system(("cp '" + io.lesConfigFileName + "' '" + io.workDirectory + "'").c_str());
		}
		// make sure the config files can be read
		std::ifstream lesConfigStream(io.lesConfigFileName.c_str());

		if (!lesConfigStream) {
			if (myid==0) {
				cout << " ERROR: Can't open LES config file \"" << io.lesConfigFileName << "\" for reading!\n";
			}
			MPI::Finalize();
			exit(0);
		}
		lesConfigStream.close();

		// LES restart file
		io.lesGridFile = "gri3dp_in.dat";
		if (command_line.search("-grid")) {
			io.lesGridFile = command_line.next(io.lesGridFile.c_str());
			if (myid==0) {
				cout<<" Using grid file '"<<io.lesGridFile<<"' for LES\n";
			}
		}

		// LES restart file
		io.lesRestartFile = "old_res_form";
		if (command_line.search("-restart")) {
			io.lesRestartFile = command_line.next(io.lesRestartFile.c_str());
			if (myid==0) {
				cout<<" Using restart file '"<<io.lesRestartFile<<"' for LES\n";
			}
		}
	}


}

void printUfo(IO& io, GetPot& command_line){

	// warn about unused parameters

	std::vector<std::string> ufos;

	ufos = command_line.unidentified_arguments();
	if (ufos.size() > 0) {
		cout << " WARNING: Unrecognized command line argument(s):"; //<< INFO
		for (unsigned int i = 0; i < ufos.size(); ++i) {
			cout << " " << ufos[i];
		}
		cout << endl;
	}

	if (lesSolve) {
		ufos=io.lesCfgFile.unidentified_variables();
		if (ufos.size() > 0) {
			cout << " WARNING: Unrecognized LES config file parameter(s):"; //<< INFO
			for (unsigned int i = 0; i < ufos.size(); ++i) {
				cout << " " << ufos[i];
			}
			cout << endl;
		}
	}

	if (demSolve) {
		ufos=io.demCfgFile.unidentified_variables();
		if (ufos.size() > 0) {
			cout << " WARNING: Unrecognized DEM config file parameter(s):"; //<< INFO
			for (unsigned int i = 0; i < ufos.size(); ++i) {
				cout << " " << ufos[i];
			}
			cout << endl;
		}

	}
}

void parseDemConfigFile(IO& io,DEM& dem, GetPot& command_line) {

	if (myid==0) {
		cout<<" Parsing DEM configuration file"<<endl;
	}
	// PROBLEM NAME //////////////
	// necessary for hard coded sections of the code
	string problemNameString;
	PARSE_CLASS_MEMBER(io.demCfgFile,problemNameString,"problemName","none");
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
	PARSE_CLASS_MEMBER(io.demCfgFile,io.screenExpTime,"screenExpTime",0.0);
	ASSERT(io.screenExpTime>= 0);
	PARSE_CLASS_MEMBER(io.demCfgFile,io.partExpTime,"partExpTime",0.0);
	ASSERT(io.partExpTime>= 0);
	PARSE_CLASS_MEMBER(io.demCfgFile,io.objectExpTime,"objectExpTime",0.0);
	ASSERT(io.objectExpTime>= 0);
	PARSE_CLASS_MEMBER(io.demCfgFile,io.recycleExpTime,"recycleExpTime",0.0);
	ASSERT(io.recycleExpTime>= 0);

	// getting DEM parameters and particle initial status
	dem.discreteElementGet(io.demCfgFile,command_line);

	// problem specific parameters
	switch (problemName) {
	case NONE: {
		break;
	}
	}
}

void parseLesConfigFile(IO& io, GetPot& command_line) {

	if (myid==0) {
		cout<<" Parsing LES configuration file"<<endl;
	}

	// AGENERALE ---------------------------------------------------------
	// Iteration,restart,print
	PARSE_EXTERNAL(io.lesCfgFile,int,niter,"niter",-1);
	ASSERT(niter>=1);
	PARSE_EXTERNAL(io.lesCfgFile,int,i_rest,"i_rest",-1);
	ASSERT(i_rest>=0 && i_rest<=3);
	PARSE_EXTERNAL(io.lesCfgFile,int,print_iter_or_time,"print_iter_or_time",-1);
	ASSERT(print_iter_or_time>=0 && print_iter_or_time<=1);
	PARSE_EXTERNAL(io.lesCfgFile,int,i_print,"i_print",-1);
	ASSERT(i_print>=0);
	PARSE_EXTERNAL(io.lesCfgFile,double,i_time,"i_time",-1);
	ASSERT(i_time>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,int,i_printfile,"i_printfile",-1);
	ASSERT(i_printfile>=0 && i_printfile<=2);
	PARSE_EXTERNAL(io.lesCfgFile,bool,paraview_aver,"paraview_aver",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,paraview_inst,"paraview_inst",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,paraview_compact,"paraview_compact",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,i_newres,"i_newres",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,i_medietempo,"i_medietempo",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,i_cumulative,"i_cumulative",false);
	PARSE_EXTERNAL(io.lesCfgFile,int,iformat_newres,"iformat_newres",-1);
	ASSERT(iformat_newres>=0 && iformat_newres<=2);
	if (iformat_newres==1) {
		PARSE_CLASS_MEMBER(io.lesCfgFile,io.restartFormatString,"string_newres_format","10e18.10");
	}
	PARSE_EXTERNAL(io.lesCfgFile,bool,inf,"inf",false);

	// Characteristic numbers
	PARSE_EXTERNAL(io.lesCfgFile,double,re,"re",-1.0);
	ASSERT(re>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,double,rich,"rich",-1.0);
	ASSERT(rich>=0.0);

	// Characteristic numbers
	PARSE_EXTERNAL(io.lesCfgFile,int,ind_cou,"ind_cou",-1);
	ASSERT(ind_cou>=0);
	PARSE_EXTERNAL(io.lesCfgFile,double,cou,"cou",-1.0);
	ASSERT(cou>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,double,dt,"dt",-1.0);
	ASSERT(dt>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,int,espl,"espl",-1);
	ASSERT(espl>=0);

	// Numerical schemes,equations to solve
	PARSE_EXTERNAL(io.lesCfgFile,bool,attiva_scal,"attiva_scal",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,potenziale,"potenziale",false);
	PARSE_EXTERNAL(io.lesCfgFile,int,insc,"insc",-1);
	ASSERT(insc>=0 && insc<=2);

	// Pressure condition,pressure solver
	PARSE_EXTERNAL(io.lesCfgFile,double,omega,"om_pre",-1.0);
	ASSERT(omega>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,int,islor,"islor",-1);
	ASSERT(islor>=0);
	PARSE_EXTERNAL(io.lesCfgFile,int,nlevmultimax,"nlevmultimax",-1);
	ASSERT(nlevmultimax>=0 && nlevmultimax<=4);
	PARSE_EXTERNAL(io.lesCfgFile,double,bbx,"bbx",-1.0);
	PARSE_EXTERNAL(io.lesCfgFile,double,bby,"bby",-1.0);
	PARSE_EXTERNAL(io.lesCfgFile,double,bbz,"bbz",-1.0);
	PARSE_EXTERNAL(io.lesCfgFile,double,eps,"eps",-1.0);
	ASSERT(eps>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,int,ficycle,"ficycle",-1);
	ASSERT(ficycle>=0);
	PARSE_EXTERNAL(io.lesCfgFile,bool,freesurface,"freesurface",false);

	// Buondary condition
	PARSE_EXTERNAL(io.lesCfgFile,bool,lett,"lett",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,ibb,"ibb",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,bodyforce,"bodyforce",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,particles,"particles",false);
	if (particles) {
		ASSERT(particles && bodyforce);
	}
	PARSE_EXTERNAL(io.lesCfgFile,bool,particles,"particles",false);
	PARSE_EXTERNAL(io.lesCfgFile,int,num_iter,"num_iter",-1);
	ASSERT(num_iter>0);
	PARSE_EXTERNAL(io.lesCfgFile,int,coef_wall,"coef_wall",-1);
	ASSERT(coef_wall>=0 && coef_wall<=2);
	PARSE_EXTERNAL(io.lesCfgFile,bool,rough,"rough",false);
	PARSE_EXTERNAL(io.lesCfgFile,double,z0,"z0",-1.0);
	ASSERT(z0>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,bool,att_wm_sgs,"att_wm_sgs",false);
	if (coef_wall==0) {
		ASSERT(coef_wall==0 && !att_wm_sgs);
	}
	if (potenziale) {
		ASSERT(coef_wall==0 && potenziale);
		ASSERT(!att_wm_sgs && potenziale);
	}

	// Turbulence model
	PARSE_EXTERNAL(io.lesCfgFile,int,nsgs,"nsgs",-1);
	ASSERT(nsgs>=0);
	PARSE_EXTERNAL(io.lesCfgFile,bool,inmod,"inmod",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,inmodrho,"inmodrho",false);
	if (nsgs==0) {
		ASSERT(nsgs==0 && !inmodrho);
		ASSERT(nsgs==0 && !inmod);
	}
	PARSE_EXTERNAL(io.lesCfgFile,bool,isotropo,"isotropo",false);
	if (!isotropo) {
		ASSERT(!isotropo && nsgs!=0);
	}
	PARSE_EXTERNAL(io.lesCfgFile,double,cost,"cost",-1.0);
	ASSERT(cost>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,double,costh,"costH",-1.0);
	ASSERT(costh>=0.0);
	PARSE_EXTERNAL(io.lesCfgFile,double,costv,"costV",-1.0);
	ASSERT(costv>=0.0);

	// Forcing
	PARSE_EXTERNAL(io.lesCfgFile,int,windyes,"windyes",-1);
	ASSERT(windyes>=0);
	PARSE_EXTERNAL(io.lesCfgFile,double,latitude,"LATITUDE",-1.0);
	ASSERT(latitude>=0.0);

	// ABOUNDARY ---------------------------------------------------------
	// Periodicity
	// a more reasonable system than the old in LES_COAST
	bool csi_periodic,zita_periodic;
	extern int ip,kp;

	PARSE_CLASS_MEMBER(io.lesCfgFile,csi_periodic,"csi_periodic",false);
	PARSE_CLASS_MEMBER(io.lesCfgFile,zita_periodic,"zita_periodic",false);

	// and the sad old confusing system
	ip=(csi_periodic)?0:1;
	kp=(zita_periodic)?0:1;

	// For open boundary,if periodicity on,these are not used!
	PARSE_EXTERNAL(io.lesCfgFile,int,infout1,"infout1",-1);
	ASSERT(infout1>=0 && infout1<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,infout2,"infout2",-1);
	ASSERT(infout2>=0 && infout2<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,infout3,"infout3",-1);
	ASSERT(infout3>=0 && infout3<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,infout4,"infout4",-1);
	ASSERT(infout4>=0 && infout4<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,infout5,"infout5",-1);
	ASSERT(infout5>=0 && infout5<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,infout6,"infout6",-1);
	ASSERT(infout6>=0 && infout6<=2);

	// Settings for solid wall,only if infout=2
	PARSE_EXTERNAL(io.lesCfgFile,int,iboun1,"iboun1",-1);
	ASSERT(iboun1>=0 && iboun1<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,iboun2,"iboun2",-1);
	ASSERT(iboun2>=0 && iboun2<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,iboun3,"iboun3",-1);
	ASSERT(iboun3>=0 && iboun3<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,iboun4,"iboun4",-1);
	ASSERT(iboun4>=0 && iboun4<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,iboun5,"iboun5",-1);
	ASSERT(iboun5>=0 && iboun5<=2);
	PARSE_EXTERNAL(io.lesCfgFile,int,iboun6,"iboun6",-1);
	ASSERT(iboun6>=0 && iboun6<=2);

	// wall function at solid wall,only wfp3 and wfp4 are currently used
	PARSE_EXTERNAL(io.lesCfgFile,bool,wfp1,"wfp1",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,wfp2,"wfp2",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,wfp3,"wfp3",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,wfp4,"wfp4",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,wfp5,"wfp5",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,wfp6,"wfp6",false);
	// if any of these is active, we must check if coef_wall is on, that infout is set on wall,
	// that iboun is set on no-slip, and that periodicity is off
	if (wfp1) {
		ASSERT(wfp1 && coef_wall!=0);
		ASSERT(wfp1 && iboun1==0);
		ASSERT(wfp1 && !csi_periodic);
		ASSERT(wfp1 && infout1==2);
	}
	if (wfp2) {
		ASSERT(wfp2 && coef_wall!=0);
		ASSERT(wfp2 && iboun2==0);
		ASSERT(wfp2 && !csi_periodic);
		ASSERT(wfp2 && infout2==2);
	}
	if (wfp3) {
		ASSERT(wfp3 && coef_wall!=0);
		ASSERT(wfp3 && iboun3==0);
		ASSERT(wfp3 && infout3==2);
	}
	if (wfp4) {
		ASSERT(wfp4 && coef_wall!=0);
		ASSERT(wfp4 && iboun4==0);
		ASSERT(wfp4 && infout4==2);
	}
	if (wfp5) {
		ASSERT(wfp5 && coef_wall!=0);
		ASSERT(wfp5 && iboun5==0);
		ASSERT(wfp5 && !zita_periodic);
		ASSERT(wfp5 && infout5==2);
	}
	if (wfp6) {
		ASSERT(wfp6 && coef_wall!=0);
		ASSERT(wfp6 && iboun6==0);
		ASSERT(wfp6 && !zita_periodic);
		ASSERT(wfp6 && infout6==2);
	}

	// bodyforce in buffer region to enhance turbulence,if ibb=1 in Agenerale.in
	PARSE_EXTERNAL(io.lesCfgFile,bool,ibodybuffer1,"ibodybuffer1",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,ibodybuffer2,"ibodybuffer2",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,ibodybuffer3,"ibodybuffer3",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,ibodybuffer4,"ibodybuffer4",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,ibodybuffer5,"ibodybuffer5",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,ibodybuffer6,"ibodybuffer6",false);

	if(!ibb) {
		ASSERT(!ibb && !ibodybuffer1);
		ASSERT(!ibb && !ibodybuffer2);
		ASSERT(!ibb && !ibodybuffer3);
		ASSERT(!ibb && !ibodybuffer4);
		ASSERT(!ibb && !ibodybuffer5);
		ASSERT(!ibb && !ibodybuffer6);
	}

	if (ibodybuffer1) ASSERT(ibodybuffer1 && ibb);
	if (ibodybuffer2) ASSERT(ibodybuffer2 && ibb);
	if (ibodybuffer3) ASSERT(ibodybuffer3 && ibb);
	if (ibodybuffer4) ASSERT(ibodybuffer4 && ibb);
	if (ibodybuffer5) ASSERT(ibodybuffer5 && ibb);
	if (ibodybuffer6) ASSERT(ibodybuffer6 && ibb);


	bool isponge,ksponge;
	extern int ispon,kspon;
	PARSE_CLASS_MEMBER(io.lesCfgFile,isponge,"isponge",false);
	PARSE_CLASS_MEMBER(io.lesCfgFile,ksponge,"ksponge",false);
	ispon=(isponge)?1:0;
	kspon=(ksponge)?1:0;

	PARSE_EXTERNAL(io.lesCfgFile,int,corr_factor,"corr_factor",-1.0);
	ASSERT(corr_factor>=0);

	// ASCALARE  ---------------------------------------------------------
	// Number of scalar equations
	PARSE_EXTERNAL(io.lesCfgFile,int,nscal,"nscal",-1);
	ASSERT(nscal>=1);
	PARSE_EXTERNAL(io.lesCfgFile,bool,re_analogy,"re_analogy",false);
	PARSE_EXTERNAL_ARRAY(io.lesCfgFile,double,nscal,c_pran,"pran",-1.0);
	PARSE_EXTERNAL_ARRAY(io.lesCfgFile,double,nscal,c_prsc,"prsc",-1.0);

	// APIANISONDE  ---------------------------------------------------------
	PARSE_EXTERNAL(io.lesCfgFile,int,npiani,"npiani",-1);
	ASSERT(npiani>=0);
	PARSE_EXTERNAL_ARRAY(io.lesCfgFile,int,npiani,c_piani,"piani",-1);

	PARSE_EXTERNAL(io.lesCfgFile,int,nsonde,"nsonde",-1);
	ASSERT(nsonde>=0);
	PARSE_EXTERNAL_ARRAY(io.lesCfgFile,int,nsonde,c_sondeindexi,"indexi",-1);
	PARSE_EXTERNAL_ARRAY(io.lesCfgFile,int,nsonde,c_sondeindexj,"indexj",-1);
	PARSE_EXTERNAL_ARRAY(io.lesCfgFile,int,nsonde,c_sondeindexk,"indexk",-1);

	// AFILTRAGGIO  ---------------------------------------------------------

	PARSE_EXTERNAL(io.lesCfgFile,bool,ifiltro,"ifiltro",false);
	PARSE_EXTERNAL(io.lesCfgFile,int,nfiltro,"nfiltro",-1);

	PARSE_EXTERNAL(io.lesCfgFile,bool,filtrou,"filtrou",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,filtrov,"filtrov",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,filtrow,"filtrow",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,filtrorho,"filtrorho",false);
	PARSE_EXTERNAL(io.lesCfgFile,bool,filtrofi,"filtrofi",false);

	if (filtrou) ASSERT(filtrou && ifiltro);
	if (filtrov) ASSERT(filtrov && ifiltro);
	if (filtrow) ASSERT(filtrow && ifiltro);
	if (filtrorho) ASSERT(filtrorho && ifiltro);
	if (filtrofi) ASSERT(filtrofi && ifiltro);


	PARSE_EXTERNAL(io.lesCfgFile,int,xstart,"xstart",-1);
	PARSE_EXTERNAL(io.lesCfgFile,int,xend,"xend",-1);
	PARSE_EXTERNAL(io.lesCfgFile,int,ystart,"ystart",-1);
	PARSE_EXTERNAL(io.lesCfgFile,int,yend,"yend",-1);
	PARSE_EXTERNAL(io.lesCfgFile,int,zstart,"zstart",-1);
	PARSE_EXTERNAL(io.lesCfgFile,int,zend,"zend",-1);


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

}

int main(int argc,char* argv[]){

	extern int niter,ktime;
	extern bool termina;

	MPI::Status status;

	MPI::Init();

	myid = MPI::COMM_WORLD.Get_rank();
	nproc = MPI::COMM_WORLD.Get_size();

	// DECLARATION OF VARIABLES - Input-Output ///////////////
	IO io;
	// DECLARATION OF VARIABLES - DEM ///////////////
	DEM dem;

    // print some info for restorability
    time_t t=time(0);   // get time now
    io.now=localtime(&t);
    if (myid==0) {
    	cout<<"Binary was compiled on "<<__DATE__<<" at "<<__TIME__<<endl;
    	cout<<"Invocation on "<<(io.now->tm_year + 1900)<<"-"<<(io.now->tm_mon+1)<<"-"<<io.now->tm_mday<<
    		", at "<<io.now->tm_hour<<":"<<io.now->tm_min<<":"<<io.now->tm_sec<<endl;
    	cout<<"Invocation line: " << argv[0];
    	for (int i=1; i<argc; ++i) {
    		cout<<" "<<argv[i];
    	}
    	cout<<endl;
    }


	// parsing command line and configuration file
	GetPot command_line(argc,argv);

	parseCommandLine(io,command_line);
	if (demSolve) {
		// parsing DEM input file
		GetPot demCfgFile(io.demConfigFileName);
		io.demCfgFile=demCfgFile;
		parseDemConfigFile(io,dem,command_line);

	}

	if (lesSolve) {
		// parsing LES input file
		GetPot lesCfgFile(io.lesConfigFileName);
		io.lesCfgFile=lesCfgFile;
		parseLesConfigFile(io,command_line);

		// send strings to Fortran
		get_cpp_strings(io.workDirectory.c_str(),io.lesGridFile.c_str(),io.lesRestartFile.c_str(),io.restartFormatString.c_str());

		MPI::COMM_WORLD.Barrier();

		//  INITIALIZATION //////////////////////////////////
		// first step: read domain size from grid
		les_domain_size();
	}

	ktime=0;

	// initializing input-output
	if (demSolve) {
		io.initialize();
	}

	// first
	if (lesSolve) {
		les_initialize_part1();
	}

	if (demSolve) {
		// initializing DEM parameters
		dem.discreteElementInit();
	}

	if (myid==0) {
		printUfo(io,command_line);
	}

	if (lesSolve) {
		les_initialize_part2();
	}
	//dem.discreteElementStep1();


	// based on the solver, choose time step and
	int maxIterations=0;
	if (lesSolve) {
		maxIterations=niter;
	}
	else if (demSolve) {
		maxIterations=dem.maxTimeSteps;
		if (maxIterations==0) {
			maxIterations=std::numeric_limits<int>::max();
		}
	}
	else {
		ASSERT(demSolve || lesSolve)
	}

	ktime=0;
	termina=false;

	// CYCLE /////////////////////////////
	while(!termina) {

		++ktime;

		if (lesSolve && ktime>maxIterations) {
			termina=true;
		}

		if (!lesSolve && dem.demTime>dem.maxTime) {
			termina=true;
		}

		extern double dt;

		MPI::COMM_WORLD.Barrier();

		if (lesSolve && demSolve && myid==0) {
			cout<<" ~~~~~~~ Enter DEM-1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
			cout<<" DEM time = "<<dem.demTime<<endl;
		}

		// advance one step in time

		if (demSolve) {
			if (lesSolve) {

				dem.demTime+=dt;
				// for the moment being, same dt for LES and DEM
				dem.deltat=dt;
			}
			else {
				dem.demTime+=dem.deltat;
			}
			// predictor
			dem.discreteElementStep1();
		}

		if (lesSolve) {
			MPI::COMM_WORLD.Barrier();
			if (lesSolve && demSolve && myid==0) {
				cout<<" ~~~~~~~ Enter LES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
			}

			les_core();
		}

		if (lesSolve && demSolve && myid==0) {
			cout<<" ~~~~~~~ Enter DEM-2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		}

		if (demSolve) {
			// corrector
			dem.discreteElementStep2();

			if (myid==0) {
				io.outputStep(dem.demTime,dem);
				++io.currentTimeStep;
			}
		}
	}

	if (lesSolve) {
		les_finalize();
	}

	MPI::Finalize();

	return 0;

}
