
#include "IO.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

void IO::initialize(){

	currentTimeStep = 0;

	partDirectory=workDirectory+"/particleData";

	if (myid==0) {
		int a;
		a=mkdir(workDirectory.c_str(),0777);
		cout<<"Work directory created = "<<workDirectory<<". Result: "<<a<<"\n";
		a=mkdir(partDirectory.c_str(),0777);
		cout<<"Work directory created = "<<partDirectory<<". Result: "<<a<<"\n";
	}
	// build data file names for shell
	partFileFormat=partDirectory+"/part%010u.vtu";
	recycleFileFormat=partDirectory+"/recycle%010u.vtu";
	objectFileFormat=partDirectory+"/object%010u.vtu";

	//  initializing output file
	exportFileName=workDirectory+"/export.dat";

	//  initializing max speed file
	maxSpeedFileName=workDirectory+"/maxVel.dat";
	maxSpeedFile.open(maxSpeedFileName.c_str(), ios::app);
	maxSpeedFile<<"time maxFluidSpeed\n";
	maxSpeedFile.close();

	//  initializing force file
	forceFileName=workDirectory+"/force.dat";
	forceFile.open(forceFileName.c_str(), ios::app);
	forceFile<<"time collisionForces hydraulicForces\n";
	forceFile.close();

	//  initializing diffusion file
	diffFileName=workDirectory+"/diffusion.dat";
	diffFile.open(diffFileName.c_str(), ios::app);
	diffFile<<"time msd msdX msdY msdZ\n";
	diffFile.close();

	lastScreenExp=0;
	lastPartExp=0;
	lastObjectExp=0;
	lastOutputExp=0;
	lastRecycleExp=0;

}

void IO::outputStep(const double& time, const DEM& dem){

	realTime=time;

	// FILE CREATION PHASE  ////////////////////////////////////////////////////////////////
	createParaviewFiles(dem);

	// PLOTTING PHASE  ////////////////////////////////////////////////////////////////
	const unsigned int screenExpCounter = (screenExpTime>0 ? static_cast<unsigned int>(realTime/screenExpTime) : 0);

	if (screenExpCounter>lastScreenExp) {

		lastScreenExp = screenExpCounter;

		exportFile.open(exportFileName.c_str(), ios::app);

		// current iteration and real time
		cout<<currentTimeStep<<"; time="<<setprecision(5)<<realTime<<"\t";
		exportFile<<currentTimeStep<<"; time="<<setprecision(5)<<realTime<<"\t";

		if (dem.elmts.size()) {
			exportMaxSpeedParticles(dem);
			exportForces(dem);
			exportDiffusion(dem);
		}

		// closing file
		cout<<"\n";
		exportFile<<"\n";
		exportFile.close();
		cout.flush();
	}
}

void IO::outputFinal(){
	// drag file closure
	exportFile.close();
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// paraview files

void IO::createParaviewFiles(const DEM& dem) {
	// write vtk at regular interval defined with the input file


	const unsigned int partExpCounter=(partExpTime>0 ? static_cast<unsigned int>(realTime/partExpTime)+1 : 0);
	if (partExpCounter>lastPartExp) {

		lastPartExp=partExpCounter;
		sprintf(filePathBuffer, partFileFormat.c_str(), currentTimeStep);
		exportParaviewParticles(dem.elmts, dem.particles, filePathBuffer);
	}
	const unsigned int recycleExpCounter=(recycleExpTime>0 ? static_cast<unsigned int>(realTime/recycleExpTime)+1 : 0);
	if (recycleExpCounter>lastRecycleExp) {
		lastRecycleExp=recycleExpCounter;
		sprintf(filePathBuffer, partFileFormat.c_str(), currentTimeStep);
		exportRecycleParticles(dem.elmts);
	}
	const unsigned int objectExpCounter=(objectExpTime>0 ? static_cast<unsigned int>(realTime/objectExpTime)+1 : 0);
	if (objectExpCounter>lastObjectExp) {
		lastObjectExp=objectExpCounter;
		sprintf(filePathBuffer, objectFileFormat.c_str(), currentTimeStep);
		exportParaviewObjects(dem.objects, filePathBuffer);
	}
}

void IO::exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile){

	const int one=1;
	const int Pnumber=particles.size();
	const tVect nx(1.0,0.0,0.0),ny(0.0,1.0,0.0),nz(0.0,0.0,1.0);
	tVect n1,n2,n3;

	// file opening
	ofstream paraviewParticleFile;
	paraviewParticleFile.open(particleFile.c_str());
	// writing on header file
	paraviewParticleFile<<"<?xml version=\"1.0\"?>\n";
	paraviewParticleFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
	paraviewParticleFile<<" <UnstructuredGrid GhostLevel=\"0\">\n";
	paraviewParticleFile<<"  <Piece NumberOfPoints=\""<<Pnumber<<"\" NumberOfCells=\""<<Pnumber<<"\">\n";
	paraviewParticleFile<<"   <PointData>\n";
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"radius\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		paraviewParticleFile<<particles[i].r<<"\n";
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"mass\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		paraviewParticleFile<<elmts[particles[i].clusterIndex].m<<"\n";
	}
	//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"particleIndex\"/>\n";
	//    for (int i=0; i<Pnumber; ++i){
	//        paraviewParticleFile<<particles[i].particleIndex<<"\n";
	//    }
	//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"clusterIndex\"/>\n";
	//    for (int i=0; i<Pnumber; ++i){
	//        paraviewParticleFile<<particles[i].clusterIndex<<"\n";
	//    }
	//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"tableCell\"/>\n";
	//    for (int i=0; i<Pnumber; ++i){
	//        paraviewParticleFile<<particles[i].tableCell<<"\n";
	//    }
	//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"Inertia\" NumberOfComponents=\"3\"/>\n";
	//    for (int i=0; i<Pnumber; ++i){
	//        elmts[particles[i].clusterIndex].I.printLine(paraviewParticleFile);
	//    }
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		particles[i].x1.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"w\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		elmts[particles[i].clusterIndex].wGlobal.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FParticle\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		elmts[particles[i].clusterIndex].FParticle.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FGrav\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		elmts[particles[i].clusterIndex].FGrav.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"MParticle\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		elmts[particles[i].clusterIndex].MParticle.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FWall\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		elmts[particles[i].clusterIndex].FWall.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FHydroPressure\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		elmts[particles[i].clusterIndex].FHydroPressure.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FHydroShear\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		elmts[particles[i].clusterIndex].FHydroShear.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"MWall\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		elmts[particles[i].clusterIndex].MWall.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"n1\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		n1=project(nx,elmts[particles[i].clusterIndex].q0);
		n1.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"n2\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		n2=project(ny,elmts[particles[i].clusterIndex].q0);
		n2.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"n3\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		n3=project(nz,elmts[particles[i].clusterIndex].q0);
		n3.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"   </PointData>\n";
	paraviewParticleFile<<"   <CellData>\n";
	paraviewParticleFile<<"   </CellData>\n";
	paraviewParticleFile<<"   <Points>\n";
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Pnumber; ++i){
		particles[i].x0.printLine(paraviewParticleFile);
	}
	paraviewParticleFile<<"   </Points>\n";
	paraviewParticleFile<<"   <Cells>\n";
	paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
	for (int i=1; i<Pnumber+1; ++i){
		paraviewParticleFile<<i-1<<"\n";
	}
	paraviewParticleFile<<"    </DataArray>\n";
	paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
	for (int i=1; i<Pnumber+1; ++i){
		paraviewParticleFile<<i<<"\n";
	}
	paraviewParticleFile<<"    </DataArray>\n";
	paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
	for (int i=0; i<Pnumber; ++i){
		paraviewParticleFile<<one<<"\n";
	}
	paraviewParticleFile<<"    </DataArray>\n";
	paraviewParticleFile<<"   </Cells>\n";
	paraviewParticleFile<<"  </Piece>\n";
	//    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
	paraviewParticleFile<<" </UnstructuredGrid>\n";
	paraviewParticleFile<<"</VTKFile>";
	// header file closing
	paraviewParticleFile.close();
}

void IO::exportParaviewObjects(const objectList& objects, const string& objectFile){

	const int one=1;
	const int Onumber=objects.size();

	// file opening
	ofstream paraviewObjectFile;
	paraviewObjectFile.open(objectFile.c_str());
	// writing on header file
	paraviewObjectFile<<"<?xml version=\"1.0\"?>\n";
	paraviewObjectFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
	paraviewObjectFile<<" <UnstructuredGrid GhostLevel=\"0\">\n";
	paraviewObjectFile<<"  <Piece NumberOfPoints=\""<<Onumber<<"\" NumberOfCells=\""<<Onumber<<"\">\n";
	paraviewObjectFile<<"   <PointData>\n";
	paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"radius\"/>\n";
	for (int i=0; i<Onumber; ++i){
		paraviewObjectFile<<objects[i].radius<<"\n";
	}
	paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"objectIndex\"/>\n";
	for (int i=0; i<Onumber; ++i){
		paraviewObjectFile<<objects[i].index<<"\n";
	}
	paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Onumber; ++i){
		objects[i].x1.printLine(paraviewObjectFile);
	}
	paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"FParticle\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Onumber; ++i){
		objects[i].FParticle.printLine(paraviewObjectFile);
	}
	paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"FHydroPressure\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Onumber; ++i){
		objects[i].FHydroPressure.printLine(paraviewObjectFile);
	}
	paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"FHydroShear\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Onumber; ++i){
		objects[i].FHydroShear.printLine(paraviewObjectFile);
	}
	paraviewObjectFile<<"   </PointData>\n";
	paraviewObjectFile<<"   <CellData>\n";
	paraviewObjectFile<<"   </CellData>\n";
	paraviewObjectFile<<"   <Points>\n";
	paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
	for (int i=0; i<Onumber; ++i){
		objects[i].x0.printLine(paraviewObjectFile);
	}
	paraviewObjectFile<<"   </Points>\n";
	paraviewObjectFile<<"   <Cells>\n";
	paraviewObjectFile<<"    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
	for (int i=1; i<Onumber+1; ++i){
		paraviewObjectFile<<i-1<<"\n";
	}
	paraviewObjectFile<<"    </DataArray>\n";
	paraviewObjectFile<<"    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
	for (int i=1; i<Onumber+1; ++i){
		paraviewObjectFile<<i<<"\n";
	}
	paraviewObjectFile<<"    </DataArray>\n";
	paraviewObjectFile<<"    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
	for (int i=0; i<Onumber; ++i){
		paraviewObjectFile<<one<<"\n";
	}
	paraviewObjectFile<<"    </DataArray>\n";
	paraviewObjectFile<<"   </Cells>\n";
	paraviewObjectFile<<"  </Piece>\n";
	//    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
	paraviewObjectFile<<" </UnstructuredGrid>\n";
	paraviewObjectFile<<"</VTKFile>";
	// header file closing
	paraviewObjectFile.close();
}

void IO::exportRecycleParticles(const elmtList& elmts){
	// exports a format readable by this code itself, allowing the re-use of computed particle data.
	// needs further work...

	int i;
	const char* charParticleFile;
	ofstream recycleParticleFile;
	string title, particleFile;
	title=partDirectory+"/recycle_";
	stringstream ss;
	ss<<currentTimeStep;
	particleFile=title+ss.str()+".dat";
	charParticleFile=particleFile.c_str();

	double zero=0.0;
	// header file opening
	recycleParticleFile.open(charParticleFile);

	// total number of element

	recycleParticleFile<<elmts.size()<<"\n";

	// ATTENTION this does not work for clusters unless also rotational velocity and position are included in the initial file :(
	for (int i=0; i<elmts.size(); ++i) {
		// import variables
		recycleParticleFile<<elmts[i].index<<" ";
		recycleParticleFile<<elmts[i].size<<" ";
		recycleParticleFile<<elmts[i].radius<<" ";
		elmts[i].x0.print(recycleParticleFile);
		recycleParticleFile<<zero<<" ";
		recycleParticleFile<<zero<<" ";
		recycleParticleFile<<zero<<"\n";
	}
	recycleParticleFile.close();
}

void IO::exportParaviewBox(const double& time){
	int i;
	const char* charParticleFile;
	ofstream paraviewParticleFile;
	string title, particleFile;
	title=partDirectory+"/box_";
	stringstream ss;
	ss<<currentTimeStep;
	particleFile=title+ss.str()+".vtu";
	charParticleFile=particleFile.c_str();

	int one=1;
	int Pnumber=1;
	double boxVel=0.15/1.5;
	// start printing all the crap required for Paraview
	// header file opening
	paraviewParticleFile.open(charParticleFile);
	// writing on header file
	paraviewParticleFile<<"<?xml version=\"1.0\"?>\n";
	paraviewParticleFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
	paraviewParticleFile<<" <UnstructuredGrid GhostLevel=\"0\">\n";
	paraviewParticleFile<<"  <Piece NumberOfPoints=\""<<Pnumber<<"\" NumberOfCells=\""<<Pnumber<<"\">\n";
	paraviewParticleFile<<"   <PointData>\n";
	paraviewParticleFile<<"   </PointData>\n";
	paraviewParticleFile<<"   <CellData>\n";
	paraviewParticleFile<<"   </CellData>\n";
	paraviewParticleFile<<"   <Points>\n";
	paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
	for (i=0; i<Pnumber; ++i){
		paraviewParticleFile<<0.17<<" "<<0.25<<" "<<0.08+boxVel*time<<"\n"; //third component to be added in the future
	}
	paraviewParticleFile<<"   </Points>\n";
	paraviewParticleFile<<"   <Cells>\n";
	paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
	for (i=1; i<Pnumber+1; ++i){
		paraviewParticleFile<<i-1<<"\n";
	}
	paraviewParticleFile<<"    </DataArray>\n";
	paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
	for (i=1; i<Pnumber+1; ++i){
		paraviewParticleFile<<i<<"\n";
	}
	paraviewParticleFile<<"    </DataArray>\n";
	paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
	for (i=0; i<Pnumber; ++i){
		paraviewParticleFile<<one<<"\n";
	}
	paraviewParticleFile<<"    </DataArray>\n";
	paraviewParticleFile<<"   </Cells>\n";
	paraviewParticleFile<<"  </Piece>\n";
	//    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
	paraviewParticleFile<<" </UnstructuredGrid>\n";
	paraviewParticleFile<<"</VTKFile>";
	// header file closing
	paraviewParticleFile.close();
}

// print export stuff

void IO::exportMaxSpeedParticles(const DEM& dem) {
	// calculating speed for real time check
	// particle max velocity
	double maxPartSpeed2=0.0;
	unsigned int velIndex=0;
	for (int n=0;n<dem.elmts.size();++n) {
		const double thisVel=dem.elmts[n].x1.norm2();
		if (maxPartSpeed2<thisVel) {
			maxPartSpeed2=thisVel;
			velIndex=dem.elmts[n].index;
		}
	}
	const double maxPartSpeed=sqrt(maxPartSpeed2);
	cout<<"MaxPSpeed= "<<maxPartSpeed<<"("<<velIndex<<")\t";
	exportFile<<"MaxPSpeed= "<<maxPartSpeed<<"("<<velIndex<<")\t";
	//    cout<<"Force= ";
	//    dem.elmts[0].FHydro.show();
	//    cout<<"\t";
}

void IO::exportDiffusion(const DEM& dem) {
	double msd=0.0;
	tVect msdComp;

	diffFile.open(diffFileName.c_str(), ios::app);
	diffusionCoefficient(dem.elmts, msd, msdComp);
	diffFile<<realTime<<" "<<msd<<" "<<msdComp.dot(X)<<" "<<msdComp.dot(Y)<<" "<<msdComp.dot(Z)<<"\n";
	diffFile.close();
	cout<<"Diff = "<<msdComp.dot(Z)/realTime<<"\t";
	exportFile<<"Diff = "<<msdComp.dot(Z)/realTime<<"\t";

	if (currentTimeStep>samplingStart) {
		msdSum+=msd;
		msdSum2+=msd*msd;
		msdVecSum+=msdComp;
		msdVecSum2+=tVect( msdComp.dot2(X), msdComp.dot2(Y), msdComp.dot2(Z) );
		++msdSumN;
	}



}

void IO::exportForces(const DEM& dem) {
	// printing force
	// hydraulic and collision forces
	forceFile.open(forceFileName.c_str(), ios::app);
	double collisionForce=collisionForceTot(dem.elmts);
	forceFile<<realTime<<" "<<collisionForce<<"\n";
	forceFile.close();
}

void IO::exportForceObstacle(const DEM& dem) {
	// printing total force on obstacle
	obstacleFile.open(obstacleFileName.c_str(), ios::app);
	tVect totObstacleForce=totForceObject(dem.objects);
	obstacleFile<<realTime<<" ";
	totObstacleForce.printLine(obstacleFile);
	obstacleFile.close();
}

// data elaboration

double IO::totParticleMass(const elmtList& elmts) const {
	// prints the total mass in the free fluid domain
	double mass=0.0;

	for (int n=0; n<elmts.size(); ++n){
		// calculate mass
		mass+=elmts[n].m;
	}

	return mass;
}

void  IO::diffusionCoefficient(const elmtList& elmts, double& msd, tVect& msdComp) const {
	// prints the mean square displacement both in module and components in the cartesian directions
	msd=0.0;
	msdComp.reset();
	tVect xDir=tVect(1.0,0.0,0.0);
	tVect yDir=tVect(0.0,1.0,0.0);
	tVect zDir=tVect(0.0,0.0,1.0);

	for (int n=0; n<elmts.size(); ++n) {
		msd+=elmts[n].x0history.norm2();
		msdComp+=tVect(elmts[n].x0history.dot2(xDir), elmts[n].x0history.dot2(yDir), elmts[n].x0history.dot2(zDir));
	}

	msd/=elmts.size();
	msdComp/=elmts.size();

}

tVect IO::centerOfMassParticles(const elmtList& elmts) const {
	// prints the total mass in the free fluid domain
	tVect center(0.0,0.0,0.0);
	double totMass(0.0);

	for (int p=0; p<elmts.size(); ++p){
		center+=elmts[p].m*elmts[p].x0;
		totMass+=elmts[p].m;
	}

	return center/totMass;
}

double IO::collisionForceTot(const elmtList& elmts) const {
	double totForce=0.0;

	for (int n=0; n<elmts.size(); ++n){
		totForce+=elmts[n].FParticle.norm();
		totForce+=elmts[n].FWall.norm();
	}

	return totForce;
}

tVect IO::totForceObject(const objectList& objects) const {
	// prints the total mass in the free fluid domain
	tVect totForce(0.0,0.0,0.0);

	for (int o=0; o<objects.size(); ++o){
		totForce+=objects[o].FParticle+objects[o].FHydroPressure+objects[o].FHydroShear;
	}

	return totForce;
}


