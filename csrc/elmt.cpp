
#include "elmt.h"

#include "myvector.h"

using namespace std;

void elmt::elmtShow() const {
	cout<<"Element number "<<index<<" with "<<size<<" particles of radius "<<radius<<"\n";
	cout<<"Position: "; x0.show();
	cout<<"; Velocity: "; x1.show(); cout<<";\n";
	cout<<"Orientation: "; q0.show(); cout<<";\n";
	cout<<"Inertia: "; I.show(); cout<<";\n";
	cout<<"Mass: "<<m<<";\n";
}

void elmt::initialize(const double& partDensity, std::vector <vecList>& prototypes, tVect& demF) {

	// translational degree of freedom
	xp0=x0;
	xp1=x1;
	xp2.reset();
	x2.reset();
	xp3.reset();
	x3.reset();
	xp4.reset();
	x4.reset();
	xp5.reset();
	x5.reset();
	x0history.reset();
	// rotational degree of freedom (for the moment initial status is supposed to be 0)
	qp0=q0=tQuat(1.0,0.0,0.0,0.0);
	wGlobal.reset();
	wLocal.reset();
	wpGlobal.reset();
	wpLocal.reset();
	tQuat wQuat=vec2quat(wGlobal);
	qp1=q1=0.5*q0.multiply(wQuat);
	qp2=q2=tQuat(0.0,0.0,0.0,0.0);
	qp3=q3=tQuat(0.0,0.0,0.0,0.0);
	qp4=q4=tQuat(0.0,0.0,0.0,0.0);
	qp5=q5=tQuat(0.0,0.0,0.0,0.0);

	// calculated variables (the element is supposed to be a sphere for the moment)
	// mass
	const double singleMass=4.0/3.0*partDensity*M_PI*radius*radius*radius;
	m=size*singleMass;
	// inertia moment (diagonal) - Huygens-Steiner theorem
	// inertia of single spheres
	I=size*2.0/5.0*singleMass*radius*radius*tVect(1.0,1.0,1.0);
	// transport components
	for (unsigned int n=0; n<size; ++n) {
		I+=singleMass*radius*radius*prototypes[size][n].transport();
	}

	// initialize forces
	FParticle.reset();
	FWall.reset();
	FGrav=demF*m;
	MParticle.reset();
	MWall.reset();

}

void elmt::resetVelocity() {
	// resets translational and rotational velocity to zero

	// translational
	xp1.reset();
	x1.reset();
	xp2.reset();
	x2.reset();
	xp3.reset();
	x3.reset();
	xp4.reset();
	x4.reset();
	xp5.reset();
	x5.reset();
	// rotational
	wGlobal.reset();
	wLocal.reset();
	wpGlobal.reset();
	wpLocal.reset();
	tQuat wQuat=vec2quat(wGlobal);
	qp1=q1=0.5*q0.multiply(wQuat);
	qp2=q2=tQuat(0.0,0.0,0.0,0.0);
	qp3=q3=tQuat(0.0,0.0,0.0,0.0);
	qp4=q4=tQuat(0.0,0.0,0.0,0.0);
	qp5=q5=tQuat(0.0,0.0,0.0,0.0);
}

void elmt::generateParticles(unsigned int& globalIndex, particleList& particles,  const std::vector<vecList>& prototypes) {
	components.resize(size);

	for (int i=0; i<size; ++i) {
		particle dummyPart;
		dummyPart.particleIndex=globalIndex;
		components[i]=dummyPart.particleIndex;
		dummyPart.clusterIndex=index;
		dummyPart.r=radius;
		dummyPart.D=dummyPart.r*2.0;
		dummyPart.protoIndex=i;
		dummyPart.updatePredicted(*this, prototypes);
		particles.push_back(dummyPart);
		++globalIndex;
	}
}

void elmt::predict(const double c[]) {

	xp0=x0+x1*c[0]+x2*c[1]+x3*c[2]+x4*c[3]+x5*c[4];
	xp1=x1+x2*c[0]+x3*c[1]+x4*c[2]+x5*c[3];
	xp2=x2+x3*c[0]+x4*c[1]+x5*c[2];
	xp3=x3+x4*c[0]+x5*c[1];
	xp4=x4+x5*c[0];
	xp5=x5;

	qp0=q0+q1*c[0]+q2*c[1]+q3*c[2]+q4*c[3]+q5*c[4];
	qp1=q1+q2*c[0]+q3*c[1]+q4*c[2]+q5*c[3];
	qp2=q2+q3*c[0]+q4*c[1]+q5*c[2];
	qp3=q3+q4*c[0]+q5*c[1];
	qp4=q4+q5*c[0];
	qp5=q5;

	const tQuat qp0adj=qp0.adjoint();
	wpGlobal=2.0*quat2vec( qp1.multiply(qp0adj) );
	wpLocal=2.0*quat2vec( qp0adj.multiply(qp1) );

}

void elmt::correct(const double coeff[]) {

	const tVect x2Corr=x2-xp2;

	const tVect newPos=xp0+x2Corr*coeff[0];
	x0history+=newPos-x0;
	x0=newPos;
	x1=xp1+x2Corr*coeff[1];
	// x2 calculated directly at the end of force routine
	x3=xp3+x2Corr*coeff[3];
	x4=xp4+x2Corr*coeff[4];
	x5=xp5+x2Corr*coeff[5];

	xp0=x0;
	xp1=x1;
	xp2=x2;
	xp3=x3;
	xp4=x4;
	xp5=x5;

	tQuat q2Corr=q2-qp2;

	q0=qp0+q2Corr*coeff[0];
	q1=qp1+q2Corr*coeff[1];
	// q2 calculated directly at the end of force routine
	q3=qp3+q2Corr*coeff[3];
	q4=qp4+q2Corr*coeff[4];
	q5=qp5+q2Corr*coeff[5];
	//normalization of q0
	q0.normalize();
	const tQuat q0adj=q0.adjoint();
	wGlobal=2.0*quat2vec( q1.multiply(q0adj) );
	wLocal=2.0*quat2vec( q0adj.multiply(q1) );

	qp0=q0;
	qp1=q1;
	qp2=q2;
	qp3=q3;
	qp4=q4;
	qp5=q5;

}

void elmt::translate(const tVect& transVec) {
	x0+=transVec;
	xp0+=transVec;
}

void particle::updatePredicted(const elmt& motherElmt, const std::vector <vecList>& prototypes) {
	// updating position
//	if (myid==0) {
//		cout<<"x0="<<x0.dot(X)<<endl;
//		cout<<"<motherElmt.index="<<motherElmt.index<<endl;
//		cout<<" motherElmt.xp0="<<motherElmt.xp0.dot(X)<<endl;
//		cout<<" r="<<r<<endl;
//		cout<<" project(prototypes[motherElmt.size][protoIndex],motherElmt.qp0)="<<project(prototypes[motherElmt.size][protoIndex],motherElmt.qp0).dot(X)<<endl;
//	}
	x0=motherElmt.xp0+r*project(prototypes[motherElmt.size][protoIndex],motherElmt.qp0);
	// updating radius (distance of particle center of mass to element center of mass)
	radiusVec=x0-motherElmt.xp0;
	// updating particle speed
	x1=motherElmt.xp1+motherElmt.wpGlobal.cross(radiusVec);
}

void particle::updateCorrected(const elmt& motherElmt, const std::vector <vecList>& prototypes) {
	// updating position
	x0=motherElmt.x0+r*project(prototypes[motherElmt.size][protoIndex],motherElmt.q0);
	// updating radius (distance of particle center of mass to element center of mass)
	radiusVec=x0-motherElmt.x0;
	// updating particle speed
	x1=motherElmt.x1+motherElmt.wGlobal.cross(radiusVec);
}

void particle::ghostUpdate(particle& originParticle, tVect& pbcVector) {
	// updating position
	x0=originParticle.x0+pbcVector;
	// updating radius (distance of particle center of mass to element center of mass)
	radiusVec=originParticle.radiusVec;
	// updating particle speed
	x1=originParticle.x1;
}
