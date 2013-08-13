///////////////////////////////////////////////////
//
//   Demo code about  
// 
//     - creating a physical system 
//     - add/remove rigid bodies
//     - create mechanical joints between bodies
//	   - perform a simulation
//   
//	 CHRONO 
//   ------
//   Multibody dinamics engine
//  
// ------------------------------------------------ 
// 	 Copyright:Alessandro Tasora / DeltaKnowledge
//             www.deltaknowledge.com
// ------------------------------------------------ 
///////////////////////////////////////////////////


#include "physics/CHapidll.h" 
#include "physics/CHsystem.h"
#include "core/CHrealtimeStep.h"
#include "omp.h"

// Use the namespace of Chrono

using namespace chrono;

double scale(50); 
double friction(0.4f);

int numParticles(15000);
double particleRad(300*scale); // micrometres 63-180 must be less than 300 to see patterns
double particleDensity(0.005); // g/cm^3
 
double cylRad(0.10*scale); // m 0.10
double cylLen(0.0010*scale); // m 0.14
double wallThickness(0.005*scale); // m (0.005)
int numStaves(128);
double rpm(8); // 8 rpm

double timestep(0.00005);	//_____\ (5e-5,75) Still not enough, and already very slow 70m...190min...400min per frame (ones,teens,twenties frame#)
int iterations(100);		//     / trying (5e-5,100) and using omp_set_num_threads (I think APGD is not parallelized, but maybe some other parts of C::E are)

double simDuration(7.5);
int everyFrame(25);


void createBodies(ChSystem & sys)
{
	// Create rigid bodies
	GetLog() << "Creating bodies...";

	// truss
	ChSharedBodyPtr truss(new ChBody);   
	sys.AddBody(truss);
	truss->SetBodyFixed(true);

	ChSharedBodyPtr backWall(new ChBody);
	backWall->SetPos(ChVector<double>(0,0,cylLen/2.0+wallThickness/2.0));
	backWall->SetBodyFixed(true);
	backWall->GetCollisionModel()->AddBox(cylRad*2.5/2.0,cylRad*2.5/2.0, wallThickness/2.0); // Boxes are defined as half dimensions
	backWall->SetCollide(true);

	ChSharedBodyPtr frontWall(new ChBody);
	frontWall->SetPos(ChVector<double>(0,0,-cylLen/2.0-wallThickness/2.0));
	frontWall->SetBodyFixed(true);
	frontWall->GetCollisionModel()->AddBox(cylRad*2.5/2.0,cylRad*2.5/2.0, wallThickness/2.0);
	frontWall->SetCollide(true);

	frontWall->GetCollisionModel()->SetFamily(1);
	frontWall->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	backWall->GetCollisionModel()->SetFamily(1);
	backWall->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	backWall ->SetFriction(0.0f);
	frontWall->SetFriction(0.0f);

	sys.AddBody(backWall);
	sys.AddBody(frontWall);


	// Particles
	double particleMass = particleDensity/1000000 * 4/3*CH_C_PI*pow(particleRad/1000000,3);
	double particleInertia = particleRad/1000000*particleRad/1000000*particleMass*0.4;

	for (int np = 0; np <numParticles/2; np++) 
	{
		ChSharedBodyPtr greenparticle;

		greenparticle=ChSharedBodyPtr(new ChBody);
		greenparticle->SetPos(ChVector<>((CHrandom()-1)*1.3/2.0*cylRad-0.004, (CHrandom()-0.9)*0.7*cylRad, (CHrandom()-0.5)*cylLen));
		greenparticle->SetMass(particleMass); // density * volume
		greenparticle->SetInertiaXX(ChVector<>(particleInertia,particleInertia,particleInertia));
		greenparticle->GetMaterialSurface()->SetFriction(friction);
		greenparticle->GetMaterialSurface()->SetRollingFriction(0.0);
		greenparticle->GetCollisionModel()->AddSphere(particleRad/1000000);
		greenparticle->SetCollide(true);
		sys.AddBody(greenparticle);
	}


	for (int np = 0; np <numParticles/2; np++) 
	{
		ChSharedBodyPtr redparticle;

		redparticle=ChSharedBodyPtr(new ChBody);
		redparticle->SetPos(ChVector<double>((CHrandom())*1.3/2.0*cylRad+0.004, (CHrandom()-0.9)*0.7*cylRad, (CHrandom()-0.5)*cylLen));
		redparticle->SetMass(particleMass); // density * volume
		redparticle->SetInertiaXX(ChVector<>(particleInertia,particleInertia,particleInertia));
		redparticle->GetMaterialSurface()->SetFriction(friction);
		redparticle->GetMaterialSurface()->SetRollingFriction(0.0);
		redparticle->GetCollisionModel()->AddSphere(particleRad/1000000);
		redparticle->SetCollide(true);
		sys.AddBody(redparticle);
	}

	GetLog() << "DONE\n";
}


int main(int argc, char* argv[])
{
	// Set max number of threads
	omp_set_num_threads(64);

	// The DLL_CreateGlobals() - DLL_DeleteGlobals(); pair is needed if
	// global functions are needed.
	DLL_CreateGlobals();

	// The physical system: it contains all physical objects.
	ChSystem mphysicalSystem;

	// Adjust gravity to scale
	mphysicalSystem.Set_G_acc(ChVector<double>(0,-9.81*scale,0));


	// Create rigid bodies.
	createBodies(mphysicalSystem);


	// barrel
	ChSharedBodyPtr barrel(new ChBody);
	mphysicalSystem.AddBody(barrel);


	for (int i = 0; i < numStaves; i++)
	{
		double ang = 2*CH_C_PI/numStaves * i;

		barrel->GetCollisionModel()->AddBox(
			sin(2*CH_C_PI/numStaves)*cylRad*1.1, wallThickness, cylLen,
			(ChVector<double>*)&((ChVector<double>(-sin(ang),cos(ang),0))*cylRad),
			(ChMatrix33<double>*)&((ChMatrix33<double>)Q_from_AngZ(ang))
			);		
	}
	barrel->SetCollide(true);
	barrel->SetFriction(friction);
	barrel->AddCollisionModelsToSystem();


	ChSharedBodyPtr dividingWall(new ChBody);
	mphysicalSystem.AddBody(dividingWall);
	dividingWall->GetCollisionModel()->AddBox(0.002/2.0, cylRad*2.0/2.0, cylLen/2.0);
	dividingWall->SetPos(VNULL);
	dividingWall->SetBodyFixed(true);
	dividingWall->GetCollisionModel()->SetFamily(1);
	dividingWall->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	dividingWall->SetFriction(0.0f);


	// Prepare the physical system for the simulation
	mphysicalSystem.SetIntegrationType(ChSystem::INT_ANITESCU);
	mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
	mphysicalSystem.SetIterLCPmaxItersSpeed(iterations);
	mphysicalSystem.SetMaxPenetrationRecoverySpeed(0.5);

	ChCollisionModel::SetDefaultSuggestedEnvelope(particleRad*0.01);
	ChCollisionModel::SetDefaultSuggestedMargin(  particleRad*0.03);


	// OK! NOW GET READY FOR THE DYNAMICAL SIMULATION!
	GetLog() << "Running Simulation...\n";

	// A very simple simulation loop..
	double simTime(0);
	int currFrame(0);
	bool dividing(true);

	while (simTime < simDuration)
	{
		simTime += timestep; 

		// PERFORM SIMULATION UP TO chronoTime
		mphysicalSystem.DoFrameDynamics(simTime);

		if (currFrame%everyFrame == 0 )
		{

			// Output data to files
			char padnumber[256];
			sprintf(padnumber, "%d", (currFrame/everyFrame + 10000));
			char filename[256];
			// filepath (/tumbledata) must exist before executing
			sprintf(filename, "%s/pos%s.txt", "tumbledata", padnumber + 1);

			// create output file
			ChStreamOutAsciiFile output(filename);

			// Keep barrel rotation constant
			barrel->SetRot(Q_from_AngZ(-rpm*2*CH_C_PI/60*mphysicalSystem.GetChTime()));
			barrel->SetRot_dt(Q_from_AngZ(-rpm*2*CH_C_PI/60));

			// Keep barrel center of gravity fixed
			barrel->SetPos(VNULL);
			barrel->SetPos_dt(VNULL);
			barrel->SetPos_dtdt(VNULL);

			if (mphysicalSystem.GetChTime() > 0.1 & dividing)
			{
				dividing = false;
				dividingWall->SetCollide(false);
			}

			// i=3 => don't save barrel data or front/back walls (size()-2 ... dividing wall)
			for (int i = 3; i < mphysicalSystem.Get_bodylist()->size()-2; i++) {
				ChBody* abody = (ChBody*) mphysicalSystem.Get_bodylist()->at(i);
				ChVector<> bodPos = abody->GetPos();
				ChQuaternion<> bodRot = abody->GetRot();
				output << (float)bodPos.x << ",";
				output << (float)bodPos.y << ",";
				output << (float)bodPos.z << ",";
				output << (float)bodRot.e0 << ",";
				output << (float)bodRot.e1 << ",";
				output << (float)bodRot.e2 << ",";
				output << (float)bodRot.e3 << ",\n";
			}	
		}
		currFrame++;
		GetLog() << "Completed Step " << currFrame << "\n";
	}

	GetLog() << "DONE\n";


	// Remember this at the end of the program, if you started
	// with DLL_CreateGlobals();
	DLL_DeleteGlobals();

	return 0;
}
