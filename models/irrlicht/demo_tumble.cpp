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
#include "physics/CHconveyor.h"
#include "irrlicht_interface/CHbodySceneNode.h"
#include "irrlicht_interface/CHbodySceneNodeTools.h" 
#include "irrlicht_interface/CHirrAppInterface.h"
#include "irrlicht_interface/CHparticlesSceneNode.h"
#include "core/CHrealtimeStep.h"
#include "parallel/CHopenMP.h"
#include "omp.h"


// Use the namespace of Chrono

using namespace chrono;

// Use the main namespaces of Irrlicht
using namespace irr;

using namespace core;
using namespace scene; 
using namespace video;
using namespace io; 
using namespace gui; 

double scale(50); 
double friction(0.4f);

int numParticles(100);
double particleRad(1500*scale); // micrometres 63-180 must be less than 300 to see patterns
double particleDensity(0.005); // g/cm^3
 
double cylRad(0.10*scale); // m 0.10
double cylLen(0.0030*2.5*scale); // m 0.14
double wallThickness(0.005*scale); // m (0.005)
int numStaves(128);
double rpm(8); // 8 rpm

double timestep(0.001);
int iterations(25+numParticles/1000*5);

void createBodies(ChSystem& mphysicalSystem, ISceneManager* msceneManager, IVideoDriver* driver)
{
	// Create rigid bodies
	GetLog() << "Creating bodies...";

	// truss
	ChSharedBodyPtr truss(new ChBody);   
	mphysicalSystem.AddBody(truss);
	truss->SetBodyFixed(true);


	ChBodySceneNode* marker = (ChBodySceneNode*)addChBodySceneNode_easyBox(
			&mphysicalSystem, msceneManager,
			0.1,
			ChVector<double>(0,-1,0)*(cylRad+wallThickness),
			QNULL,
			ChVector<double>(sin(2*CH_C_PI/numStaves)*cylRad*1.1, wallThickness, cylLen) );
	marker->GetBody()->SetBodyFixed(true);
	marker->setMaterialTexture(0, driver->getTexture("../data/red.png"));
	marker->GetBody()->SetCollide(false);
	
	ChBodySceneNode* backWall = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem, msceneManager,
		1.0,
		ChVector<double>(0,0,cylLen/2.0+wallThickness/2.0),
		QNULL,
		ChVector<double>(cylRad*2.5,cylRad*2.5, wallThickness) );

	ChBodySceneNode* frontWall = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem, msceneManager,
		1.0,
		ChVector<double>(0,0,-cylLen/2.0-wallThickness/2.0),
		QNULL,
		ChVector<double>(cylRad*2.5,cylRad*2.5, wallThickness) );

	backWall->GetBody()->SetBodyFixed(true);
	frontWall->GetBody()->SetBodyFixed(true);

	frontWall->GetBody()->GetCollisionModel()->SetFamily(1);
	frontWall->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	backWall->GetBody()->GetCollisionModel()->SetFamily(1);
	backWall->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	
	backWall ->GetBody()->SetFriction(0.0f);
	frontWall->GetBody()->SetFriction(0.0f);

	frontWall->setVisible(false);
	backWall->setVisible(false);
	
	
	// Particles
	video::ITexture* red = driver->getTexture("../data/red.png");
	video::ITexture* green = driver->getTexture("../data/green.png");
	double particleMass = particleDensity/1000000 * 4/3*CH_C_PI*pow(particleRad/1000000,3);
	double particleInertia = particleRad/1000000*particleRad/1000000*particleMass;

	for (int np = 0; np <numParticles/2; np++) 
	{
		ChBodySceneNode* greenparticle = (ChBodySceneNode*)addChBodySceneNode_easySphere(
			&mphysicalSystem, msceneManager,
			particleMass,
			ChVector<>((CHrandom()-1)*1.3/2.0*cylRad-0.004, (CHrandom()-0.9)*0.7*cylRad, (CHrandom()-0.5)*cylLen), // position
			particleRad/1000000 // radius
			);
		greenparticle->setMaterialTexture(0,green);
		greenparticle->GetBody()->SetInertiaXX(ChVector<>(particleInertia,particleInertia,particleInertia));
		greenparticle->GetBody()->GetMaterialSurface()->SetFriction(friction);
		greenparticle->GetBody()->GetMaterialSurface()->SetRollingFriction(0.0);
	}


	for (int np = 0; np <numParticles/2; np++) 
	{
		ChBodySceneNode* redparticle = (ChBodySceneNode*)addChBodySceneNode_easySphere(
			&mphysicalSystem, msceneManager,
			particleMass,
			ChVector<double>((CHrandom())*1.3/2.0*cylRad+0.004, (CHrandom()-0.9)*0.7*cylRad, (CHrandom()-0.5)*cylLen), // position
			particleRad/1000000 // radius
			);
		redparticle->setMaterialTexture(0,red);
		redparticle->GetBody()->SetInertiaXX(ChVector<>(particleInertia,particleInertia,particleInertia));
		redparticle->GetBody()->GetMaterialSurface()->SetFriction(friction);
		redparticle->GetBody()->GetMaterialSurface()->SetRollingFriction(0.0);
	}

	GetLog() << "DONE\n";
}


int main(int argc, char* argv[])
{
	// Set max number of threads
	omp_set_num_threads(64);

	// In CHRONO engine, The DLL_CreateGlobals() - DLL_DeleteGlobals(); pair is needed if
	// global functions are needed. 
	DLL_CreateGlobals();


	// Create a ChronoENGINE physical system
	ChSystem mphysicalSystem;
	
	mphysicalSystem.Set_G_acc(ChVector<double>(0,-9.81*scale,0));

	// Create the Irrlicht visualization (open the Irrlicht device, 
	// bind a simple user interface, etc. etc.)
	ChIrrAppInterface application(&mphysicalSystem, L"Collisions between objects",core::dimension2d<u32>(1280,720),false);

	// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	ChIrrWizard::add_typical_Logo(application.GetDevice());
	ChIrrWizard::add_typical_Sky(application.GetDevice());
	ChIrrWizard::add_typical_Lights(application.GetDevice());
	ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(0,-cylRad/3.0,-cylRad*1.5), core::vector3df(0,-cylRad/3.0,0));


	// Create all the rigid bodies.
	createBodies(mphysicalSystem, application.GetSceneManager(), application.GetVideoDriver());


		// barrel
	ChSharedBodyPtr barrel(new ChBody);
	mphysicalSystem.AddBody(barrel);

	/*
	ChSharedPtr<ChLinkEngine> eng(new ChLinkEngine);
	eng->Initialize(barrel, truss, ChCoordsys<>(VNULL) );
	eng->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
	eng->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng->Get_spe_funct()))
		mfun->Set_yconst(-rpm*2*CH_C_PI/60); // rad/s  angular speed	:8
	mphysicalSystem.AddLink(eng);
	*/

	for (int i = 0; i < numStaves; i++) {
		double ang = 2*CH_C_PI/numStaves * i;
		
		/*
		ChSharedPtr<ChConveyor> mconveyor (new ChConveyor(sin(2*CH_C_PI/numStaves)*cylRad*1.1, wallThickness, cylLen));
		mconveyor->SetBodyFixed(true);
		mconveyor->SetFriction(friction);
		
		mconveyor->SetConveyorSpeed(rpm*2*CH_C_PI/60*cylRad*x); // !!! one half of belts are backward (not just first or last half, maybe first and last quarter?)!!!
		mconveyor->SetRot(Q_from_AngAxis(ang, VECT_Z) ); // ang-.0001 stuck pushed back in, !!!cannot slide down/ratchet
		mconveyor->SetPos( ChVector<double>(-sin(ang),cos(ang),0)*cylRad);

		mphysicalSystem.Add(mconveyor);
		*/

		/*
		double cl = cylLen;
		
		if ((i*4)%numStaves)
		{
			cl = cylLen*1.01;
		}
		*/
		
		/*
		barrel = (ChBodySceneNode*)addChBodySceneNode_easyBox(
			&mphysicalSystem, msceneManager,
			0.1,
			ChVector<double>(-sin(ang),cos(ang),0)*cylRad,
			Q_from_AngAxis(ang,VECT_Z), 
			ChVector<double>(sin(2*CH_C_PI/numStaves)*cylRad*1.1, wallThickness, cl) );
		*/
		
		barrel->GetCollisionModel()->AddBox(
			sin(2*CH_C_PI/numStaves)*cylRad*1.1, wallThickness, cylLen,
			(ChVector<double>*)&((ChVector<double>(-sin(ang),cos(ang),0))*cylRad),
			(ChMatrix33<double>*)&((ChMatrix33<double>)Q_from_AngZ(ang))
			);

		/*
		if ((i*4)%numStaves==0)
		{
			barrel->setMaterialTexture(0, driver->getTexture("../data/red.png"));
		}
		*/

		/*
		//barrel->GetCollisionModel()->AddBox(sin(ang)*cylRad, 0.01, cylLen, &(ChVector<double>(sin(ang),cos(ang),0)*cylRad), &ChMatrix33<double>(Q_from_AngAxis(ang,VECT_Z)));
		// custom collision callback
		barrel->GetBody()->GetCollisionModel()->SetFamily(1);
		barrel->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
		barrel->GetBody()->SetCollide(true);	
		barrel->GetBody()->SetFriction(friction);
		barrel->GetBody()->SetRollingFriction(0.0);
		*/
	}
	barrel->SetCollide(true);
	barrel->SetFriction(friction);
	barrel->AddCollisionModelsToSystem();

	ChBodySceneNode* dividingWall = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem, application.GetSceneManager(),
		1.0,
		VNULL,
		QNULL,
		ChVector<double>(0.002,cylRad*2.0, cylLen) );
		dividingWall->GetBody()->SetBodyFixed(true);
	dividingWall->GetBody()->GetCollisionModel()->SetFamily(1);
	dividingWall->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	dividingWall->setVisible(false);
	dividingWall->GetBody()->SetFriction(0.0f);


	// Modify some setting of the physical system for the simulation, if you want
	mphysicalSystem.SetIntegrationType(ChSystem::INT_ANITESCU);
	mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
	mphysicalSystem.SetIterLCPmaxItersSpeed(iterations);
	mphysicalSystem.SetMaxPenetrationRecoverySpeed(0.5);

	//mphysicalSystem.SetUseSleeping(true);

	application.SetStepManage(true);
	application.SetTryRealtime(true);
	application.SetTimestep(timestep);

	ChCollisionModel::SetDefaultSuggestedEnvelope(particleRad*0.01);
	ChCollisionModel::SetDefaultSuggestedMargin(  particleRad*0.03);

	// 
	// THE SOFT-REAL-TIME CYCLE
	//

	int loop(0);
	bool dividing = true;

	while(application.GetDevice()->run()) 
	{
		application.GetVideoDriver()->beginScene(true, true, SColor(255,140,161,192));

		application.DrawAll();

		application.DoStep();
		

		barrel->SetRot(Q_from_AngZ(-rpm*2*CH_C_PI/60*mphysicalSystem.GetChTime()));
		barrel->SetRot_dt(Q_from_AngZ(-rpm*2*CH_C_PI/60));
		barrel->SetPos(VNULL);
		barrel->SetPos_dt(VNULL);
		barrel->SetPos_dtdt(VNULL);
		
		
		if (mphysicalSystem.GetChTime() > 0.1 & dividing)
		{
			dividing = false;
			dividingWall->remove();
		}

		application.GetVideoDriver()->endScene(); 
		loop++;
	}



	// Remember this at the end of the program, if you started
	// with DLL_CreateGlobals();
	DLL_DeleteGlobals();

	return 0;
}