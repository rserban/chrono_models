///////////////////////////////////////////////////
//
//   Demo code about  
//
//     - gear constraint (ChLinkGear) as a method
//       to impose a transmission ratio between two
//       shafts as they were connected by gears, without
//       the need of performing collision detection between
//       gear teeth geometries (which would be inefficient)
// 
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
#include "irrlicht_interface/CHbodySceneNode.h"
#include "irrlicht_interface/CHbodySceneNodeTools.h" 
#include "irrlicht_interface/CHirrAppInterface.h"
#include "core/CHrealtimeStep.h"

#include <irrlicht.h>



// Use the namespace of Chrono

using namespace chrono;

// Use the main namespaces of Irrlicht
using namespace irr;

using namespace core;
using namespace scene; 
using namespace video;
using namespace io; 
using namespace gui; 


int main(int argc, char* argv[])
{

	// In CHRONO engine, The DLL_CreateGlobals() - DLL_DeleteGlobals(); pair is needed if
	// global functions are needed. 
	DLL_CreateGlobals();

	// Create a ChronoENGINE physical system
	ChSystem mphysicalSystem;

	// Set up solver settings
	mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR);
	mphysicalSystem.SetMaxiter(100);
	mphysicalSystem.SetIterLCPmaxItersSpeed(100);
	mphysicalSystem.SetIterLCPmaxItersStab(100);

	// Create the Irrlicht visualization (open the Irrlicht device, 
	// bind a simple user interface, etc. etc.)
	ChIrrAppInterface application(&mphysicalSystem, L"Walker on Rigid Terrain",core::dimension2d<u32>(800,600),false); 


	// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	ChIrrWizard::add_typical_Logo(application.GetDevice());
	ChIrrWizard::add_typical_Sky(application.GetDevice());
	ChIrrWizard::add_typical_Lights(application.GetDevice());
	ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(12, 15,-20));


	// Create a floor
	ChBodySceneNode* floor = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		1.0,
		ChVector<>(0, -2.5, 400/2.-20),
		QUNIT,
		ChVector<>(40, 1, 400)
		);
	floor->GetBody()->SetBodyFixed(true);
	floor->GetBody()->SetCollide(true);
	floor->setMaterialTexture(0, application.GetVideoDriver()->getTexture("../data/concrete.jpg"));

	// chassis
	ChBodySceneNode* chassis = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		4.,
		VNULL,
		QUNIT,
		ChVector<>(2.,1.,15.) );
	
	// CREATE AXLES
	// front
	ChBodySceneNode* axle_F = (ChBodySceneNode*)addChBodySceneNode_easyCylinder(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		Vector(0,0,7),
		Q_from_AngZ(CH_C_PI/2.),
		ChVector<>(1.,3.,1.) );
	axle_F->GetBody()->SetInertiaXX(ChVector<>(0.1*0.5*0.5/2., 0.1*0.5*0.5/4 + 0.1/3., 0.1*0.5*0.5/4 + 0.1/3.));
	// center
	ChBodySceneNode* axle_C = (ChBodySceneNode*)addChBodySceneNode_easyCylinder(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		VNULL,
		Q_from_AngZ(CH_C_PI/2),
		ChVector<>(1.,3.,1.) );
	axle_C->GetBody()->SetInertiaXX(ChVector<>(0.1*0.5*0.5/2., 0.1*0.5*0.5/4 + 0.1/3., 0.1*0.5*0.5/4 + 0.1/3.));
	// rear
	ChBodySceneNode* axle_R = (ChBodySceneNode*)addChBodySceneNode_easyCylinder(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		Vector(0,0,-7),
		Q_from_AngZ(CH_C_PI/2),
		ChVector<>(1.,3.,1.) );
	axle_R->GetBody()->SetInertiaXX(ChVector<>(0.1*0.5*0.5/2., 0.1*0.5*0.5/4 + 0.1/3., 0.1*0.5*0.5/4 + 0.1/3.));
	
		// no axle collision
		axle_F->GetBody()->SetCollide(false);
		axle_C->GetBody()->SetCollide(false);
		axle_R->GetBody()->SetCollide(false);
	

	// CREATE LEGS
	// front right
	ChBodySceneNode* leg_FR = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(2., -1., 7),
		QUNIT,
		ChVector<>(1., 2., 0.5) );
	// front left
	ChBodySceneNode* leg_FL = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(-2., 1., 7),
		QUNIT,
		ChVector<>(1., 2., 0.5) );
	// center right
	ChBodySceneNode* leg_CR = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(-2., -1., 0),
		QUNIT,
		ChVector<>(1., 2., 0.5) );
	// center left
	ChBodySceneNode* leg_CL = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(2., 1., 0),
		QUNIT,
		ChVector<>(1., 2., 0.5) );
	// rear right
	ChBodySceneNode* leg_RR = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(2., -1., -7),
		QUNIT,
		ChVector<>(1., 2., 0.5) );
	// rear left
	ChBodySceneNode* leg_RL = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(-2., 1., -7),
		QUNIT,
		ChVector<>(1., 2., 0.5) );

	// CREATE FEET
	// front right
	ChBodySceneNode* foot_FR = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(2.25, -1.75, 8),
		QUNIT,
		ChVector<>(1.5, 0.5, 2.) );
	// front left
	ChBodySceneNode* foot_FL = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(-2.25, 1.75, 6),
		QUNIT,
		ChVector<>(1.5, 0.5, 2.) );
	// center right
	ChBodySceneNode* foot_CR = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(-2.25, -1.75, 1),
		QUNIT,
		ChVector<>(1.5, 0.5, 2.) );
	// center left
	ChBodySceneNode* foot_CL = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(2.25, 1.75, -1),
		QUNIT,
		ChVector<>(1.5, 0.5, 2.) );
	// rear right
	ChBodySceneNode* foot_RR = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(2.25, -1.75, -6),
		QUNIT,
		ChVector<>(1.5, 0.5, 2.) );
	// rear left
	ChBodySceneNode* foot_RL = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&mphysicalSystem,
		application.GetSceneManager(),
		0.1,
		ChVector<>(-2.25, 1.75, -8),
		QUNIT,
		ChVector<>(1.5, 0.5, 2.) );

	// ATTACH FEET TO LEGS
	ChSharedPtr<ChLinkLockLock> link_FR(new ChLinkLockLock); 
	link_FR->Initialize(foot_FR->GetBody(), leg_FR->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(link_FR);

	ChSharedPtr<ChLinkLockLock> link_FL(new ChLinkLockLock); 
	link_FL->Initialize(foot_FL->GetBody(), leg_FL->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(link_FL);

	ChSharedPtr<ChLinkLockLock> link_CR(new ChLinkLockLock); 
	link_CR->Initialize(foot_CR->GetBody(), leg_CR->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(link_CR);

	ChSharedPtr<ChLinkLockLock> link_CL(new ChLinkLockLock); 
	link_CL->Initialize(foot_CL->GetBody(), leg_CL->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(link_CL);

	ChSharedPtr<ChLinkLockLock> link_RR(new ChLinkLockLock); 
	link_RR->Initialize(foot_RR->GetBody(), leg_RR->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(link_RR);

	ChSharedPtr<ChLinkLockLock> link_RL(new ChLinkLockLock); 
	link_RL->Initialize(foot_RL->GetBody(), leg_RL->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(link_RL);
	
	// ATTACH LEGS TO AXLES
	ChSharedPtr<ChLinkLockLock> axle_FR(new ChLinkLockLock); 
	axle_FR->Initialize(leg_FR->GetBody(), axle_F->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(axle_FR);

	ChSharedPtr<ChLinkLockLock> axle_FL(new ChLinkLockLock); 
	axle_FL->Initialize(leg_FL->GetBody(), axle_F->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(axle_FL);

	ChSharedPtr<ChLinkLockLock> axle_CR(new ChLinkLockLock); 
	axle_CR->Initialize(leg_CR->GetBody(), axle_C->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(axle_CR);

	ChSharedPtr<ChLinkLockLock> axle_CL(new ChLinkLockLock); 
	axle_CL->Initialize(leg_CL->GetBody(), axle_C->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(axle_CL);

	ChSharedPtr<ChLinkLockLock> axle_RR(new ChLinkLockLock); 
	axle_RR->Initialize(leg_RR->GetBody(), axle_R->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(axle_RR);

	ChSharedPtr<ChLinkLockLock> axle_RL(new ChLinkLockLock); 
	axle_RL->Initialize(leg_RL->GetBody(), axle_R->GetBody(), ChCoordsys<>(VNULL));
	mphysicalSystem.AddLink(axle_RL);
	
	// CREATE ENGINE BETWEEN AXLES AND CHASSIS
	ChSharedPtr<ChLinkEngine> eng_F(new ChLinkEngine);
		eng_F->Initialize(axle_F->GetBody(), chassis->GetBody(), 
			ChCoordsys<>(ChVector<>(0, 0, 7) , Q_from_AngY(CH_C_PI/2) ) );
		eng_F->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
		eng_F->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_F->Get_spe_funct()))
			mfun->Set_yconst(1); // rad/s  angular speed
	mphysicalSystem.AddLink(eng_F);

	ChSharedPtr<ChLinkEngine> eng_C(new ChLinkEngine);
		eng_C->Initialize(axle_C->GetBody(), chassis->GetBody(), 
			ChCoordsys<>(ChVector<>(0, 0, 0) , Q_from_AngY(CH_C_PI/2) ) );
		eng_C->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
		eng_C->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_C->Get_spe_funct()))
			mfun->Set_yconst(1); // rad/s  angular speed
	mphysicalSystem.AddLink(eng_C);

	ChSharedPtr<ChLinkEngine> eng_R(new ChLinkEngine);
		eng_R->Initialize(axle_R->GetBody(), chassis->GetBody(), 
			ChCoordsys<>(ChVector<>(0, 0, -7) , Q_from_AngY(CH_C_PI/2) ) );
		eng_R->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
		eng_R->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_R->Get_spe_funct()))
			mfun->Set_yconst(1); // rad/s  angular speed
	mphysicalSystem.AddLink(eng_R);



	// Prepare the physical system for the simulation

	mphysicalSystem.SetIntegrationType(ChSystem::INT_ANITESCU); 


	// 
	// THE SOFT-REAL-TIME CYCLE
	//

	
	application.SetStepManage(true);
	application.SetTimestep(0.01);
	application.SetTryRealtime(true);
	

	while(application.GetDevice()->run()) 
	{
		application.GetVideoDriver()->beginScene(true, true, SColor(255,140,161,192));

		application.DrawAll();

		// ADVANCE THE SIMULATION FOR ONE STEP

		application.DoStep();

		application.GetVideoDriver()->endScene();  
	}


	// Remember this at the end of the program, if you started
	// with DLL_CreateGlobals();
	DLL_DeleteGlobals();

	return 0;
}

