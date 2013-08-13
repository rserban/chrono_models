///////////////////////////////////////////////////
//
//   Demo code about  
//
//     - creating a pendulum
//     - applying torque to a body
//     - applying changes to a body during
//       simulation loop
//     - 3D viewing with the Irrlicht library
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
#include "core/CHtimer.h"
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

// Some parameters to play around with
double adjWeightHeight(0.5);	// height of adjustable weight from center of base (not where pin is)
double spacing(1.5);			// space scaling between metronomes (and rollers)
int numMetros(3);				// number of metronomes

double freq(10);	// changes frequency of metronome by increasing spring constant and restoring force
// !!! Too large causes significant moment in table leading to bouncing !!!


class metronome
{

public:
	ChBodySceneNode* metroArm;

	metronome(ChIrrAppInterface& application, int i, ChBodySceneNode* table, ISceneNode* metroParent, ISceneNode* armParent) 
	{
		// driver used for textures
		IVideoDriver* driver = application.GetVideoDriver();
		
		// texture names
		video::ITexture* goldMap = driver->getTexture("../data/gold.png");
		video::ITexture* woodMap = driver->getTexture("../data/wood.jpg");

		//
		// building a metronome
		//

		// .. metronome base
		ChBodySceneNode* base = (ChBodySceneNode*)addChBodySceneNode_easyBox(
			application.GetSystem(), application.GetSceneManager(),
			0.0,
			ChVector<>(i*(6*spacing),0.251,0),
			QUNIT, 
			ChVector<>(2.5,6.249,1),
			metroParent); 
			
		base->setMaterialTexture(0, woodMap);
		base->GetBody()->SetCollide(true);
		base->GetBody()->SetFriction(0.2f);
		base->addShadowVolumeSceneNode(); 


		// .. metronome arm
		metroArm = (ChBodySceneNode*)addChBodySceneNode_easyBox(
			application.GetSystem(), application.GetSceneManager(),
			0.00,
			ChVector<>(i*(6*spacing),0.5,-0.625),
			QUNIT, 
			ChVector<>(0.25,6,0.2),
			armParent );
			
		metroArm->GetBody()->SetCollide(true);
		metroArm->GetBody()->SetFriction(0.0f);
		metroArm->addShadowVolumeSceneNode(); 
		// initial angular velocity
		//	metroArm->GetBody()->SetRot_dt(ChQuaternion<>(1,0,0,8));


		// .. metronome static weight
		ChBodySceneNode* staticWeight = (ChBodySceneNode*)addChBodySceneNode_easyCylinder(
			application.GetSystem(), application.GetSceneManager(),
			2.0,
			ChVector<>(i*(6*spacing),-2,-0.75),
			chrono::Q_from_AngAxis(CH_C_PI/2, VECT_X),
			ChVector<>(1,0.25,1),
			metroParent );

		staticWeight->setMaterialTexture(0, goldMap);
		staticWeight->GetBody()->SetCollide(false);
		staticWeight->GetBody()->SetFriction(0.2f);
		staticWeight->addShadowVolumeSceneNode();


		// .. metronome adjustable weight
		ChBodySceneNode* adjWeight = (ChBodySceneNode*)addChBodySceneNode_easyBox(
			application.GetSystem(), application.GetSceneManager(),
			0.05,
			ChVector<>(i*(6*spacing),adjWeightHeight,-0.75),
			QUNIT, 
			ChVector<>(0.75,0.5,0.25),
			metroParent );

		adjWeight->setMaterialTexture(0, goldMap);
		adjWeight->GetBody()->SetCollide(false);
		adjWeight->GetBody()->SetFriction(0.2f);
		adjWeight->addShadowVolumeSceneNode();

		/*
		//	Create Stoppers
		// .. stoppers to limit magnitude of metronome action

		// .. Left Stopper
		ChBodySceneNode* stopper1 = (ChBodySceneNode*)addChBodySceneNode_easyCylinder(
		application.GetSystem(), application.GetSceneManager(),
			0.0,
			ChVector<>(i*(6*spacing)-1,1,-0.5),
			chrono::Q_from_AngAxis(CH_C_PI/2, VECT_X),
			ChVector<>(0.05,0.25,0.05),
			metroParent );

		stopper1->setMaterialTexture(0, goldMap);
		stopper1->GetBody()->SetCollide(true);
		stopper1->GetBody()->SetFriction(0.0f);
		stopper1->addShadowVolumeSceneNode();
		stopper1->GetBody()->SetImpactC(1.0f);

		// .. Right Stopper
		ChBodySceneNode* stopper2 = (ChBodySceneNode*)addChBodySceneNode_easyCylinder(
			application.GetSystem(), application.GetSceneManager(),
			0.0,
			ChVector<>(i*(6*spacing)+1,1,-0.5),
			chrono::Q_from_AngAxis(CH_C_PI/2, VECT_X),
			ChVector<>(0.05,0.25,0.05),
			metroParent );

		stopper2->setMaterialTexture(0, goldMap);
		stopper2->GetBody()->SetCollide(true);
		stopper2->GetBody()->SetFriction(0.0f);
		stopper2->addShadowVolumeSceneNode();
		stopper2->GetBody()->SetImpactC(1.0f);
		*/


		// 
		// Create the links between bodies (declare, specify, add)
		// 

		// .. pin arm to base
		ChSharedPtr<ChLinkLockRevolute> my_link_armbase(new ChLinkLockRevolute);
		my_link_armbase->Initialize(metroArm->GetBody(), base->GetBody(), ChCoordsys<>(ChVector<>(i*(6*spacing),-0.5,0.5)));
		application.GetSystem()->AddLink(my_link_armbase);

		// .. fix static weight to arm
		ChSharedPtr<ChLinkLockLock> my_link_staticarm(new ChLinkLockLock);
		my_link_staticarm->Initialize(staticWeight->GetBody(), metroArm->GetBody(), ChCoordsys<>(ChVector<>(i*(6*spacing),-2,-0.625)));
		application.GetSystem()->AddLink(my_link_staticarm);

		// .. fix adjustable weight to arm
		ChSharedPtr<ChLinkLockLock> my_link_adjarm(new ChLinkLockLock);
		my_link_adjarm->Initialize(adjWeight->GetBody(), metroArm->GetBody(), ChCoordsys<>(ChVector<>(i*(6*spacing),adjWeightHeight,-0.625)));
		application.GetSystem()->AddLink(my_link_adjarm);

		// .. fix metronome base to table
		ChSharedPtr<ChLinkLockLock> my_link_basetable(new ChLinkLockLock);
		my_link_basetable->Initialize(base->GetBody(), table->GetBody(), ChCoordsys<>(ChVector<>(i*(6*spacing),-3,0)));
		application.GetSystem()->AddLink(my_link_basetable);

		// .. add spring to metronome (increases frequency and adds damping)
		ChSharedPtr<ChLinkSpring> my_spring = ChSharedPtr<ChLinkSpring>(new ChLinkSpring);
		my_spring->Initialize(staticWeight->GetBody(), base->GetBody(), false, ChVector<>(i*(6*spacing),-2,-0.625), ChVector<>(i*(6*spacing)+0.0,-2,-0.5));
		my_spring->Set_SpringK(10*freq);	// Spring force
		my_spring->Set_SpringR(5);		// Spring resistance (damping) 5:CS
		my_spring->Set_SpringRestLenght(sqrt(pow(0.125,2)+pow(0.0,2))); // pythagorean using distance from static to base and distance offset from center(0)
		application.GetSystem()->AddLink(my_spring);

		/*
		//	attach stoppers
		// .. fix stoppers to base
		ChSharedPtr<ChLinkLockLock> my_link_stopper1(new ChLinkLockLock);
		my_link_stopper1->Initialize(stopper1->GetBody(), base->GetBody(), ChCoordsys<>(ChVector<>(i*(6*spacing)-1.2,-0.125,-0.75)));
		application.GetSystem()->AddLink(my_link_stopper1);

		ChSharedPtr<ChLinkLockLock> my_link_stopper2(new ChLinkLockLock);
		my_link_stopper2->Initialize(stopper2->GetBody(), base->GetBody(), ChCoordsys<>(ChVector<>(i*(6*spacing)+1.2,-0.125,-0.75)));
		application.GetSystem()->AddLink(my_link_stopper2);
		*/
	}

};

int main(int argc, char* argv[])
{

	// In CHRONO engine, The DLL_CreateGlobals() - DLL_DeleteGlobals(); pair is needed if
	// global functions are needed.
	DLL_CreateGlobals();

	// Create a ChronoENGINE physical system
	ChSystem my_system;

	// Specify solver and number of iterations
	my_system.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR);
	my_system.SetMaxiter(75*(int)pow(numMetros,1.5));
	my_system.SetIterLCPmaxItersSpeed(75*(int)pow(numMetros,1.5));

	// Create the Irrlicht visualization (open the Irrlicht device, 
	// bind a simple user interface, etc. etc.)
	ChIrrAppInterface application(&my_system, L"Metronome on a Sliding Table",core::dimension2d<u32>(1280,720),false,false, video::EDT_OPENGL); 

	// Easy shortcuts to add logo, camera, lights and sky in Irrlicht scene:
	ChIrrWizard::add_typical_Logo(application.GetDevice());
	ChIrrWizard::add_typical_Sky(application.GetDevice());
	ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(5,15,-5));
	ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(0,14,-20));

	// driver used for textures
	IVideoDriver* driver = application.GetVideoDriver();

	// texture names
	video::ITexture* concreteMap = driver->getTexture("../data/concrete.jpg");
	video::ITexture* ceilingMap = driver->getTexture("../data/blu.png");
	video::ITexture* cylinderMap = driver->getTexture("../data/bluwhite.png");

	// Create Irrlicht 'directories' where metronomes will be put during the simulation loop
	ISceneNode* metronomes = application.GetSceneManager()->addEmptySceneNode();
	ISceneNode* arms = application.GetSceneManager()->addEmptySceneNode();


	//
	// Create the rigid bodies
	//

	// .. floor
	ChBodySceneNode* floor = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&my_system, application.GetSceneManager(),
		1.0,
		ChVector<>(3*spacing*(numMetros-1),-9.6,0),
		QUNIT, 
		ChVector<>(6*spacing*(numMetros+7),1,20) ); 
		
	floor->GetBody()->SetBodyFixed(true);
	floor->GetBody()->SetCollide(true);
	floor->GetBody()->SetFriction(0.4f);
	// rolling friction: resisting torque generated. Lesser of two contact surfaces used.
	floor->GetBody()->SetRollingFriction(0.002f);
	floor->addShadowVolumeSceneNode();
	floor->setMaterialTexture(0, concreteMap);


	// .. bumper walls to limit table translation
	ChBodySceneNode* wall1 = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&my_system, application.GetSceneManager(),
		1.0,
		ChVector<>(3*spacing*(numMetros-1)-6*spacing*(numMetros+7)/2,-5.6,0),
		QUNIT, 
		ChVector<>(1,9,20) ); 
		
	wall1->GetBody()->SetBodyFixed(true);
	wall1->GetBody()->SetCollide(true);
	wall1->addShadowVolumeSceneNode();
	wall1->setMaterialTexture(0, concreteMap);

	ChBodySceneNode* wall2 = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&my_system, application.GetSceneManager(),
		1.0,
		ChVector<>(3*spacing*(numMetros-1)+6*spacing*(numMetros+7)/2,-5.6,0),
		QUNIT, 
		ChVector<>(1,9,20) ); 
		
	wall2->GetBody()->SetBodyFixed(true);
	wall2->GetBody()->SetCollide(true);
	wall2->addShadowVolumeSceneNode();
	wall2->setMaterialTexture(0, concreteMap);


	// .. sliding table
	ChBodySceneNode* table = (ChBodySceneNode*)addChBodySceneNode_easyBox(
		&my_system, application.GetSceneManager(),
		0.0,
		ChVector<>(3*spacing*(numMetros-1),-3,0),
		QUNIT, 
		ChVector<>(6*spacing*(numMetros+3),0.25,6) ); 
	table->GetBody()->SetCollide(true);
	
	table->GetBody()->SetFriction(0.4f);
	table->GetBody()->SetRollingFriction(0.002f);
	table->addShadowVolumeSceneNode();
	// fix the location of the table
	table->GetBody()->SetBodyFixed(true);


	//	create rollers
	// .. left table roller
	ChBodySceneNode* rollerLeft = (ChBodySceneNode*)addChBodySceneNode_easyCylinder(
		&my_system, application.GetSceneManager(),
		0.1,
		ChVector<>(-3*spacing,-6.5,0),
		chrono::Q_from_AngAxis(CH_C_PI/2, VECT_X),
		ChVector<>(5,12,5) );
		
	rollerLeft->GetBody()->SetCollide(true);
	rollerLeft->GetBody()->SetFriction(0.4f);
	rollerLeft->GetBody()->SetRollingFriction(0.002f);
	rollerLeft->addShadowVolumeSceneNode();
	rollerLeft->setMaterialTexture(0, driver->getTexture("../data/bluwhite.png"));


	// .. right table roller
	ChBodySceneNode* rollerRight = (ChBodySceneNode*)addChBodySceneNode_easyCylinder(
		&my_system, application.GetSceneManager(),
		0.1,
		ChVector<>(3*spacing*(2*numMetros),-6.5,0),
		chrono::Q_from_AngAxis(CH_C_PI/2, VECT_X),
		ChVector<>(5,12,5) );
		
	rollerRight->GetBody()->SetCollide(true);
	rollerRight->GetBody()->SetFriction(0.4f);
	rollerRight->GetBody()->SetRollingFriction(0.002f);
	rollerRight->addShadowVolumeSceneNode();
	rollerRight->setMaterialTexture(0, cylinderMap);


	// .. create metronomes
	for (int i = 0; i < numMetros; i++)
	{
		metronome metronome(application, i, table, metronomes, arms);
	}


	// Create a GUI button (later assigned to free the table)
	IGUIButton* freeButton = application.GetDevice()->getGUIEnvironment()->addButton(rect<s32>(150, 10, 250, 50), 0, 0, L"Release Table");


	//
	// THE SOFT-REAL-TIME CYCLE
	//

	// set time step settings
	application.SetStepManage(true);
	application.SetTimestep(0.02/numMetros);
	application.SetTryRealtime(true);

	while(application.GetDevice()->run())
	{
		// Irrlicht must prepare frame to draw
		application.GetVideoDriver()->beginScene(true, true, SColor(255,140,161,192));

		// Irrlicht application draws all 3D objects and all GUI items
		application.DrawAll();

		// HERE CHRONO INTEGRATION IS PERFORMED: THE 
		// TIME OF THE SIMULATION ADVANCES FOR A SINGLE
		// STEP:
		application.DoStep();

		// release table if button is ever pressed
		if (freeButton->isPressed())
		{
			table->GetBody()->SetBodyFixed(false);
		}

		// Modifying metronome during simulation
		ISceneNodeList::ConstIterator it = arms->getChildren().begin();
		ChBodySceneNode* anArm(0);
		for (; it != arms->getChildren().end(); it++)
		{
			ISceneNode* arm = (*it);
			anArm = (ChBodySceneNode*)arm;
			double rot = anArm->GetBody()->GetRot().GetVector().z;
			double vel = anArm->GetBody()->GetRot_dt().GetVector().z;

			// Escapement force
			// applies just after passing bottom of swing
			if (rot > -0.2 && rot < -0.01 && vel > 0)
			{ 
				anArm->GetBody()->Accumulate_torque(ChVector<>(0,0,freq/10.0), 1);
			}
			else if (rot > 0.01 && rot < 0.2 && vel < 0)
			{
				anArm->GetBody()->Accumulate_torque(ChVector<>(0,0,-freq/10.0), 1);
			}
			// clear accumualted torque while outside of range
			else
			{
				anArm->GetBody()->Empty_forces_accumulators();
			}

			// Color change to represent 'ticking'
			if (rot > 0.01 && fabs(vel) < 0.1)
			{
				anArm->setMaterialTexture(0, driver->getTexture("../data/red.png"));
			}
			else if (rot < 0.01 && fabs(vel) < 0.1)
			{
				// default texture: white
				anArm->setMaterialTexture(0,0);
			}
		}

		application.GetVideoDriver()->endScene(); 

	}

	// Remember this at the end of the program, if you started
	// with DLL_CreateGlobals();
	DLL_DeleteGlobals();

	return 0;
}