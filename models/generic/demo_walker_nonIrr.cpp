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
#include "../../common/common.h"
using namespace chrono;

double timestep(0.01);
double simDuration(20);
int everyFrame(10);
real gravity = -9.80665;
double chassisL(15.0);
double chassisW(2.0);
double legW(1.0);
double legL(2.0);
double footH(0.5);
double footW(1.5);
double footL(2.0);
double axleL(3.0);
template<class T>
void RunTimeStep(T* mSys, const int frame) {
	double TIME = mSys->GetChTime();
	double STEP = mSys->GetTimerStep();
	double BROD = mSys->GetTimerCollisionBroad();
	double NARR = mSys->GetTimerCollisionNarrow();
	double LCP = mSys->GetTimerLcp();
	double UPDT = mSys->GetTimerUpdate();
	int BODS = mSys->GetNbodies();
	int CNTC = mSys->GetNcontacts();
	int REQ_ITS = ((ChLcpSolverGPU*) (mSys->GetLcpSolverSpeed()))->GetTotalIterations();

	printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS);
}
int main(int argc, char* argv[]) {
	DLL_CreateGlobals();

	ChSystem * system = new ChSystem();
	system->SetMaxiter(100);
	system->SetIterLCPmaxItersSpeed(100);
	system->SetTol(0);
	system->SetTolSpeeds(0);

	ChLcpSystemDescriptor *mdescriptor = new ChLcpSystemDescriptor();
	ChContactContainer *mcontactcontainer = new ChContactContainer();
	//ChCollisionSystemBulletGPU *mcollisionengine = new ChCollisionSystemBulletGPU();
	ChCollisionSystemBullet *mcollisionengine = new ChCollisionSystemBullet();
	system->ChangeLcpSystemDescriptor(mdescriptor);
	system->ChangeContactContainer(mcontactcontainer);
	system->ChangeCollisionSystem(mcollisionengine);
	system->SetIntegrationType(ChSystem::INT_ANITESCU);
	system->SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
	system->SetStep(timestep);
	system->Set_G_acc(ChVector<>(0, gravity, 0));

	ChSharedBodyPtr floor(new ChBody);

	ChSharedBodyPtr chassis(new ChBody);

	ChSharedBodyPtr axle_F(new ChBody);
	ChSharedBodyPtr axle_C(new ChBody);
	ChSharedBodyPtr axle_R(new ChBody);

	ChSharedBodyPtr leg_FR(new ChBody);
	ChSharedBodyPtr leg_FL(new ChBody);
	ChSharedBodyPtr leg_CR(new ChBody);
	ChSharedBodyPtr leg_CL(new ChBody);
	ChSharedBodyPtr leg_RR(new ChBody);
	ChSharedBodyPtr leg_RL(new ChBody);

	InitObject(floor, 1.0, ChVector<>(0, -3.5, 0), Quaternion(1, 0, 0, 0), 1, 1, 0, true, true, -20, -20);
	InitObject(chassis, 1, ChVector<>(0, 0, 0), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, 0, 1);
	InitObject(axle_F, 1, ChVector<>(0, 0, chassisL / 2), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, false, false, -2, -2);
	InitObject(axle_C, 1, ChVector<>(0, 0, 0), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, false, false, -2, -2);
	InitObject(axle_R, 1, ChVector<>(0, 0, -chassisL / 2), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, false, false, -2, -2);
	InitObject(leg_FR, 1, ChVector<>((axleL + legW) / 2.0, -legL / 2.0, chassisL / 2), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_FL, 1, ChVector<>(-(axleL + legW) / 2.0, legL / 2.0, chassisL / 2), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_CR, 1, ChVector<>(-(axleL + legW) / 2.0, -legL / 2.0, 0), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_CL, 1, ChVector<>((axleL + legW) / 2.0, legL / 2.0, 0), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_RR, 1, ChVector<>((axleL + legW) / 2.0, -legL / 2.0, -chassisL / 2), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_RL, 1, ChVector<>(-(axleL + legW) / 2.0, legL / 2.0, -chassisL / 2), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, 2, 2);

	AddCollisionGeometry(floor, BOX, ChVector<>(10, 1 / 2.0, 10), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(chassis, ELLIPSOID, ChVector<>(1, 1, 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(axle_F, ELLIPSOID, ChVector<>(0.5 / 2.0, axleL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(axle_C, ELLIPSOID, ChVector<>(0.5 / 2.0, axleL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(axle_R, ELLIPSOID, ChVector<>(0.5 / 2.0, axleL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(leg_FR, ELLIPSOID, ChVector<>(legW / 2.0, legL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(leg_FL, ELLIPSOID, ChVector<>(legW / 2.0, legL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(leg_CR, ELLIPSOID, ChVector<>(legW / 2.0, legL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(leg_CL, ELLIPSOID, ChVector<>(legW / 2.0, legL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(leg_RR, ELLIPSOID, ChVector<>(legW / 2.0, legL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(leg_RL, ELLIPSOID, ChVector<>(legW / 2.0, legL / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	floor->SetInertiaXX(Vector(1,1,1));
	chassis->SetInertiaXX(Vector(1,1,1));
	axle_F->SetInertiaXX(Vector(1,1,1));
	axle_C->SetInertiaXX(Vector(1,1,1));
	axle_R->SetInertiaXX(Vector(1,1,1));
	leg_FR->SetInertiaXX(Vector(1,1,1));
	leg_FL->SetInertiaXX(Vector(1,1,1));
	leg_CR->SetInertiaXX(Vector(1,1,1));
	leg_CL->SetInertiaXX(Vector(1,1,1));
	leg_RR->SetInertiaXX(Vector(1,1,1));
	leg_RL->SetInertiaXX(Vector(1,1,1));

	FinalizeObject(floor, (ChSystem *) system);
	FinalizeObject(chassis, (ChSystem *) system);
	FinalizeObject(axle_F, (ChSystem *) system);
	FinalizeObject(axle_C, (ChSystem *) system);
	FinalizeObject(axle_R, (ChSystem *) system);
	FinalizeObject(leg_FR, (ChSystem *) system);
	FinalizeObject(leg_FL, (ChSystem *) system);
	FinalizeObject(leg_CR, (ChSystem *) system);
	FinalizeObject(leg_CL, (ChSystem *) system);
	FinalizeObject(leg_RR, (ChSystem *) system);
	FinalizeObject(leg_RL, (ChSystem *) system);

	ChSharedBodyPtr chassis_ptr = ChSharedBodyPtr(chassis);
	ChSharedBodyPtr axle_F_ptr = ChSharedBodyPtr(axle_F);
	ChSharedBodyPtr axle_C_ptr = ChSharedBodyPtr(axle_C);
	ChSharedBodyPtr axle_R_ptr = ChSharedBodyPtr(axle_R);
	ChSharedBodyPtr leg_FR_ptr = ChSharedBodyPtr(leg_FR);
	ChSharedBodyPtr leg_FL_ptr = ChSharedBodyPtr(leg_FL);
	ChSharedBodyPtr leg_CR_ptr = ChSharedBodyPtr(leg_CR);
	ChSharedBodyPtr leg_CL_ptr = ChSharedBodyPtr(leg_CL);
	ChSharedBodyPtr leg_RR_ptr = ChSharedBodyPtr(leg_RR);
	ChSharedBodyPtr leg_RL_ptr = ChSharedBodyPtr(leg_RL);

//	ChSharedPtr<ChLinkLockLock> axle_FR(new ChLinkLockLock);
//	axle_FR->Initialize(leg_FR_ptr, axle_F_ptr, ChCoordsys<>(VNULL));
//	system->AddLink(axle_FR);
////
//	ChSharedPtr<ChLinkLockLock> axle_FL(new ChLinkLockLock);
//	axle_FL->Initialize(leg_FL_ptr, axle_F_ptr, ChCoordsys<>(VNULL));
//	system->AddLink(axle_FL);
//
//	ChSharedPtr<ChLinkLockLock> axle_CR(new ChLinkLockLock);
//	axle_CR->Initialize(leg_CR_ptr, axle_C_ptr, ChCoordsys<>(VNULL));
//	system->AddLink(axle_CR);
//
//	ChSharedPtr<ChLinkLockLock> axle_CL(new ChLinkLockLock);
//	axle_CL->Initialize(leg_CL_ptr, axle_C_ptr, ChCoordsys<>(VNULL));
//	system->AddLink(axle_CL);
//
//	ChSharedPtr<ChLinkLockLock> axle_RR(new ChLinkLockLock);
//	axle_RR->Initialize(leg_RR_ptr, axle_R_ptr, ChCoordsys<>(VNULL));
//	system->AddLink(axle_RR);
//
//	ChSharedPtr<ChLinkLockLock> axle_RL(new ChLinkLockLock);
//	axle_RL->Initialize(leg_RL_ptr, axle_R_ptr, ChCoordsys<>(VNULL));
//	system->AddLink(axle_RL);

	// Create engine between axles and chassis
	GetLog() << "Creating Motors\n";

//	ChSharedPtr<ChLinkEngine> eng_F(new ChLinkEngine);
//	eng_F->Initialize(axle_F_ptr, chassis_ptr, ChCoordsys<>(ChVector<>(0, 0, chassisL / 2), Q_from_AngY(CH_C_PI / 2)));
//	eng_F->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
//	eng_F->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
//	if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_F->Get_spe_funct()))
//		mfun->Set_yconst(1); // rad/s  angular speed
//	system->AddLink(eng_F);

	ChSharedPtr<ChLinkEngine> eng_C(new ChLinkEngine);
	eng_C->Initialize(axle_C_ptr, chassis_ptr, ChCoordsys<>(ChVector<>(0, 0, 0), Q_from_AngY(CH_C_PI / 2)));
	eng_C->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
	eng_C->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_C->Get_spe_funct()))
		mfun->Set_yconst(1); // rad/s  angular speed
	system->AddLink(eng_C);

//	ChSharedPtr<ChLinkEngine> eng_R(new ChLinkEngine);
//	eng_R->Initialize(axle_R_ptr, chassis_ptr, ChCoordsys<>(ChVector<>(0, 0, -chassisL / 2), Q_from_AngY(CH_C_PI / 2)));
//	eng_R->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
//	eng_R->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
//	if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_R->Get_spe_funct()))
//		mfun->Set_yconst(1); // rad/s  angular speed
//	system->AddLink(eng_R);

	system->DoStepDynamics(timestep);
	system->DoStepDynamics(timestep);

	exit(0);

	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
	openGLView.render_camera->look_at = Vector(0, -5, 0);
	openGLView.render_camera->mScale = .5;
	openGLView.SetCustomCallback(RunTimeStep);
	openGLView.StartSpinning(window_manager);
	window_manager->CallGlutMainLoop();

	DLL_DeleteGlobals();

	return 0;
}

