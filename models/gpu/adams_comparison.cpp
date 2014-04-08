#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 1;
real tolerance = 1e-8;
real recovery_speed = 10000;
//#define USEPARALLEL

//#ifdef USEPARALLEL
//#define CHBODY new ChBody(new ChCollisionModelParallel)
//#define CHSYSTEM ChSystemParallel
//#else

//#define CHBODY new ChBody()
//#define CHSYSTEM ChSystem
//#endif
int max_iter = 100;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(6, 6, 6);
real container_thickness = .4;
real container_height = 0;
real container_friction = 1;

real particle_radius = .2;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = 0;
Vector particle_initial_vel = Vector(0, -5.5, 0);     //initial velocity
ChSharedBodyPtr R, REF;
ChSharedPtr<ChLinkEngine> engine_anchor;
template<class T>
void RunTimeStep(T* mSys, const int frame) {

	Vector force = engine_anchor->Get_react_force();
	Vector torque = engine_anchor->Get_react_torque();
	double motor_torque = engine_anchor->Get_mot_torque();
	cout << force.x << " " << force.y << " " << force.z << " " << torque.x << " " << torque.y << " " << torque.z << " " << motor_torque << endl;

	R->SetPos(REF->GetPos());
	R->SetPos_dt(REF->GetPos_dt());
	R->SetRot(REF->GetRot());
	R->SetWvel_loc(REF->GetWvel_loc());
}

int main(int argc, char* argv[]) {
	omp_set_num_threads(4);
	bool useparallel = false;

	useparallel = atoi(argv[1]);
	ChSystem * system_gpu;
//=========================================================================================================
	if (useparallel) {
		system_gpu = new ChSystemParallelDVI;
	} else {
		system_gpu = new ChSystem;
		system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);
		system_gpu->SetLcpSolverType(ChSystem::LCP_ITERATIVE_SYMMSOR);
		//system_gpu->GetLcpSolverSpeed()->SetVerbose(true);
		((ChLcpIterativeSolver*) system_gpu->GetLcpSolverSpeed())->SetRecordViolation(true);
		//((ChLcpIterativeSolver*) system_gpu->GetLcpSolverSpeed())->SetDiagonalPreconditioning(false);
		system_gpu->SetIterLCPwarmStarting(false);
	}
//=========================================================================================================
	if (useparallel) {
		((ChSystemParallelDVI*) system_gpu)->DoThreadTuning(false);
		((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(0);
		((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(30);
		((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
		((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(0);
		((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(tolerance);
		((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(tolerance);
		((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(1);
		((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
		//((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->DoStabilization(true);
		((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(.1 * .05);
		((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(50, 50, 50));
		((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(100, 50);
	}
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	system_gpu->SetIterLCPmaxItersStab(max_iter);
	system_gpu->SetTol(tolerance);
	system_gpu->SetTolSpeeds(tolerance);

	system_gpu->SetMaxPenetrationRecoverySpeed(1e9);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(0);
	material->SetRollingFriction(0);
	material->SetSpinningFriction(0);
	material->SetCompliance(0);
	material->SetCohesion(0);

	ChSharedBodyPtr L, BB;
	if (useparallel) {
		L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
		R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
		REF = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
		BB = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	} else {

		L = ChSharedBodyPtr(new ChBody());
		R = ChSharedBodyPtr(new ChBody());
		REF = ChSharedBodyPtr(new ChBody());
		BB = ChSharedBodyPtr(new ChBody());
	}

	InitObject(L, 1, Vector(0, 0, 0), Quaternion(1, 0, 0, 0), material, false, true);
	InitObject(R, 1000, Vector(0, 0, 0), Quaternion(1, 0, 0, 0), material, true, false);
	InitObject(REF, 1000, Vector(0, 0, 0), Quaternion(1, 0, 0, 0), material, false, false);
	InitObject(BB, 1000, Vector(2.5, 0, 0), Quaternion(1, 0, 0, 0), material, true, false);

	AddCollisionGeometry(L, SPHERE, Vector(.1, 1., 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(.1, 5, 2), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(REF, BOX, Vector(.1, 5, 2), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(BB, SPHERE, Vector(1, 3, 3), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	if (useparallel) {
		FinalizeObject(L, (ChSystemParallelDVI *) system_gpu);
		FinalizeObject(R, (ChSystemParallelDVI *) system_gpu);
		FinalizeObject(REF, (ChSystemParallelDVI *) system_gpu);
		//FinalizeObject(BB, (ChSystemParallelDVI *) system_gpu);

	} else {
		FinalizeObject(L, (ChSystem *) system_gpu);
		FinalizeObject(R, (ChSystem *) system_gpu);
		FinalizeObject(REF, (ChSystem *) system_gpu);
		FinalizeObject(BB, (ChSystem *) system_gpu);

		Vector r(.1, 5, 2);
		real mass = 1000;

		real3 inr = R3(1 / 12.0 * mass * (r.y * r.y + r.z * r.z), 1 / 12.0 * mass * (r.x * r.x + r.z * r.z), 1 / 12.0 * mass * (r.x * r.x + r.y * r.y));

		R->SetInertiaXX(ChVector<>(inr.x, inr.y, inr.z));
		REF->SetInertiaXX(ChVector<>(inr.x, inr.y, inr.z));
		real radius = 1;
		mass = 1000;
		BB->SetInertiaXX(ChVector<>(2 / 5.0 * mass * radius * radius, 2 / 5.0 * mass * radius * radius, 2 / 5.0 * mass * radius * radius));

	}
	 engine_anchor = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	engine_anchor->Initialize(L, R, ChCoordsys<>(R->GetPos(), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X)));
	engine_anchor->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);     // also works as revolute support
	engine_anchor->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);

	system_gpu->AddLink(engine_anchor);
	if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(engine_anchor->Get_spe_funct())) {
		mfun->Set_yconst(1);     // rad/s  angular speed
	}

	ChSharedPtr<ChLinkEngine> engine_ref = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	engine_ref->Initialize(L, REF, ChCoordsys<>(R->GetPos(), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X)));
	engine_ref->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);     // also works as revolute support
	engine_ref->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);

	system_gpu->AddLink(engine_ref);
	if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(engine_ref->Get_spe_funct())) {
		mfun->Set_yconst(1);     // rad/s  angular speed
	}

	ChSharedBodyPtr CONTAINER;

	if (useparallel) {
		CONTAINER = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	} else {
		CONTAINER = ChSharedBodyPtr(new ChBody());
	}

	InitObject(CONTAINER, 1, Vector(0, 0, 0), Quaternion(1, 0, 0, 0), material, true, true);

	AddCollisionGeometry(
			CONTAINER,
			BOX,
			Vector(container_thickness, container_size.y, container_size.z),
			Vector(-container_size.x + container_thickness, container_height - container_thickness, 0),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(
			CONTAINER,
			BOX,
			Vector(container_thickness, container_size.y, container_size.z),
			Vector(container_size.x - container_thickness, container_height - container_thickness, 0),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(
			CONTAINER,
			BOX,
			Vector(container_size.x, container_size.y, container_thickness),
			Vector(0, container_height - container_thickness, -container_size.z + container_thickness),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(
			CONTAINER,
			BOX,
			Vector(container_size.x, container_size.y, container_thickness),
			Vector(0, container_height - container_thickness, container_size.z - container_thickness),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(
			CONTAINER,
			BOX,
			Vector(container_size.x, container_thickness, container_size.z),
			Vector(0, container_height - container_size.y, 0),
			Quaternion(1, 0, 0, 0));

	if (useparallel) {
		FinalizeObject(CONTAINER, (ChSystemParallelDVI *) system_gpu);
	} else {
		FinalizeObject(CONTAINER, (ChSystem *) system_gpu);
	}
	if (useparallel) {
//		//ReadAllObjectsWithGeometryChrono(system_gpu, "bilateral_sphere_set.dat");
		ParticleGenerator<ChSystemParallelDVI> *layer_gen;
		layer_gen = new ParticleGenerator<ChSystemParallelDVI>((ChSystemParallelDVI *) system_gpu);
		layer_gen->SetDensity(particle_density);
		layer_gen->SetRadius(R3(particle_radius, particle_radius * .5, particle_radius));
		layer_gen->material->SetFriction(.1);
		layer_gen->material->SetCohesion(0);
		layer_gen->material->SetRollingFriction(0);
		layer_gen->material->SetSpinningFriction(0);
		layer_gen->AddMixtureType(MIX_SPHERE);

		layer_gen->addPerturbedVolumeMixture(R3(2, -2, 0), I3(8, 10, 20), R3(0, 0, 0), R3(0, 0, 0));
		layer_gen->addPerturbedVolumeMixture(R3(-2, -2, 0), I3(8, 10, 20), R3(0, 0, 0), R3(0, 0, 0));
	} else {
//		ParticleGenerator<ChSystem> *layer_gen;
//		layer_gen = new ParticleGenerator<ChSystem>((ChSystem *) system_gpu, false);
//		layer_gen->SetDensity(particle_density);
//		layer_gen->SetRadius(R3(particle_radius, particle_radius * .5, particle_radius));
//		layer_gen->material->SetFriction(.1);
//		layer_gen->material->SetCohesion(0);
//		layer_gen->material->SetRollingFriction(0);
//		layer_gen->material->SetSpinningFriction(0);
//		layer_gen->AddMixtureType(MIX_SPHERE);
//
//		layer_gen->addPerturbedVolumeMixture(R3(2, -2, 0), I3(8, 10, 20), R3(0, 0, 0), R3(0, 0, 0));
//		layer_gen->addPerturbedVolumeMixture(R3(-2, -2, 0), I3(8, 10, 20), R3(0, 0, 0), R3(0, 0, 0));

	}
//=========================================================================================================
//Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_position = glm::vec3(0, -5, -10);
	openGLView.render_camera->camera_look_at = glm::vec3(0, 0, 0);
	openGLView.render_camera->camera_scale = .5;
	openGLView.SetCustomCallback(RunTimeStep);
	openGLView.StartSpinning(window_manager);
	window_manager->CallGlutMainLoop();

//=========================================================================================================
	int file = 0;
	for (int i = 0; i < num_steps; i++) {
		system_gpu->DoStepDynamics(timestep);
		double TIME = system_gpu->GetChTime();
		double STEP = system_gpu->GetTimerStep();
		double BROD = system_gpu->GetTimerCollisionBroad();
		double NARR = system_gpu->GetTimerCollisionNarrow();
		double LCP = system_gpu->GetTimerLcp();
		double UPDT = system_gpu->GetTimerUpdate();
		std::vector<double> violation = ((ChLcpIterativeSolver*) system_gpu->GetLcpSolverSpeed())->GetViolationHistory();
		int REQ_ITS = violation.size();
		double RESID = 0;
		if (REQ_ITS != 0) {
			RESID = violation.at(violation.size() - 1);
		}
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

		RunTimeStep(system_gpu, i);
	}
	//DumpAllObjectsWithGeometryChrono(system_gpu, "bilateral_sphere_set.dat");
}
