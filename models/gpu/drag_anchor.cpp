#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 30;

int max_iter = 10;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(10, 3, 3);
real container_thickness = .2;
real container_height = 0;
real container_friction = 1;

real particle_radius = .025;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = .5;
real particle_cohesion = .1;
real rolling_fric = 1;
real spinning_fric = 1;

ParticleGenerator<ChSystemParallel>* layer_gen;
vector<ChSharedBodyPtr> string_vector(1000);
string data_folder = "data/drag_anchor";
real chain_radius = .1;
int fibers = 0;
ChSharedPtr<ChLinkEngine> eng_roller;
ChSharedBodyPtr Roller;
ChSharedPtr<ChMaterialSurface> material_fiber;
template<class T>
void CreateSegment(T* mSys, ChVector<> position, ChVector<> velocity) {
	real mass = 1;
	Quaternion q(1, 0, 0, 0);
	//q.Q_from_AngZ(PI / 2.0);

	string_vector[fibers] = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(string_vector[fibers], mass, Vector(0, 0, 0) + position, q, material_fiber, true, false, -1, -2);
	AddCollisionGeometry(string_vector[fibers], SPHERE, Vector(chain_radius, chain_radius, chain_radius), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(string_vector[fibers], (ChSystemParallel *) mSys);
	string_vector[fibers]->SetPos_dt(velocity);
	fibers++;
}
template<class T>
void JoinSegments(T* mSys, ChSharedBodyPtr& A, ChSharedBodyPtr& B) {

	ChCoordsys<> pos1, pos2;
	pos1.pos = Vector(-chain_radius, 0, 0);
	pos2.pos = Vector(chain_radius, 0, 0);
	ChSharedPtr<ChLinkLockRevolute> joint(new ChLinkLockRevolute);
	joint->Initialize(A, B, true, pos1, pos2);
	mSys->AddLink(joint);

}

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	if (frame * timestep > 1) {
		Roller->SetBodyFixed(false);
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_roller->Get_spe_funct())) {
			mfun->Set_yconst(2);     // rad/s  angular speed
		}
	}
}

int main(int argc, char* argv[]) {
	//omp_set_num_threads(8);

	if (argc == 6) {
		particle_cohesion = atof(argv[1]);
		particle_friction = atof(argv[2]);
		rolling_fric = atof(argv[3]);
		spinning_fric = atof(argv[4]);
		data_folder = argv[5];
	}

//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;

	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	//system_gpu->SetMaxiter(max_iter);
	//system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(50);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(20);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(20);
	system_gpu->SetTol(0);
	system_gpu->SetTolSpeeds(0);
	system_gpu->SetMaxPenetrationRecoverySpeed(1000);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(.0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(100);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->DoStabilization(false);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(40, 12, 12));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(200, 100);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	((ChSystemParallel*) system_gpu)->SetAABB(-container_size * 1.5, container_size * 1.5);
//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);
	//material->SetRollingFriction(.5);
	//material->SetSpinningFriction(.5);
	material->SetCompliance(0);
	material->SetCohesion(0);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);

	real roller_rad = 1;
	real roller_width = 2;
	real anchor_dim = .5;
	real system_height = -1;
	Vector chain_start(container_size.x - 2, system_height - roller_rad - chain_radius, 0);
	Vector chain_end(-container_size.x + 2, system_height - roller_rad - chain_radius, 0);

	material_fiber = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_fiber->SetFriction(0);
	material_fiber->SetRollingFriction(0);
	material_fiber->SetSpinningFriction(0);
	material_fiber->SetCompliance(0);
	material_fiber->SetCohesion(-1000);

	real rx = anchor_dim;
	real ry = anchor_dim;
	real rz = anchor_dim;
	real mass_anchor = 50;

	ChSharedBodyPtr Anchor = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(Anchor, mass_anchor, Vector(-container_size.x + 2 + chain_radius * 2, system_height - roller_rad - chain_radius, 0), Quaternion(1, 0, 0, 0), material_fiber, true, false, -11, -11);
	AddCollisionGeometry(Anchor, BOX, Vector(rx, ry, rz), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Anchor, SPHERE, Vector(.01, .01, .01), Vector(-rx, -ry, -rz), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Anchor, SPHERE, Vector(.01, .01, .01), Vector(-rx, -ry, rz), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Anchor, SPHERE, Vector(.01, .01, .01), Vector(rx, -ry, rz), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Anchor, SPHERE, Vector(.01, .01, .01), Vector(rx, -ry, -rz), Quaternion(1, 0, 0, 0));
	FinalizeObject(Anchor, (ChSystemParallel *) system_gpu);

	//Anchor->SetInertiaXX( ChVector<>(1 / 12.0 * mass_anchor * (ry * ry + rz * rz), 1 / 12.0 * mass_anchor * (rx * rx + rz * rz), 1 / 12.0 * mass_anchor * (rx * rx + ry * ry)));

	Quaternion roller_quat;
	roller_quat.Q_from_AngX(PI / 2.0);

	Roller = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(Roller, 1000, Vector(container_size.x - 2, system_height, 0), Quaternion(1, 0, 0, 0), material_fiber, true, false, -5, -5);
	AddCollisionGeometry(Roller, CYLINDER, Vector(roller_rad, roller_width, roller_rad), Vector(0, 0, 0), roller_quat);
	FinalizeObject(Roller, (ChSystemParallel *) system_gpu);
	//Roller->SetBodyFixed(true);

	eng_roller = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	eng_roller->Initialize(Roller, Bottom, ChCoordsys<>(ChVector<>(container_size.x - 2, system_height, 0), Q_from_AngZ(CH_C_PI / 2.0)));
	eng_roller->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);     // also works as revolute support
	eng_roller->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	//system_gpu->AddLink(eng_roller);

//	int chain_segments = (chain_start.x - chain_end.x - anchor_dim) / (chain_radius * 2);
//
//	for (int i = 0; i < chain_segments; i++) {
//		CreateSegment(system_gpu, chain_start - ChVector<>(i * chain_radius * 2, 0, 0), Vector(0, 0, 0));
//	}
//
//	for (int i = 2; i < chain_segments + 1; i++) {
//		JoinSegments(system_gpu, string_vector[i - 2], string_vector[i - 1]);
//	}
//
//	ChCoordsys<> pos1, pos2;
//	pos1.pos = Vector(0, chain_radius, 0);
//	pos2.pos = Vector(0, -roller_rad, 0);
//	//pos2.rot = roller_quat;
//	ChSharedPtr<ChLinkLockLock> joint_roller(new ChLinkLockLock);
//	joint_roller->Initialize(string_vector[0], Roller, true, pos1, pos2);
//	system_gpu->AddLink(joint_roller);
////
//	pos1.pos = Vector(-chain_radius, 0, 0);
//	pos2.pos = Vector(anchor_dim, 0, 0);
//	ChSharedPtr<ChLinkLockRevolute> joint_anchor(new ChLinkLockRevolute);
//	joint_anchor->Initialize(string_vector[fibers - 1], Anchor, true, pos1, pos2);
//	system_gpu->AddLink(joint_anchor);

	int3 num_per_dir;
	num_per_dir = I3(200, 30, 80);
	//num_per_dir = I3(50, 8, 20);
	layer_gen = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
	layer_gen->SetDensity(50);
	layer_gen->SetRadius(R3(particle_radius, particle_radius, particle_radius));
	layer_gen->material->SetFriction(particle_friction);
	layer_gen->material->SetCohesion(particle_cohesion);
	layer_gen->material->SetRollingFriction(rolling_fric);
	layer_gen->material->SetSpinningFriction(spinning_fric);
	layer_gen->material->SetCompliance(0);
	layer_gen->AddMixtureType(MIX_SPHERE);
	//layer_gen.SetNormalDistribution(rad.x, rad.x/4.0);
	//layer_gen->UseNormalCohesion(particle_cohesion, 1);

	//layer_gen->addPerturbedVolumeMixture(R3(0 , -1, 0 ), I3(num_per_dir.x, num_per_dir.y, num_per_dir.z), R3(.01, .01, .01), R3(0, 0, 0));
	//layer_gen.addPerturbedVolume(R3(2, -2, 0), SPHERE, num_per_dir, R3(.1, .1, .1), R3(-22, 0, 0), false);

//=========================================================================================================
//Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	//openGLView.render_camera->camera_pos = Vector(0, -5, -10);
	//openGLView.render_camera->look_at = Vector(0, -5, 0);
	//openGLView.render_camera->mScale = .4;
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
		double RESID = ((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

		int save_every = 1.0 / timestep / 60.0;     //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
