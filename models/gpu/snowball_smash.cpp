#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .0005;
real seconds_to_simulate = 5;

int max_iter = 10;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(12, 3, 6);
real container_thickness = .2;
real container_height = 0;
real container_friction = .1;

real particle_radius = .02;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = 1;
real particle_cohesion = .5;
real rolling_fric = 1;
real spinning_fric = 1;
Vector particle_initial_vel = Vector(0, -5.5, 0);     //initial velocity

int particle_grid_x = 2;
int particle_grid_z = 2;
real start_height = 1;

ParticleGenerator<ChSystemParallel>* layer_gen;

real3 mass = R3(1, 1, 1);
real3 friction = R3(0, .1, 0);
real cohesion = .1;
ChSharedBodyPtr Bunny;
string data_folder = "data/snowball_smash";
template<class T>
void RunTimeStep(T* mSys, const int frame) {
}

int main(int argc, char* argv[]) {
	//omp_set_num_threads(8);

	if (argc == 6) {
		cohesion = atof(argv[1]);
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
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(max_iter * 2);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	system_gpu->SetTol(1e-4);
	system_gpu->SetTolSpeeds(1e-4);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(1e-4);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, .0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(25);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(120, 30, 60));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(200, 100);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	((ChSystemParallel*) system_gpu)->SetAABB(R3(-12, -3, -6), R3(12, 3, 6));
//=========================================================================================================


	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));


	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);
	//material->SetRollingFriction(.5);
	//material->SetSpinningFriction(.5);
	material->SetCompliance(0);
	material->SetCohesion(-1000);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Top, 100000, Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	Quaternion qqq;
	qqq.Q_from_AngZ(15);


	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));


	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);


	Quaternion quat;
	quat.Q_from_AngAxis(90, Vector(0, 1, 0));

	ChSharedPtr<ChMaterialSurface> material_monkey;
	material_monkey = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_monkey->SetFriction(container_friction);
	material_monkey->SetCompliance(0);
	material_monkey->SetCohesion(10);


	real3 rad = R3(particle_radius, particle_radius, particle_radius);

	layer_gen = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
	layer_gen->SetDensity(1000);
	layer_gen->SetRadius(rad);
	layer_gen->material->SetFriction(particle_friction);
	layer_gen->material->SetCohesion(cohesion);
	layer_gen->material->SetRollingFriction(rolling_fric);
	layer_gen->material->SetSpinningFriction(spinning_fric);
	layer_gen->material->SetCompliance(0);

	layer_gen->addSnowball(R3(-2, 0, 0), SPHERE, 1, R3(10, 0, 0));
	layer_gen->addSnowball(R3(2, 0, 0), SPHERE, 1, R3(-10, 0, 0));

//=========================================================================================================
////Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	//openGLView.render_camera->camera_pos = Vector(0, -5, -10);
//	//openGLView.render_camera->look_at = Vector(0, -5, 0);
//	//openGLView.render_camera->mScale = .4;
//	openGLView.SetCustomCallback(RunTimeStep);
//	openGLView.StartSpinning(window_manager);
//	window_manager->CallGlutMainLoop();
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
			DumpAllObjects(system_gpu, ss.str(), ",", true);
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
