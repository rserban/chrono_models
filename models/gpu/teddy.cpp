#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 5;

int max_iter = 10;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(10, 2, 10);
real container_thickness = .2;
real container_height = -4;
real container_friction = .1;

real particle_radius = .015;
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

ChSharedBodyPtr impactor;
ParticleGenerator<ChSystemParallel>* bear_outside;
ParticleGenerator<ChSystemParallel>* bear_inside;
real3 mass = R3(1, 1, 1);
real3 friction = R3(0, .1, 0);
real cohesion_bear = 10000*timestep;
real cohesion_goop = 10*timestep;
ChSharedBodyPtr Bunny;
string data_folder = "data/teddy";
template<class T>
void RunTimeStep(T* mSys, const int frame) {
}

int main(int argc, char* argv[]) {
	//omp_set_num_threads(8);

	if (argc == 4) {

		cohesion_bear = atof(argv[1]);
		cohesion_goop = atof(argv[2]);
		data_folder=argv[3];
		cout<<cohesion_bear<<" "<<cohesion_goop<<" "<<data_folder<<endl;
	}

//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	//system_gpu->SetMaxiter(max_iter);
	//system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(max_iter*2);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	system_gpu->SetTol(.5);
	system_gpu->SetTolSpeeds(.5);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(.5);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(.0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(1);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(100, 100, 100));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(200, 100);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	//((ChSystemParallel*) system_gpu)->SetAABB(R3(-6, -6, -6), R3(6, 6, 6));
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================
//
//	if (stream) {
//		Bunny = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
//		//
//		InitObject(Bunny, 1, Vector(0, -3, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
//		//	//AddCollisionGeometry(Bunny, BOX, Vector(1, 1, 1),  Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
//		AddCollisionGeometryTriangleMesh(Bunny, "bunny_low.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
//		FinalizeObject(Bunny, (ChSystemParallel *) system_gpu);
//
//		//	real mass = 1;
//		//	Vector r = ChVector<>(1, 1, 1);
//		//	Vector inertia = Vector((1 / 5.0 * mass * (r.y * r.y + r.z * r.z), 1 / 5.0 * mass * (r.x * r.x + r.z * r.z), 1 / 5.0 * mass * (r.x * r.x + r.y * r.y)));
//
//		//Bunny->SetInertiaXX(inertia);
//
//	}
//
//	if (!stream) {

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	ChSharedBodyPtr Float = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

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
	InitObject(Top, 100000, Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	Quaternion qqq;
	qqq.Q_from_AngZ(15);

	InitObject(Float, 100000, Vector(0, -3, 0), qqq, material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	//AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	//AddCollisionGeometry(Float, BOX, Vector(1, 1, 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);

	//FinalizeObject(Float, (ChSystemParallel *) system_gpu);

	Quaternion quat;
	quat.Q_from_AngAxis(90, Vector(0, 1, 0));



	real3 rad = R3(particle_radius, particle_radius, particle_radius);

	bear_outside = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
	bear_outside->SetDensity(1000);
	bear_outside->SetRadius(rad);
	bear_outside->material->SetFriction(particle_friction);
	bear_outside->material->SetCohesion(cohesion_bear);
	bear_outside->material->SetRollingFriction(0);
	bear_outside->material->SetSpinningFriction(0);
	bear_outside->material->SetCompliance(0);
	//bear_outside->SetActive(true);

	bear_outside->loadAscii("bear_015.txt", R3(0, -2, 0), SPHERE, R3(.015*(2), 0, 0), R3(0, 0, 0), R3(4,4,4));

	rad = R3(.02, .02, .02);
	bear_inside = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
	bear_inside->SetDensity(1000);
	bear_inside->SetRadius(rad);
	bear_inside->material->SetFriction(.5);
	bear_inside->material->SetCohesion(cohesion_goop);
	bear_inside->material->SetRollingFriction(0);
	bear_inside->material->SetSpinningFriction(0);
	bear_inside->material->SetCompliance(0);
	bear_inside->AddMixtureType(MIX_SPHERE);
	bear_inside->SetCylinderRadius(.015*(2)*5);
	//bear_outside->SetActive(true);
	bear_inside->loadAscii("bear_inside_015.txt", R3(0, -2, 0), SPHERE, R3(.015*(2), 0, 0), R3(0, 0, 0), R3(4,4,4));
	//bear_inside->loadAscii("bear_in_trim_.02.txt", R3(0, 0, 0), SPHERE, R3(.035, 0, 0), R3(0, 0, 0), R3(1,1,1));

//=========================================================================================================
//Rendering specific stuff:
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

		int save_every = 1.0 / timestep / 600.0;     //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}


	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
