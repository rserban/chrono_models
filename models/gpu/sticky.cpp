#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .0005;
real seconds_to_simulate = 30;

int max_iter = 10;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(6, 6, 6);
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

ChSharedBodyPtr impactor;
ParticleGenerator* layer_gen;

real3 mass = R3(1, 1, 1);
real3 friction = R3(0, .1, 0);
real cohesion = .25;
ChSharedBodyPtr Bunny;
string data_folder = "data/sticky";
template<class T>
void RunTimeStep(T* mSys, const int frame) {

	if (frame % 1000 == 0) {
		layer_gen->addSnowball(R3(0, 2, 0), SPHERE, 1, R3(0, -10, 0));
	}
}

int main(int argc, char* argv[]) {
	//omp_set_num_threads(8);

	if (argc == 6) {
		cohesion = atof(argv[1]);
		particle_friction = atof(argv[2]);
		rolling_fric = atof(argv[3]);
		spinning_fric = atof(argv[4]);
		data_folder=argv[5];
	}

//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	ChCollisionSystemParallel *mcollisionengine = new ChCollisionSystemParallel();
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	//system_gpu->SetMaxiter(max_iter);
	//system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(max_iter*2);
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
	mcollisionengine->setBinsPerAxis(I3(100, 100, 100));
	mcollisionengine->setBodyPerBin(200, 100);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	((ChSystemParallel*) system_gpu)->SetAABB(R3(-6, -6, -6), R3(6, 6, 6));
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
	material->SetCohesion(-1000);

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

	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	AddCollisionGeometry(Float, BOX, Vector(1, 1, 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);

	FinalizeObject(Float, (ChSystemParallel *) system_gpu);

	Quaternion quat;
	quat.Q_from_AngAxis(90, Vector(0, 1, 0));

	ChSharedPtr<ChMaterialSurface> material_monkey;
	material_monkey = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_monkey->SetFriction(container_friction);
	material_monkey->SetCompliance(0);
	material_monkey->SetCohesion(10);

	//ChSharedBodyPtr Monkey = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	//InitObject(Monkey, 100000, Vector(0,-2,0), quat, material_monkey, true, true, -20, -20);
	//AddCollisionGeometryTriangleMesh(Monkey, "monkey.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//FinalizeObject(Monkey, (ChSystemParallel *) system_gpu);

//		ifstream ifile("bunny_init.txt");
//		string temp;
//
//		for (int i = 0; i < 27000; i++) {
//			if (i > 0) {
//				//cout<<i<<endl;
//				getline(ifile, temp);
//				stringstream ss(temp);
//				Vector pos;
//				Vector vel;
//				ss >> pos.x >> pos.y >> pos.z >> vel.x >> vel.y >> vel.z;
//				ChSharedBodyGPUPtr sphere = ChSharedBodyGPUPtr(new ChBody(new ChCollisionModelGPU));
//				InitObject(sphere, 1, pos, Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, 1, 0);
//				AddCollisionGeometry(sphere, SPHERE, ChVector<>(particle_radius, particle_radius, particle_radius), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
//				FinalizeObject(sphere, (ChSystemParallel *) system_gpu);
//				sphere->SetCohesion(cohesion);
//				sphere->SetPos_dt(vel);
//			}
//		}

	real3 rad = R3(particle_radius, particle_radius, particle_radius);
//		real3 size = container_size;
//		size.y = container_size.y / 3.0;
//
	int3 num_per_dir;
//		cout << num_per_dir.x * num_per_dir.y * num_per_dir.z * 3 << endl;
	//num_per_dir = I3(1, size.y / rad.y * .85, 1);
//
	num_per_dir = I3(40, 40, 40);

	layer_gen = new ParticleGenerator((ChSystemParallel *) system_gpu);
	layer_gen->SetDensity(1000);
	layer_gen->SetRadius(rad);
	layer_gen->material->SetFriction(particle_friction);
	layer_gen->material->SetCohesion(cohesion);
	layer_gen->material->SetRollingFriction(rolling_fric);
	layer_gen->material->SetSpinningFriction(spinning_fric);
	layer_gen->material->SetCompliance(0);
	//layer_gen.SetNormalDistribution(rad.x, rad.x/4.0);
	//layer_gen->UseNormalCohesion(particle_cohesion, 1);

	//layer_gen.addPerturbedVolume(R3(0 , 0, 0 ), SPHERE, I3(num_per_dir.x, 1, num_per_dir.z), R3(.1, .1, .1), R3(0, 0, 0), false);
	//layer_gen.addPerturbedVolume(R3(2, -2, 0), SPHERE, num_per_dir, R3(.1, .1, .1), R3(-22, 0, 0), false);
	//addPerturbedLayer(R3(-2, -2, 0), SPHERE, rad, num_per_dir, R3(.1,.1,.1), .333, .01, rand()%5, R3(44, 0, 0), system_gpu);
//		addPerturbedLayer(R3(0, 0, 0), SPHERE, rad, num_per_dir, R3(.1, .1, .1), .666, 0, 0, R3(0, 0, 0), system_gpu);
//		addPerturbedLayer(R3(0, 2, 0), SPHERE, rad, num_per_dir, R3(.1, .1, .1), .999, 0, 0, R3(0, 0, 0), system_gpu);
//	}

	//layer_gen->loadAscii("bunny.txt", R3(0, 2, 0), SPHERE, R3(.0235, 0, 0), R3(0, -10, 0), R3(2, 2, 2));

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
