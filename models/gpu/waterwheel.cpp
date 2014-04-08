#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 30;

int max_iter = 10;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(3, 5, 10);
real container_thickness = .2;
real container_height = -6;
real container_friction = .1;

real particle_radius = .04;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = .001;
real particle_cohesion = .5;
real rolling_fric = 1;
real spinning_fric = 1;
Vector particle_initial_vel = Vector(0, -5.5, 0);     //initial velocity

int particle_grid_x = 2;
int particle_grid_z = 2;
real start_height = 1;

ChSharedBodyPtr impactor;
ParticleGenerator<ChSystemParallel>* water_generator;
ParticleGenerator<ChSystemParallel>* bear_inside;
real3 mass = R3(1, 1, 1);
real3 friction = R3(0, .1, 0);
real cohesion_water = 1 * timestep;
real cohesion_goop = 1 * timestep;
ChSharedBodyPtr Wheel;
string data_folder = "data/waterwheel";
template<class T>
void RunTimeStep(T* mSys, const int frame) {
	if (frame % 20 == 0) {
		water_generator->addPerturbedVolumeMixture(R3(0, .4, 9.5), I3(10, 1, 10), R3(.05, .05, 0), R3(0, -5, 0));
	}

	//Wheel->SetPos(Vector(0, -6, 4) );
	//Wheel->SetPos_dt(Vector(0, 0, 0) );
	//Wheel->SetWvel_loc(Vector(Wheel->GetWvel_loc().x, 0, 0) );

}

int main(int argc, char* argv[]) {
	//omp_set_num_threads(8);

	if (argc == 4) {

		cohesion_water = atof(argv[1]);
		cohesion_goop = atof(argv[2]);
		data_folder = argv[3];
		cout << cohesion_water << " " << cohesion_goop << " " << data_folder << endl;
	}

//=========================================================================================================
	ChSystemParallelDVI * system_gpu = new ChSystemParallelDVI;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	//system_gpu->SetMaxiter(max_iter);
	//system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(max_iter);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(max_iter);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(max_iter * 2);
	system_gpu->SetTol(.5);
	system_gpu->SetTolSpeeds(.5);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(.5);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(.0);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(1);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(100, 100, 100));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(200, 100);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	//((ChSystemParallel*) system_gpu)->SetAABB(R3(-6, -6, -6), R3(6, 6, 6));

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

	ChSharedPtr<ChMaterialSurface> material_spout;
	material_spout = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_spout->SetFriction(0);
	//material->SetRollingFriction(.5);
	//material->SetSpinningFriction(.5);
	material_spout->SetCompliance(0);
	material_spout->SetCohesion(-1000);

	ChSharedBodyPtr Spout = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(Spout, 100000, Vector(0, 0, 8), Quaternion(1, 0, 0, 0), material_spout, true, true, -20, -20);

	real width = 1;
	real thick = .1;
	real depth = 2;
	real height = .5;

	AddCollisionGeometry(Spout, BOX, Vector(width, thick, depth), Vector(0, height, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Spout, BOX, Vector(width, thick, depth), Vector(0, -height, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Spout, BOX, Vector(thick, height, depth), Vector(width, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Spout, BOX, Vector(thick, height, depth), Vector(-width, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Spout, BOX, Vector(width, height, thick), Vector(0, 0, depth), Quaternion(1, 0, 0, 0));

	FinalizeObject(Spout, (ChSystemParallel *) system_gpu);

	real wheel_rad = 3;
	real paddle_len = .5;
	real thickness = .05;
	real wheel_width = 1;
	Wheel = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(Wheel, 100, Vector(0, -6, 4), Quaternion(1, 0, 0, 0), material, true, false, -20, -20);

	AddCollisionGeometry(Wheel, CYLINDER, Vector(wheel_rad, wheel_width, wheel_rad), Vector(0, 0, 0), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Z));
	AddCollisionGeometry(Wheel, CYLINDER, Vector(wheel_rad + paddle_len * 2, thickness, wheel_rad + paddle_len * 2), Vector(-wheel_width, 0, 0), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Z));
	AddCollisionGeometry(Wheel, CYLINDER, Vector(wheel_rad + paddle_len * 2, thickness, wheel_rad + paddle_len * 2), Vector(wheel_width, 0, 0), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Z));

	for (int i = 0; i < 360; i += 30) {

		real angle = i * PI / 180.0;
		real x = cos(angle) * (wheel_rad + paddle_len);
		real z = sin(angle) * (wheel_rad + paddle_len);
		Quaternion q;
		q.Q_from_AngAxis(-angle, Vector(1, 0, 0));

		AddCollisionGeometry(Wheel, BOX, Vector(wheel_width, thickness, paddle_len), Vector(0, z, x), q);

	}

	FinalizeObject(Wheel, (ChSystemParallel *) system_gpu);

	ChSharedPtr<ChLinkLockRevolute> my_link1(new ChLinkLockRevolute);
	my_link1->Initialize(Bottom, Wheel, ChCoordsys<>(Vector(0, -6, 4), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Y)));
	system_gpu->AddLink(my_link1);

	water_generator = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
	water_generator->SetDensity(1000);
	water_generator->SetRadius(R3(particle_radius, particle_radius, particle_radius));
	water_generator->material->SetFriction(particle_friction);
	water_generator->material->SetCohesion(cohesion_water);
	water_generator->material->SetRollingFriction(0);
	water_generator->material->SetSpinningFriction(0);
	water_generator->material->SetCompliance(0);
	water_generator->AddMixtureType(MIX_ELLIPSOID);
//=========================================================================================================
////Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	//openGLView.render_camera->camera_position = glm::vec3(0, -5, -10);
//	//openGLView.render_camera->camera_look_at = glm::vec3(0, -5, 0);
//	//openGLView.render_camera->camera_scale = .4;
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
		double RESID = ((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverParallelDVI*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

		int save_every = 1.0 / timestep / 60.0;     //save data every n steps
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
