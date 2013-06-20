#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 20;

int max_iter = 10;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(4, 4, 4);
real container_thickness = .2;
real container_height = 0;
real container_friction = 1;

real particle_radius = .025;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = 0;
Vector particle_initial_vel = Vector(0, -5.5, 0); //initial velocity

int particle_grid_x = 2;
int particle_grid_z = 2;
real start_height = 1;

ChSharedBodyGPUPtr axle_F;
ChSharedBodyGPUPtr axle_R;

real3 mass = R3(1, 1, 1);
real3 friction = R3(0, .1, 0);
real cohesion = 0;
real ang = 0;

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	if (mSys->GetNbodies() < 1071630) {
		ChSharedBodyGPUPtr sphere;
		real3 rad = R3(particle_radius, particle_radius, particle_radius);
		real3 size = container_size;
		size.y = container_size.y / 3.0;

		int3 num_per_dir = I3(10, 1, 20);

		if (frame % 10 == 0) {
			//addPerturbedLayer(R3(-2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.x, friction.x, cohesion.x, R3(0, 5, 0), (ChSystemGPU*) mSys);
			addPerturbedLayer(R3(0, 3, 0), SPHERE, rad, num_per_dir, R3(0, 0, 0), mass.y, 1, .2, R3(0, -4, 0), (ChSystemGPU*) mSys, 0);
			//addPerturbedLayer(R3(2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.z, friction.z, cohesion.z, R3(0, 5, 0), (ChSystemGPU*) mSys);
		}
	}

	axle_F->SetPos(Vector(0, 0, .6));
	axle_R->SetPos(Vector(0, 0, -.6));
	axle_F->SetPos_dt(Vector(0, 0, 0));
	axle_R->SetPos_dt(Vector(0, 0, 0));

	ang += CH_C_PI * timestep * 8;
	if (ang >= 2 * CH_C_PI) {
		ang = 0;
	}

	Quaternion q1;
	q1.Q_from_AngY(ang);
	Quaternion q2;
	q1 = Q_from_AngX(-ang);
	q2 = Q_from_AngX(ang); //*Q_from_AngZ(CH_C_PI / 2.0);

	axle_F->SetRot(q1);
	axle_R->SetRot(q2);

	axle_F->SetWvel_loc(Vector(-16, 0, 0));
	axle_R->SetWvel_loc(Vector(16, 0, 0));

}

int main(int argc, char* argv[]) {
	omp_set_num_threads(8);
//=========================================================================================================
	ChSystemGPU * system_gpu = new ChSystemGPU;
	ChLcpSystemDescriptorGPU *mdescriptor = new ChLcpSystemDescriptorGPU();
	ChContactContainerGPU *mcontactcontainer = new ChContactContainerGPU();
//ChCollisionSystemBulletGPU *mcollisionengine = new ChCollisionSystemBulletGPU();
	ChCollisionSystemGPU *mcollisionengine = new ChCollisionSystemGPU();
	system_gpu->ChangeLcpSystemDescriptor(mdescriptor);
	system_gpu->ChangeContactContainer(mcontactcontainer);
	system_gpu->ChangeCollisionSystem(mcollisionengine);
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(0);
	system_gpu->SetTolSpeeds(0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(1);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	mcollisionengine->setBinsPerAxis(R3(30, 30, 30));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedBodyGPUPtr L = ChSharedBodyGPUPtr(new ChBodyGPU);
	ChSharedBodyGPUPtr R = ChSharedBodyGPUPtr(new ChBodyGPU);
	ChSharedBodyGPUPtr F = ChSharedBodyGPUPtr(new ChBodyGPU);
	ChSharedBodyGPUPtr B = ChSharedBodyGPUPtr(new ChBodyGPU);
	ChSharedBodyGPUPtr Bottom = ChSharedBodyGPUPtr(new ChBodyGPU);
	ChSharedBodyGPUPtr Top = ChSharedBodyGPUPtr(new ChBodyGPU);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(Top, 100000, Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	L->SetCohesion(5);
	R->SetCohesion(5);
	F->SetCohesion(5);
	B->SetCohesion(5);
	Bottom->SetCohesion(5);
	Top->SetCohesion(5);

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);
	FinalizeObject(Top, (ChSystemGPU *) system_gpu);

	ChSharedBodyGPUPtr chassis(new ChBodyGPU);
	axle_F = ChSharedBodyGPUPtr(new ChBodyGPU);
	axle_R = ChSharedBodyGPUPtr(new ChBodyGPU);

	real chassisL = .8;
	InitObject(chassis, 1.0, ChVector<>(0, 0, 0), Quaternion(1, 0, 0, 0), 0, 0, 0, false, true, 0, 1);
	InitObject(axle_F, 100, ChVector<>(0, 0, chassisL / 2), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, -2, -2);
	InitObject(axle_R, 100, ChVector<>(0, 0, -chassisL / 2), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, -2, -2);

	AddCollisionGeometry(axle_F, ELLIPSOID, ChVector<>(4 / 2.0, 1 / 2.0, 1 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(axle_R, ELLIPSOID, ChVector<>(4 / 2.0, 1 / 2.0, 1 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	axle_F->SetCohesion(-.1);
	axle_R->SetCohesion(-.1);
	FinalizeObject(chassis, (ChSystemGPU *) system_gpu);
	FinalizeObject(axle_F, (ChSystemGPU *) system_gpu);
	FinalizeObject(axle_R, (ChSystemGPU *) system_gpu);

	ChSharedBodyPtr chassis_ptr = ChSharedBodyPtr(chassis);
	ChSharedBodyPtr axle_F_ptr = ChSharedBodyPtr(axle_F);
	ChSharedBodyPtr axle_R_ptr = ChSharedBodyPtr(axle_R);

	ChSharedBodyGPUPtr wall(new ChBodyGPU);
	InitObject(wall, 1.0, ChVector<>(0, 1, 0), Quaternion(1, 0, 0, 0), 0, 0, 0, true, true, 0, 1);
	AddCollisionGeometry(wall, BOX, ChVector<>(.1, 3, 1), Vector(-.75, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(wall, BOX, ChVector<>(.1, 3, 1), Vector(.75, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(wall, BOX, ChVector<>(1, 3, .1), Vector(0, 0, -.75), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(wall, BOX, ChVector<>(1, 3, .1), Vector(0, 0, .75), Quaternion(1, 0, 0, 0));

	FinalizeObject(wall, (ChSystemGPU *) system_gpu);


//=========================================================================================================
//Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
//	openGLView.render_camera->look_at = Vector(0, -5, 0);
//	openGLView.render_camera->mScale = .5;
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
		double RESID = ((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverGPU*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

		int save_every = 1.0 / timestep / 60.0; //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << "data/rollers/" << "/" << file << ".txt";
			//DumpAllObjects(system_gpu, ss.str(), ",", true);
			DumpAllObjectsWithGeometry(system_gpu, ss.str(), ",");
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
