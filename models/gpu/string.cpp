//#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .0004;
real seconds_to_simulate = 30;

int max_iter = 40;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(2, 2, 2);
real container_thickness = .1;
real container_height = 0;
real container_friction = 1;

real particle_radius = .025;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = 0;
Vector particle_initial_vel = Vector(0, -5.5, 0);     //initial velocity

int particle_grid_x = 2;
int particle_grid_z = 2;
real start_height = 1;

real3 mass = R3(1, 1, 1);
real3 friction = R3(0, .1, 0);
real cohesion = 0;
real ang = 0;

real segment_length = .05;
real segment_thickness = .02;
vector<ChSharedBodyPtr> string_vector(1000);

int fibers = 0;
ChSharedPtr<ChMaterialSurface> material_fiber;

template<class T>
void CreateSegment(T* mSys, ChVector<> position, ChVector<> velocity) {
	real mass = .1;
	Quaternion q;
	q.Q_from_AngZ(PI / 2.0);

	string_vector[fibers] = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(string_vector[fibers], mass, Vector(0, 0, 0) + position, q, material_fiber, true, false, -1, -2);
	AddCollisionGeometry(string_vector[fibers], CYLINDER, Vector(segment_thickness, segment_length, segment_length), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(string_vector[fibers], SPHERE, Vector(segment_thickness, segment_length, segment_length), Vector(0, segment_length, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(string_vector[fibers], SPHERE, Vector(segment_thickness, segment_length, segment_length), Vector(0, -segment_length, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(string_vector[fibers], (ChSystemParallel *) mSys);
	string_vector[fibers]->SetPos_dt(velocity);
	fibers++;
}
template<class T>
void JoinSegments(T* mSys, ChSharedBodyPtr& A, ChSharedBodyPtr& B) {

	real mass = .1;
	Quaternion q;
	q.Q_from_AngZ(PI / 2.0);

	ChCoordsys<> pos1;
	ChCoordsys<> pos2;
	pos1.pos = Vector(0,segment_length+segment_thickness, 0);
	pos2.pos = Vector(0,-segment_length-segment_thickness, 0);

	ChCoordsys<> pos;
	pos.pos = Vector(0, 0, 0);
	ChSharedPtr<ChLinkLockRevolute> joint(new ChLinkLockRevolute);
	joint->Initialize(A, B, true,pos1,pos2);

	mSys->AddLink(joint);

}

template<class T>
void RunTimeStep(T* mSys, const int frame) {

	int inter_frame_time = 100;
	real len = segment_length * 2 + segment_thickness * 2;
	real half_len = len * .5;
	real total_time = inter_frame_time * timestep;
	ChSharedBodyPtr A, B;
	int num_seg = 0;

	if (frame == 0) {
		CreateSegment(mSys, Vector(half_len, 0, 0), Vector(half_len / total_time * 2, 0, 0));

//		JoinSegments(mSys, string_vector[num_seg - 2], string_vector[num_seg - 1]);
//		CreateSegment(mSys, Vector(-half_len*3, 0, 0), Vector(half_len / total_time * 2, 0, 0));
//		num_seg += 1;
//		JoinSegments(mSys, string_vector[num_seg - 2], string_vector[num_seg - 1]);
	}

	if (frame == 0) {
		//CreateSegment(mSys, Vector(-half_len, 0, 0), Vector(half_len / total_time * 2, 0, 0));

		//JoinSegments(mSys, string_vector[fibers - 2], string_vector[fibers - 1]);
	}
if(fibers>5){

	string_vector[fibers - 1]->SetPos_dt(ChVector<>(half_len / total_time * 2,string_vector[fibers - 1]->GetPos_dt().y,string_vector[fibers - 1]->GetPos_dt().z));
	string_vector[fibers - 2]->SetPos_dt(ChVector<>(half_len / total_time * 2,string_vector[fibers - 2]->GetPos_dt().y,string_vector[fibers - 2]->GetPos_dt().z));
	string_vector[fibers - 3]->SetPos_dt(ChVector<>(half_len / total_time * 2,string_vector[fibers - 3]->GetPos_dt().y,string_vector[fibers - 3]->GetPos_dt().z));
	string_vector[fibers - 4]->SetPos_dt(ChVector<>(half_len / total_time * 2,string_vector[fibers - 4]->GetPos_dt().y,string_vector[fibers - 4]->GetPos_dt().z));

}
	if (frame % inter_frame_time == 0 && frame * timestep < 20) {
		CreateSegment(mSys, Vector(-half_len, 0, 0), Vector(half_len / total_time * 2, 0, 0));
		//num_seg += 1;
		JoinSegments(mSys, string_vector[fibers - 2], string_vector[fibers - 1]);
////		for (int i = 0; i < 40; i++) {
////			CreateFiber(mSys, Vector(2 - 1.6 * 0, 3, i / 8.0));
////			CreateFiber(mSys, Vector(2 - 1.6 * 1, 3, i / 8.0));
////			CreateFiber(mSys, Vector(2 - 1.6 * 2, 3, i / 8.0));
////			CreateFiber(mSys, Vector(2 - 1.6 * 3, 3, i / 8.0));
////			fibers += 4;
////		}
////		cout << "Fibers: " << fibers << endl;
	}
}

int main(int argc, char* argv[]) {
	int threads = 8;

	if (argc > 1) {
		threads = (atoi(argv[1]));
	}

	omp_set_num_threads(threads);
//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	ChCollisionSystemParallel *mcollisionengine = new ChCollisionSystemParallel();
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	system_gpu->SetParallelThreadNumber(threads);
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(10);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(10);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(20);
	system_gpu->SetTol(.1);
	system_gpu->SetTolSpeeds(.1);
	system_gpu->SetMaxPenetrationRecoverySpeed(100);

	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(.1);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(10);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .01);
	mcollisionengine->setBinsPerAxis(I3(50, 50, 50));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
	((ChSystemParallel*) system_gpu)->SetAABB(R3(-6, -3, -12), R3(6, 6, 12));
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Tube = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(.1);
	material->SetRollingFriction(0);
	material->SetSpinningFriction(0);
	material->SetCompliance(0);
	material->SetCohesion(-100);

	Quaternion q;
	q.Q_from_AngX(-.1);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Top, 100000, Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	ChSharedPtr<ChMaterialSurface> material_tube;
	material_tube = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_tube->SetFriction(0);
	material_tube->SetRollingFriction(0);
	material_tube->SetSpinningFriction(0);
	material_tube->SetCompliance(0);
	material_tube->SetCohesion(-100);

	InitObject(Tube, 100000, Vector(0, 0, 0), Quaternion(1, 0, 0, 0), material_tube, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	real tube_size = .05;
	AddCollisionGeometry(Tube, BOX, Vector(tube_size * 10, container_thickness / 4.0, tube_size), Vector(0, tube_size, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Tube, BOX, Vector(tube_size * 10, container_thickness / 4.0, tube_size), Vector(0, -tube_size, 0), Quaternion(1, 0, 0, 0));

	AddCollisionGeometry(Tube, BOX, Vector(tube_size * 10, tube_size, container_thickness / 4.0), Vector(0, 0, -tube_size), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Tube, BOX, Vector(tube_size * 10, tube_size, container_thickness / 4.0), Vector(0, 0, tube_size), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);
//FinalizeObject(Top, (ChSystemParallel *) system_gpu);
	FinalizeObject(Tube, (ChSystemParallel *) system_gpu);

	material_fiber = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_fiber->SetFriction(.4);
	material_fiber->SetRollingFriction(1);
	material_fiber->SetSpinningFriction(1);
	material_fiber->SetCompliance(0);
	material_fiber->SetCohesion(0);

//=========================================================================================================
//Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
////openGLView.render_camera->camera_pos = Vector(0, -5, -10);
////openGLView.render_camera->look_at = Vector(0, -5, 0);
//	openGLView.render_camera->mScale = .1;
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
			ss << "data/fiber/" << "/" << file << ".txt";
			//DumpAllObjects(system_gpu, ss.str(), ",", true);
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
