#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 60;

int max_iter = 30;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(6, 1, 6);
real container_thickness = .1;
real container_height = 0;
real container_friction = 1;

real particle_radius = .05;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = 0;
Vector particle_initial_vel = Vector(0, -5.5, 0);     //initial velocity

int particle_grid_x = 2;
int particle_grid_z = 2;
real start_height = 1;

ChSharedBodyPtr slicer1, slicer2, spinner;

real3 mass = R3(1, 1, 1);
real3 friction = R3(0, .1, 0);
real cohesion = 0;
real ang = 0;
ParticleGenerator *layer_gen;
ChSharedPtr<ChMaterialSurface> material_fiber;
template<class T>
void CreateFiber(T* mSys, ChVector<> position) {



	real length = .1;
	real thickness = .05;

	ChSharedBodyPtr Fiber1 = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

	ChSharedBodyPtr Fiber2;
	vector<ChSharedBodyPtr> fibers(20);
	fibers[0] = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	Quaternion q;
	q.Q_from_AngZ(PI / 2.0);
	InitObject(fibers[0], 1, Vector(0, -2, 0)+position, q, material_fiber, true, false, -1, -2);
	AddCollisionGeometry(fibers[0], CYLINDER, Vector(thickness, length, length), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(fibers[0], (ChSystemGPU *) mSys);
	fibers[0]->SetPos_dt(Vector(0, -1, 0));

	for (int i = 1; i < 20; i++) {
		fibers[i] = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
		InitObject(fibers[i], 1, Vector(i * length * 2.1, -2, 0)+position, q, material_fiber, true, false, -1, -2);
		AddCollisionGeometry(fibers[i], CYLINDER, Vector(thickness, length, length), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
		FinalizeObject(fibers[i], (ChSystemGPU *) mSys);
		fibers[i]->SetPos_dt(Vector(0, -1, 0));
		ChCoordsys<> pos;
		pos.pos = Vector((i - 1) * length * 2.1 + length, -2, 0)+position;
		ChSharedPtr<ChLinkLockRevolute> joint(new ChLinkLockRevolute);
		joint->Initialize(fibers[i - 1], fibers[i], pos);
		//		joint->Initialize(fibers[i - 1], fibers[i], true, Vector(length,0,0),Vector(-length,0,0),true);
		//		joint->Set_SpringF(100);
		//		joint->Set_SpringK(100);
		//		joint->Set_SpringR(100);

		mSys->AddLink(joint);

	}

}

template<class T>
void RunTimeStep(T* mSys, const int frame) {

	if (frame % 100 == 0) {
		for(int i=0; i<20; i++)
		CreateFiber(mSys,Vector(0,5,i/5.0));

	}
}

int main(int argc, char* argv[]) {
	omp_set_num_threads(8);
//=========================================================================================================
	ChSystemGPU * system_gpu = new ChSystemGPU;
	ChCollisionSystemGPU *mcollisionengine = new ChCollisionSystemGPU();
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(.1);
	system_gpu->SetTolSpeeds(.1);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(.1);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(50);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .01);
	mcollisionengine->setBinsPerAxis(R3(50, 50, 50));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr Tube = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(.4);
	material->SetRollingFriction(.4);
	material->SetSpinningFriction(.4);
	material->SetCompliance(0);
	material->SetCohesion(-100);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Top, 100000, Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Tube, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	AddCollisionGeometry(Tube, BOX, Vector(2, container_thickness / 6.0, 1), Vector(0, container_size.y / 2.0 + .6 + .4, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Tube, BOX, Vector(2, container_thickness / 6.0, 1), Vector(0, container_size.y / 2.0 - .6 + .4, 0), Quaternion(1, 0, 0, 0));

	AddCollisionGeometry(Tube, BOX, Vector(2, .6, container_thickness / 6.0), Vector(0, container_size.y / 2.0 + .4, -1), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Tube, BOX, Vector(2, .6, container_thickness / 6.0), Vector(0, container_size.y / 2.0 + .4, 1), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);
	//FinalizeObject(Top, (ChSystemGPU *) system_gpu);
	//FinalizeObject(Tube, (ChSystemGPU *) system_gpu);

	material_fiber = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_fiber->SetFriction(.4);
	material_fiber->SetRollingFriction(.4);
	material_fiber->SetSpinningFriction(.4);
	material_fiber->SetCompliance(0);
	material_fiber->SetCohesion(0);


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
