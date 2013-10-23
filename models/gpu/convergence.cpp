#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"

ChVector<> lpos(0, 0, 0);
ChQuaternion<> quat(1, 0, 0, 0);

//all dimensions are in millimeters, milligrams
real plate_height = 0;
real plate_thickness = 1;
real plate_radius = 7;
real plate_friction = 1;

real particle_radius = .2;
real particle_mass = .5;
real particle_density = 1.14;
real particle_friction = 1.0;
Vector particle_initial_vel = Vector(0, -5, 0);     //initial velocity

real container_width = 3.0;
real container_thickness = .25;
real container_height = 6.0;
real wscale = 1;

int max_iter = 150;
real tolerance = 0;

real gravity = -9.810;
real timestep = .002;
real seconds_to_simulate = 3;
int num_steps = seconds_to_simulate / timestep;

int particle_grid_x = 10;
int particle_grid_z = 10;
int particles_every = 50;     //add particles every n steps
int save_every = 1.0 / timestep / 100.0;     //save data every n steps
ChSharedBodyPtr BLOCK;
int particle_configuration = 0;
//0: single sphere
//1: two spheres joined together
template<class T>
void RunTimeStep(T* mSys, const int frame) {

	BLOCK->SetRot(ChQuaternion<>(1, 0, 0, 0));
	BLOCK->SetWvel_loc(ChVector<>(0, 0, 0));
	BLOCK->SetPos(ChVector<>(0, BLOCK->GetPos().y, 0));
	BLOCK->SetPos_dt(ChVector<>(0, BLOCK->GetPos_dt().y, 0));
	if (frame % particles_every == 0) {
		//addHCPSheet(10, 10, 0, particle_mass, particle_radius, particle_friction, true, 0, 0, particle_initial_vel, (ChSystemParallel*) mSys);
	}


	real residual = ((ChLcpSolverParallel *) (mSys->GetLcpSolverSpeed()))->GetResidual()*.5;
	mSys->SetTol(residual);
	mSys->SetTolSpeeds(residual);
	((ChLcpSolverParallel *) (mSys->GetLcpSolverSpeed()))->SetTolerance(residual);


//	if(frame >300){
//
//		BLOCK->SetPos(ChVector<>(0, 5, 0));
//
//	}

}
int main(int argc, char* argv[]) {
	omp_set_num_threads(4);

	//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	ChCollisionSystemParallel *mcollisionengine = new ChCollisionSystemParallel();
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

	//=========================================================================================================
	//ReadInputFile("convergence_config.txt",system_gpu);
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(tolerance);
	system_gpu->SetTolSpeeds(tolerance);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(tolerance);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(10);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->DoStabilization(false);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .01);
	mcollisionengine->setBinsPerAxis(R3(30, 30, 30));
	mcollisionengine->setBodyPerBin(100, 50);
	//=========================================================================================================
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	//=========================================================================================================
//	Quaternion plate_quat;
//	plate_quat.Q_from_AngAxis(0, Vector(1, 0, 0));
//
//	ChSharedBodyGPUPtr PLATE = ChSharedBodyGPUPtr(new ChBody(new ChCollisionModelParallel));
//	InitObject(PLATE, 1, ChVector<>(0, plate_height, 0), plate_quat, plate_friction, plate_friction, 0, true, true, -1000, -20000);
//	AddCollisionGeometry(PLATE, BOX, ChVector<>(plate_radius, plate_thickness, plate_radius), lpos, quat);
//	FinalizeObject(PLATE, (ChSystemParallel *) system_gpu);

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr BTM = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	BLOCK = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(plate_friction);
	material->SetCompliance(0);

	InitObject(L, 100000, Vector(-container_width + container_thickness, plate_height, 0), quat, material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_width - container_thickness, plate_height, 0), quat, material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, plate_height, -container_width + container_thickness), quat, material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, plate_height, container_width - container_thickness), quat, material, true, true, -20, -20);
	InitObject(BTM, 100000, Vector(0, plate_height - container_height, 0), quat, material, true, true, -20, -20);
	InitObject(BLOCK, 1000, Vector(0, 5, 0), quat, material, true, false, -1, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
	AddCollisionGeometry(F, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
	AddCollisionGeometry(B, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
	AddCollisionGeometry(BTM, BOX, Vector(container_width * wscale, container_thickness, container_width), lpos, quat);
	AddCollisionGeometry(BLOCK, BOX, Vector(container_width, container_thickness, container_width), lpos, quat);

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(BTM, (ChSystemParallel *) system_gpu);
	FinalizeObject(BLOCK, (ChSystemParallel *) system_gpu);

	ParticleGenerator layer_gen(system_gpu);
	layer_gen.SetMass(.1);
	layer_gen.SetRadius(R3(particle_radius ));

	layer_gen.material->SetFriction(1);
	layer_gen.material->SetRollingFriction(1);
	layer_gen.material->SetSpinningFriction(1);
	layer_gen.material->SetCohesion(0);
	layer_gen.material->SetCompliance(0);
	int3 num_per_dir = I3(10 , 20 , 10 );
	//layer_gen.addPerturbedVolume(R3(0 + container_pos.x, container_pos.y, 0 + container_pos.z), SPHERE, I3(num_per_dir.x, 1, num_per_dir.z), R3(.1, .1, .1), R3(0, 0, 0), false);
	//layer_gen.SetNormalDistribution(particle_radius - particle_radius / 6.0, particle_radius / 6.0);
	layer_gen.addPerturbedVolume(R3(0, -1.5, 0), SPHERE, num_per_dir, R3(.1, .1, .1), R3(0, -4, 0), false);

	//addHCPCube(23, 1, 23, particle_mass, particle_radius, particle_friction, true, 0, -15, 0, Vector(0, 0, 0), system_gpu);

	//=========================================================================================================
	//////Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_pos = Vector(0, -5, -40);
	openGLView.render_camera->look_at = Vector(0, -5, 0);
	openGLView.SetCustomCallback(RunTimeStep);
	openGLView.StartSpinning(window_manager);
	window_manager->CallGlutMainLoop();
	//=========================================================================================================
	//timestep
	//iterations
	//warm start
	//mass
	ofstream ofile("convergence.txt");
	ChTimer<double> timer;
	timer.start();
	int file = 0;
	for (int i = 0; i < num_steps; i++) {

		cout << "step " << i << endl;
		system_gpu->DoStepDynamics(timestep);
		RunTimeStep(system_gpu, i);



		if (i % save_every == 0) {
		ofile<<((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->GetResidual()<<endl;
//			stringstream ss;
//			cout << "Frame: " << file << endl;
//			DumpObjects(system_gpu, ss.str());
//			//output.ExportData(ss.str());
//			file++;
		}
		timer.stop();

	}
	cout << "TIME: " << timer() << endl;
	ofile.close();
	return 0;
}

