#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
ChVector<> lpos(0, 0, 0);
ChQuaternion<> quat(1, 0, 0, 0);

//all dimensions are in millimeters, milligrams
real plate_height = 0;
real plate_thickness = 1;
real plate_radius = 7;
real plate_friction = 1;

real particle_radius = .15;
real particle_mass = .5;
real particle_density = 1.14;
real particle_friction = 1.0;
Vector particle_initial_vel = Vector(0, -5, 0);     //initial velocity

real container_width = 3.0;
real container_thickness = .25;
real container_height = 6.0;
real wscale = 1;

int max_iter = 80;
real tolerance = 1e-2;

real gravity = -9.810;
real timestep = .001;
real seconds_to_simulate = 1;
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

		mSys->Get_bodylist()->at(5)->SetRot(ChQuaternion<>(1, 0, 0, 0));
		mSys->Get_bodylist()->at(5)->SetWvel_loc(ChVector<>(0, 0, 0));
		mSys->Get_bodylist()->at(5)->SetPos(ChVector<>(0, mSys->Get_bodylist()->at(5)->GetPos().y, 0));
		mSys->Get_bodylist()->at(5)->SetPos_dt(ChVector<>(0, mSys->Get_bodylist()->at(5)->GetPos_dt().y, 0));
	if (frame % particles_every == 0) {
		//addHCPSheet(10, 10, 0, particle_mass, particle_radius, particle_friction, true, 0, 0, particle_initial_vel, (ChSystemParallel*) mSys);
	}

	real residual = ((ChLcpSolverParallel *) (mSys->GetLcpSolverSpeed()))->GetResidual() * .5;
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
	omp_set_num_threads(8);

	//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

	//=========================================================================================================
	//ReadInputFile("convergence_config.txt",system_gpu);
	//system_gpu->SetMaxiter(max_iter);
	//system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(20);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(20);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	system_gpu->SetTol(tolerance);
	system_gpu->SetTolSpeeds(tolerance);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(tolerance);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance( 0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(10);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	//((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->solver.SetAPGDParams(10, .9, 2.0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->DoStabilization(false);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .01);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(20, 20, 20));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(100, 50);
	//=========================================================================================================
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
	ReadAllObjectsWithGeometryChrono(system_gpu, "dump.txt");
/*
	//=========================================================================================================

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
	InitObject(BLOCK, 100, Vector(0, 5-.5, 0), quat, material, true, false, -1, -20);

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
	BLOCK->SetPos_dt(ChVector<>(0,-4,0));
	ParticleGenerator layer_gen(system_gpu);
	layer_gen.SetDensity(1000);
	layer_gen.SetRadius(R3(particle_radius));

	layer_gen.material->SetFriction(1);
	layer_gen.material->SetRollingFriction(0);
	layer_gen.material->SetSpinningFriction(0);
	layer_gen.material->SetCohesion(0);
	layer_gen.material->SetCompliance(1e-5);
	int3 num_per_dir = I3(2.0 / particle_radius, 4.0 / particle_radius, 2.0 / particle_radius);

	layer_gen.AddMixtureType(MIX_SPHERE);
	layer_gen.AddMixtureType(MIX_ELLIPSOID);
	//layer_gen.AddMixtureType(MIX_CONE);
	layer_gen.AddMixtureType(MIX_CUBE);
	layer_gen.AddMixtureType(MIX_CYLINDER);

	//layer_gen.SetNormalDistribution(rad.x, rad.x/4.0);
	//layer_gen->UseNormalCohesion(particle_cohesion, 1);

	layer_gen.addPerturbedVolumeMixture(R3(0, -.8, 0), num_per_dir, R3(0, .0, 0), R3(0, -4, 0));

*/
	//=========================================================================================================
	//////Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");

	//openGLView.render_camera->camera_pos = Vector(0, -5, -40);
	//openGLView.render_camera->look_at = Vector(0, -5, 0);
	//openGLView.render_camera->mScale = .1;

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
			ofile << ((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->GetResidual() << endl;
//			stringstream ss;
//			cout << "Frame: " << file << endl;
//			DumpObjects(system_gpu, ss.str());
//			//output.ExportData(ss.str());
//			file++;
		}
		timer.stop();

	}
	DumpAllObjectsWithGeometryChrono(system_gpu, "dump.txt");

	cout << "TIME: " << timer() << endl;
	ofile.close();

	return 0;
}

