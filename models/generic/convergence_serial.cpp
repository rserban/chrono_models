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

real particle_radius = .2;
real particle_mass = 1;
real particle_friction = 1.0;

real container_width = 3.0;
real container_thickness = .25;
real container_height = 6.0;
real wscale = 1;

int max_iter = 1000;
real tolerance = .0001;

real gravity = -9.810;
real timestep = .001;
real seconds_to_simulate = 3;
int num_steps = seconds_to_simulate / timestep;

int save_every = 1.0 / timestep / 100.0;     //save data every n steps
ChSharedBodyPtr BLOCK;

string data_folder = "data/convergence";

template<class T>
void RunTimeStep(T* mSys, const int frame) {

	BLOCK->SetRot(ChQuaternion<>(1, 0, 0, 0));
	BLOCK->SetWvel_loc(ChVector<>(0, 0, 0));
	BLOCK->SetPos(ChVector<>(0, BLOCK->GetPos().y, 0));
	BLOCK->SetPos_dt(ChVector<>(0, BLOCK->GetPos_dt().y, 0));
	//if (frame % particles_every == 0) {
	//addHCPSheet(10, 10, 0, particle_mass, particle_radius, particle_friction, true, 0, 0, particle_initial_vel, (ChSystemParallel*) mSys);
	//}

//	std::vector<double> violation = ((ChLcpIterativeSolver*) mSys->GetLcpSolverSpeed())->GetViolationHistory();
//	for (int i = 0; i < violation.size(); i++) {
//		cout << violation[i] << endl;
//	}
	//double residual =violation[violation.size()-1];
	//cout<<residual<<endl;
	//real residual = ((ChLcpSolverParallel *) (mSys->GetLcpSolverSpeed()))->GetResidual() * .5;
	//mSys->SetTol(residual);
	//mSys->SetTolSpeeds(residual);
	//((ChLcpSolverParallel *) (mSys->GetLcpSolverSpeed()))->SetTolerance(residual);

//	if(frame >300){
//
//		BLOCK->SetPos(ChVector<>(0, 5, 0));
//
//	}

}
int main(int argc, char* argv[]) {
	omp_set_num_threads(8);
	string solver = argv[1];

	max_iter = atof(argv[2]);
	tolerance = atof(argv[3]);
	real recovery_speed = atof(argv[4]);
	particle_radius = atof(argv[5]);
	real block_mass = atof(argv[6]);
	data_folder = argv[7];

	cout << solver << " " << max_iter << " " << tolerance << " " << particle_radius << " " << data_folder << endl;

	//=========================================================================================================
	ChSystem * system = new ChSystem;
	system->SetIntegrationType(ChSystem::INT_ANITESCU);
	if (solver == "APGD") {
		system->SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
	} else if (solver == "JACOBI") {
		system->SetLcpSolverType(ChSystem::LCP_ITERATIVE_JACOBI);
	} else if (solver == "SOR") {
		system->SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR);
	}
	//=========================================================================================================
	//ReadInputFile("convergence_config.txt",system_gpu);
	system->SetMaxiter(max_iter);
	system->SetIterLCPmaxItersSpeed(max_iter);
	system->SetTol(tolerance);
	system->SetTolSpeeds(tolerance);
	system->SetMaxPenetrationRecoverySpeed(recovery_speed);

	//=========================================================================================================
	system->Set_G_acc(ChVector<>(0, gravity, 0));
	system->SetStep(timestep);
	((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetRecordViolation(true);
	//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody());
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody());
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody());
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody());
	ChSharedBodyPtr BTM = ChSharedBodyPtr(new ChBody());
	BLOCK = ChSharedBodyPtr(new ChBody());

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(plate_friction);
	material->SetCompliance(0);

	InitObject(L, 100000, Vector(-container_width + container_thickness, plate_height, 0), quat, material, true, true, 1, 1);
	InitObject(R, 100000, Vector(container_width - container_thickness, plate_height, 0), quat, material, true, true, 1, 1);
	InitObject(F, 100000, Vector(0, plate_height, -container_width + container_thickness), quat, material, true, true, 1, 1);
	InitObject(B, 100000, Vector(0, plate_height, container_width - container_thickness), quat, material, true, true, 1, 1);
	InitObject(BTM, 100000, Vector(0, plate_height - container_height, 0), quat, material, true, true, 1, 1);
	InitObject(BLOCK, block_mass, Vector(0, 5, 0), quat, material, true, false, 1, 1);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
	AddCollisionGeometry(F, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
	AddCollisionGeometry(B, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
	AddCollisionGeometry(BTM, BOX, Vector(container_width * wscale, container_thickness, container_width), lpos, quat);
	AddCollisionGeometry(BLOCK, BOX, Vector(container_width, container_thickness, container_width), lpos, quat);

	FinalizeObject(L, system);
	FinalizeObject(R, system);
	FinalizeObject(F, system);
	FinalizeObject(B, system);
	FinalizeObject(BTM, system);
	FinalizeObject(BLOCK, system);

	ParticleGeneratorCPU layer_gen(system);
	layer_gen.SetMass(particle_mass);
	layer_gen.SetRadius(R3(particle_radius));

	layer_gen.material->SetFriction(1);
	layer_gen.material->SetRollingFriction(0);
	layer_gen.material->SetSpinningFriction(0);
	layer_gen.material->SetCohesion(0);
	layer_gen.material->SetCompliance(0);
	int3 num_per_dir = I3(2.0 / particle_radius, 4.0 / particle_radius, 2.0 / particle_radius);

	layer_gen.addPerturbedVolume(R3(0, -.2, 0), SPHERE, num_per_dir, R3(0.0, 0.0, 0.0), R3(0, -4, 0), false);

	//=========================================================================================================
	//////Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system, 800, 600, 0, 0, "Test_Solvers");
//
//	//openGLView.render_camera->camera_pos = Vector(0, -5, -40);
//	//openGLView.render_camera->look_at = Vector(0, -5, 0);
//	//openGLView.render_camera->mScale = .1;
//
//	openGLView.SetCustomCallback(RunTimeStep);
//	openGLView.StartSpinning(window_manager);
//	window_manager->CallGlutMainLoop();
	//=========================================================================================================
	ofstream ofile("convergence.txt");
	ChTimer<double> timer;
	timer.start();
	int file = 0;
	for (int i = 0; i < num_steps; i++) {

		cout << "step " << i << endl;
		system->DoStepDynamics(timestep);
		RunTimeStep(system, i);

		//if (i % save_every == 0) {
		stringstream ss;
		cout << "Frame: " << file << endl;
		ss << data_folder << "/" << file << ".txt";
		DumpResidualHist(system, ss.str());
		file++;
		//}
		timer.stop();

	}
	cout << "TIME: " << timer() << endl;
	ofile.close();
	return 0;
}

