#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real3 particle_radius = R3(.25, 1, .25);
real seconds_to_simulate = 5;

int max_iter = 10;
ChVector<> lpos(0, 0, 0);
ChQuaternion<> quat(1, 0, 0, 0);
int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(50, 10, 50);
real container_thickness = .5;
real container_height = 10;
real container_friction = 1;
real current_time = 0;
int save_every = 1.0 / timestep / 60.0; //save data every n steps

int3 num_per_dir = I3((container_size - particle_radius * 2) / particle_radius);

double frequency = 20;
double amplitude = 4 * fabs(gravity) / (4 * PI * PI * frequency);

int particles_every = 200;
int particle_grid_x = 10;
int particle_grid_z = 10;
real particle_friction = 1;
real start_height = 30;
ChSharedBodyPtr impactor;
ChSharedBodyPtr Bottom;
template<class T>
void RunTimeStep(T* mSys, const int frame) {

	if (frame % particles_every == 0) {
		ChSharedBodyPtr sphere;
		for (int i = 0; i < particle_grid_x; i++) {
			for (int j = 0; j < particle_grid_z; j++) {
				sphere = ChSharedBodyPtr(new ChBodyGPU);
				Quaternion q;
				q.Q_from_NasaAngles(Vector(rand() % 1000 / 1000.0, rand() % 1000 / 1000.0, rand() % 1000 / 1000.0));
				q.Normalize();

				ChVector<> position(i * particle_radius.x * 2 - particle_grid_x * .5 * particle_radius.x * 2, start_height, j * particle_radius.x * 2 - particle_grid_z * .5 * particle_radius.x * 2);
				position.x*=2;
				position.z*=2;
				position.x += rand() % 1000 / 100000.0;
				//position.y += rand() % 1000 / 1000.0;
				position.z += rand() % 1000 / 100000.0;

				InitObject(sphere, 1, position, quat, particle_friction, particle_friction, 0, true, false, -1, i);
				AddCollisionGeometry(sphere, CYLINDER, ChVector<>(particle_radius.x, particle_radius.y, particle_radius.z), Vector(0, 0, 0), quat);
				//AddCollisionGeometry(sphere, ELLIPSOID, ChVector<>(particle_radius.x, particle_radius.x*2, particle_radius.z), Vector(0, particle_radius.y, 0), quat);
				//AddCollisionGeometry(sphere, SPHERE, ChVector<>(particle_radius.x, 0, 0), Vector(particle_radius.x, -particle_radius.y,particle_radius.x ), quat);
				//AddCollisionGeometry(sphere, SPHERE, ChVector<>(particle_radius.x, 0, 0), Vector(-particle_radius.x, -particle_radius.y,particle_radius.x ), quat);

				real3 r = particle_radius;
				real mass = 1;

				Vector inertia = Vector((1 / 5.0 * mass * (r.y * r.y + r.z * r.z), 1 / 5.0 * mass * (r.x * r.x + r.z * r.z), 1 / 5.0 * mass * (r.x * r.x + r.y * r.y)));
				sphere->SetInertiaXX(inertia);
				sphere->SetPos_dt(Vector(0,-15,0));
				FinalizeObject(sphere, (ChSystemGPU *) mSys);
			}
		}
	}

}

int main(int argc, char* argv[]) {
	omp_set_num_threads(1);
	//=========================================================================================================
	ChSystemGPU * system_gpu = new ChSystemGPU;
	ChCollisionSystemGPU *mcollisionengine = new ChCollisionSystemGPU();
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

	//=========================================================================================================
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(1e-8);
	system_gpu->SetTolSpeeds(1e-8);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(1e-8);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(1e-5, 1e-5, .2);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(0.6);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(CONJUGATE_GRADIENT);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius.x * .05);
	mcollisionengine->setBinsPerAxis(R3(num_per_dir.x , num_per_dir.y , num_per_dir.z ));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBodyGPU);
	Bottom = ChSharedBodyPtr(new ChBodyGPU);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true,
			-20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20,
			-20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true,
			-20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20,
			-20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness * 4, container_size.y * 1.5, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness * 4, container_size.y * 1.5, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y * 1.5, container_thickness * 4), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y * 1.5, container_thickness * 4), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);

//	impactor = ChSharedBodyGPUPtr(new ChBodyGPU);
//	InitObject(impactor, 1500, Vector(-container_size.x, container_height + container_size.y * 2, 0), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, -1, -2);
//	AddCollisionGeometry(impactor, SPHERE, ChVector<>(.5, 0, 0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
//	FinalizeObject(impactor, (ChSystemGPU *) system_gpu);
//	impactor->SetPos_dt(Vector(2.5, 0, 0));
	//=========================================================================================================
	ChSharedBodyPtr A1 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A2 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A3 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A4 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A5 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A6 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A7 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A8 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A9 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A10 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A11 = ChSharedBodyPtr(new ChBodyGPU);
	ChSharedBodyPtr A12 = ChSharedBodyPtr(new ChBodyGPU);

	Quaternion q7, q8, q10, q11;
	q7.Q_from_AngZ(10*PI/180.0);
	q8.Q_from_AngZ(-10*PI/180.0);
	q10.Q_from_AngZ(15*PI/180.0);
	q11.Q_from_AngZ(107.5*PI/180.0);

	q7.Normalize();
	q8.Normalize();
	q10.Normalize();
	q11.Normalize();

	InitObject(A1, 100000, Vector(0, 11.75, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A2, 100000, Vector(2.07, 12.5, -2.88), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A3, 100000, Vector(-2.07, 12.5, -2.88), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A4, 100000, Vector(5.2, 13.5, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A5, 100000, Vector(-5.2, 13.5, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A6, 100000, Vector(0, 6.3, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A7, 100000, Vector(2.8, 5.75, -0.45), q7, container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A8, 100000, Vector(-2.8, 5.75, -0.45), q8, container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A9, 100000, Vector(-0.15, 15.75, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A10, 100000, Vector(-0.45, 16.6, 0), q10, container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A11, 100000, Vector(-2.85, 17.9, 0), q11, container_friction, container_friction, 0, true, true, -20, -20);
	InitObject(A12, 100000, Vector(-4.2, 17.5, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);

	AddCollisionGeometry(A1, ELLIPSOID, Vector(6, 4, 3.5), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A2, ELLIPSOID, Vector(0.9, 0.9, 0.9), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A3, ELLIPSOID, Vector(0.9, 0.9, 0.9), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A4, ELLIPSOID, Vector(1.5, 1.5, 0.5), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A5, ELLIPSOID, Vector(1.5, 1.5, 0.5), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A6, ELLIPSOID, Vector(3, 6, 1.8), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A7, ELLIPSOID, Vector(1, 2.4, 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A8, ELLIPSOID, Vector(1, 2.4, 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A9, ELLIPSOID, Vector(1, 0.25, 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A10, ELLIPSOID, Vector(0.4, 2, 0.4), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A11, ELLIPSOID, Vector(0.4, 2, 0.4), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(A12, ELLIPSOID, Vector(1, 1, 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(A1, (ChSystemGPU *) system_gpu);
	FinalizeObject(A2, (ChSystemGPU *) system_gpu);
	FinalizeObject(A3, (ChSystemGPU *) system_gpu);
	FinalizeObject(A4, (ChSystemGPU *) system_gpu);
	FinalizeObject(A5, (ChSystemGPU *) system_gpu);
	FinalizeObject(A6, (ChSystemGPU *) system_gpu);
	FinalizeObject(A7, (ChSystemGPU *) system_gpu);
	FinalizeObject(A8, (ChSystemGPU *) system_gpu);
	FinalizeObject(A9, (ChSystemGPU *) system_gpu);
	FinalizeObject(A10, (ChSystemGPU *) system_gpu);
	FinalizeObject(A11, (ChSystemGPU *) system_gpu);
	FinalizeObject(A12, (ChSystemGPU *) system_gpu);
	//=========================================================================================================
	//////Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_pos = Vector(0, 0, -10);
	openGLView.render_camera->look_at = Vector(0, 0, 0);
	openGLView.render_camera->mScale = 1;
	openGLView.SetCustomCallback(RunTimeStep);
	openGLView.StartSpinning(window_manager);
	window_manager->CallGlutMainLoop();
	//=========================================================================================================
	int file = 0;

	stringstream ss_m;
	string data_folder = "data";

//	ss_m << data_folder << "/" << "timing.txt";
//	string timing_file_name = ss_m.str();
//	ofstream ofile(timing_file_name.c_str());
//	ofile.close();

	for (int i = 0; i < num_steps; i++) {
		system_gpu->DoStepDynamics(timestep);
		cout << "step " << i;
		cout << " Residual: " << ((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		cout << " ITER: " << ((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();
		cout << " OUTPUT STEP: Time= " << current_time << " bodies= " << system_gpu->GetNbodies() << " contacts= " << system_gpu->GetNcontacts() << " step time=" << system_gpu->GetTimerStep()
				<< " lcp time=" << system_gpu->GetTimerLcp() << " CDbroad time=" << system_gpu->GetTimerCollisionBroad() << " CDnarrow time=" << system_gpu->GetTimerCollisionNarrow() << " Iterations="
				<< ((ChLcpSolverGPU*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations() << "\n";
		//TimingFile(system_gpu, timing_file_name, current_time);
		current_time +=timestep;
//		if (i % save_every == 0) {
//			stringstream ss;
//			cout << "Frame: " << file << endl;
//			ss << data_folder << "/" << file << ".txt";
//			DumpObjects(system_gpu, ss.str());
//			//output.ExportData(ss.str());
//			file++;
//		}
		RunTimeStep(system_gpu, i);
	}

	DumpAllObjects(system_gpu, "Cylinder_Drop_State.txt");

}
