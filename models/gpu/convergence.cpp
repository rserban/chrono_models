#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"

ChVector<> lpos(0, 0, 0);
ChQuaternion<> quat(1, 0, 0, 0);

//all dimensions are in millimeters, milligrams
real plate_height = -10;
real plate_thickness = 1;
real plate_radius = 7;
real plate_friction = 1;

real particle_radius = .1;
real particle_mass = .5263;
real particle_density = 1.14;
real particle_friction = 1.0;
Vector particle_initial_vel = Vector(0, -5, 0); //initial velocity

real container_width = 3.0;
real container_thickness = .25;
real container_height = 6.0;
real wscale = 1;

real gravity = -9.810;
real timestep = .001;
int num_steps = 5000000;

int particle_grid_x = 10;
int particle_grid_z = 10;
int particles_every = 50; //add particles every n steps
int save_every = 100; //save data every n steps

int particle_configuration = 0;
//0: single sphere
//1: two spheres joined together
template<class T>
void RunTimeStep(T* mSys, const int frame) {
	cout << frame << endl;
	if (frame % particles_every == 0) {
		//addHCPSheet(10, 10, 0, particle_mass, particle_radius, particle_friction, true, 0, 0, particle_initial_vel, (ChSystemGPU*) mSys);

	}
}
int main(int argc, char* argv[]) {
	//omp_set_num_threads(1);

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
	ReadInputFile("convergence_config.txt",system_gpu);
	//=========================================================================================================
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	//system_gpu->SetStep(timestep);
	//=========================================================================================================
//	Quaternion plate_quat;
//	plate_quat.Q_from_AngAxis(0, Vector(1, 0, 0));
//
//	ChSharedBodyGPUPtr PLATE = ChSharedBodyGPUPtr(new ChBody(new ChCollisionModelGPU));
//	InitObject(PLATE, 1, ChVector<>(0, plate_height, 0), plate_quat, plate_friction, plate_friction, 0, true, true, -1000, -20000);
//	AddCollisionGeometry(PLATE, BOX, ChVector<>(plate_radius, plate_thickness, plate_radius), lpos, quat);
//	FinalizeObject(PLATE, (ChSystemGPU *) system_gpu);

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr BTM = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

	InitObject(L, 100000, Vector(-container_width + container_thickness, plate_height, 0), quat, plate_friction, plate_friction, 0, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_width - container_thickness, plate_height, 0), quat, plate_friction, plate_friction, 0, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, plate_height, -container_width + container_thickness), quat, plate_friction, plate_friction, 0, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, plate_height, container_width - container_thickness), quat, plate_friction, plate_friction, 0, true, true, -20, -20);
	InitObject(BTM, 100000, Vector(0, plate_height-container_height, 0), quat, plate_friction, plate_friction, 0, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
	AddCollisionGeometry(F, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
	AddCollisionGeometry(B, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
	AddCollisionGeometry(BTM, BOX, Vector(container_width * wscale, container_thickness, container_width), lpos, quat);

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(BTM, (ChSystemGPU *) system_gpu);

	addHCPCube(23, 1, 23, particle_mass, particle_radius, particle_friction, true, 0, -15, 0, Vector(0, 0, 0), system_gpu);

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

//	int file = 0;
//	for (int i = 0; i < num_steps; i++) {
//		cout << "step " << i << endl;
//		system_gpu->DoStepDynamics(timestep);
//		RunTimeStep(system_gpu, i);
//		if (i % save_every == 0) {
//			stringstream ss;
//			cout << "Frame: " << file << endl;
//			DumpObjects(system_gpu, ss.str());
//			//output.ExportData(ss.str());
//			file++;
//		}
//
//	}
	return 0;
}

