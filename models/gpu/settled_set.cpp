#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .0005;
real3 particle_radius = R3(.1, .05, .1);
real seconds_to_simulate = 5.0;

int max_iter = 25;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(10, 2, 30);
real container_thickness = .25;
real container_height = 4;
real container_friction = .5;

int3 num_per_dir = I3((container_size - particle_radius * 2) / particle_radius);

template<class T>
void RunTimeStep(T* mSys, const int frame) {

}

int main(int argc, char* argv[]) {
	if (argc == 2) {
		omp_set_num_threads(atoi(argv[1]));

	} else {

		omp_set_num_threads(1);
	}

	//=========================================================================================================
	ChSystemGPU * system_gpu = new ChSystemGPU;
	ChCollisionSystemGPU *mcollisionengine = new ChCollisionSystemGPU();
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);
	//=========================================================================================================
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(1e-3);
	system_gpu->SetTolSpeeds(1e-3);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(1e-3);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(.6);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius.x * .05);
	mcollisionengine->setBinsPerAxis(R3(num_per_dir.x * 2, num_per_dir.y * 2, num_per_dir.z * 2));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
	//=========================================================================================================
	cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
	if (particle_radius.x == particle_radius.y == particle_radius.z) {
		addPerturbedLayer(R3(0, -5 + particle_radius.y + container_thickness, 0), SPHERE, particle_radius, num_per_dir, R3(.01, .01, .01), 4, .1, system_gpu);

	} else {
		addPerturbedLayer(R3(0, -5 + particle_radius.y + container_thickness, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 4, .1, system_gpu);
	}
	//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true,
			-20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20,
			-20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true,
			-20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20,
			-20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);

	//=========================================================================================================
	//////Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_position = glm::vec3(0, -5, -10);
	openGLView.render_camera->camera_look_at = glm::vec3(0, -5, 0);
	openGLView.render_camera->camera_scale = .1;
	openGLView.SetCustomCallback(RunTimeStep);
	openGLView.StartSpinning(window_manager);
	window_manager->CallGlutMainLoop();
	//=========================================================================================================

	for (int i = 0; i < num_steps; i++) {
		system_gpu->DoStepDynamics(timestep);
		double TIME = system_gpu->GetChTime();
		double STEP = system_gpu->GetTimerStep();
		double BROD = system_gpu->GetTimerCollisionBroad();
		double NARR = system_gpu->GetTimerCollisionNarrow();
		double LCP = system_gpu->GetTimerLcp();
		double UPDT = system_gpu->GetTimerUpdate();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverGPU*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS);
		if (i % 1000 == 0) {
			cout << "SAVED STATE" << endl;
			DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");
		}
		RunTimeStep(system_gpu, i);
	}

	DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
