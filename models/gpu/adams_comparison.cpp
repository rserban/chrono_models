#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = 0;
real timestep = .001;
real seconds_to_simulate = timestep*5;

int max_iter = 100;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(6, 6, 6);
real container_thickness = .4;
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
ParticleGenerator<ChSystemParallel> *layer_gen;

template<class T>
void RunTimeStep(T* mSys, const int frame) {


}

int main(int argc, char* argv[]) {
	int threads = 8;


	if (argc > 1) {
		threads=(atoi(argv[1]));
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
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(0);
	system_gpu->SetTolSpeeds(0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(1000);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .02);
	mcollisionengine->setBinsPerAxis(I3(50, 50, 50));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(0);
	material->SetRollingFriction(0);
	material->SetSpinningFriction(0);
	material->SetCompliance(0);
	material->SetCohesion(0);

	InitObject(L, 1, Vector(-.09,0,0), Quaternion(1, 0, 0, 0), material, true, false, -20, 1);
	InitObject(R, 1, Vector(.09,0,0), Quaternion(1, 0, 0, 0), material, true, false, -20, 1);

	AddCollisionGeometry(L, SPHERE, Vector(.1 ,1.,1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, SPHERE, Vector(.1,.1,.1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));


	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	//L->SetPos_dt(Vector(.5,0,0));
	//R->SetPos_dt(Vector(-.5,0,0));


//=========================================================================================================
//Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
//	openGLView.render_camera->look_at = Vector(0, 0, 0);
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
		double RESID = ((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

//		int save_every = 1.0 / timestep / 60.0;     //save data every n steps
//		if (i % save_every == 0) {
//			stringstream ss;
//			cout << "Frame: " << file << endl;
//			ss << "data/extruder/" << "/" << file << ".txt";
//			//DumpAllObjects(system_gpu, ss.str(), ",", true);
//			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
//			//output.ExportData(ss.str());
//			file++;
//		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
