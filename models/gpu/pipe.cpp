#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .0005;
real seconds_to_simulate = 15;

int max_iter = 15 * 3;
real tolerance = .001;


real3 container_size = R3(1.75, 1, 4.7);
real container_thickness = .04;
real container_height = 0;
Vector container_pos = Vector(0, container_height,0);
real container_friction = 0;


int num_steps = seconds_to_simulate / timestep;
float particle_radius = .01;
float cohesion = 0;
int read_file = 0;

string data_folder = "data/trough";

ParticleGenerator *layer_gen;

template<class T>
void RunTimeStep(T* mSys, const int frame) {

	double TIME = mSys->GetChTime();
	double STEP = mSys->GetTimerStep();
	double BROD = mSys->GetTimerCollisionBroad();
	double NARR = mSys->GetTimerCollisionNarrow();
	double LCP = mSys->GetTimerLcp();
	double UPDT = mSys->GetTimerUpdate();
	double RESID = ((ChLcpSolverParallel *) (mSys->GetLcpSolverSpeed()))->GetResidual();
	int BODS = ((ChSystemParallel*) mSys)->GetNbodies();
	int CNTC = ((ChSystemParallel*) mSys)->GetNcontacts();
	int REQ_ITS = ((ChLcpSolverParallel*) (mSys->GetLcpSolverSpeed()))->GetTotalIterations();

	printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

}

int main(int argc, char* argv[]) {
	if (argc == 2) {
		omp_set_num_threads(atoi(argv[1]));

	}

	if (argc == 4) {
		omp_set_num_threads(atoi(argv[1]));
		cohesion = atof(argv[2]);
		data_folder = argv[3];
	}
//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	ChCollisionSystemParallel *mcollisionengine = new ChCollisionSystemParallel();
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(tolerance);
	system_gpu->SetTolSpeeds(tolerance);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(tolerance);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(8);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->DoStabilization(false);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .02);
	mcollisionengine->setBinsPerAxis(I3(100, 100, 100));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	real3 container_size = R3(1.75, 1, 4.7);
	real container_thickness = .04;
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);
	material->SetCompliance(0);

	InitObject(Bottom, 100000, Vector(0, container_height, 0) + container_pos, Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometryTriangleMesh(Bottom, "pipe.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);

	layer_gen = new ParticleGenerator(system_gpu);
	layer_gen->SetDensity(1000);
	layer_gen->SetRadius(R3(particle_radius));
	layer_gen->SetCylinderRadius(1.5);
	//layer_gen->SetNormalDistribution(particle_radius, .005);
	layer_gen->material->SetFriction(0);
	layer_gen->material->SetCohesion(0);
	layer_gen->material->SetRollingFriction(0);
	layer_gen->material->SetSpinningFriction(0);
	//layer_gen->AddMixtureType(MIX_TYPE1);
	//layer_gen->AddMixtureType(MIX_TYPE2);
	//layer_gen->AddMixtureType(MIX_TYPE3);
	//layer_gen->AddMixtureType(MIX_TYPE4);
	layer_gen->AddMixtureType(MIX_SPHERE);
	//layer_gen->AddMixtureType(MIX_ELLIPSOID);
	//layer_gen->AddMixtureType(MIX_DOUBLESPHERE);
	int3 num_per_dir = I3(120, 140, 120);
	//120 140 120
	layer_gen->addPerturbedVolumeMixture(R3(0, container_height+2.5 , 0 ), num_per_dir, R3(.01, .01, .01), R3(0));

//=========================================================================================================
//Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
//	openGLView.render_camera->look_at = Vector(0, -5, 0);
//	openGLView.render_camera->mScale = .1;
//	openGLView.SetCustomCallback(RunTimeStep);
//	openGLView.StartSpinning(window_manager);
//	window_manager->CallGlutMainLoop();
//=========================================================================================================
	int file = 0;
	for (int i = 0; i < num_steps; i++) {
		system_gpu->DoStepDynamics(timestep);

		int save_every = 1.0 / timestep / 60.0;     //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
