#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"

real gravity = -9806.65;
real timestep = .0005;
real seconds_to_simulate = 1.5;
real tolerance = .0005;

int max_iter = 10;
int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(200, 50, 200);
real container_thickness = 10;
real container_height = 0;
real container_friction = .1;
real container_cohesion = -1000;

real particle_radius = 1.5;
real particle_density = .00265;
real particle_slide_friction = .3;
real particle_roll_friction = .3;
real particle_cohesion = 0;
real particle_std_dev = particle_radius / 6.0;

ChSharedBodyPtr BLOCK, CONTAINER;
ParticleGenerator<ChSystemParallel>* layer_gen;
real amplitude = particle_radius;
real frequency = 10;

template<class T>
void RunTimeStep(T* mSys, const int frame) {

	real t = frame * timestep * PI * 2 * frequency;

	BLOCK->SetRot(ChQuaternion<>(1, 0, 0, 0));
	BLOCK->SetWvel_loc(ChVector<>(0, 0, 0));
	BLOCK->SetPos(ChVector<>(sin(t) * amplitude, BLOCK->GetPos().y, 0));
	BLOCK->SetPos_dt(ChVector<>(cos(t) * amplitude * 2 * PI * frequency, BLOCK->GetPos_dt().y, 0));

	CONTAINER->SetPos(ChVector<>(sin(t) * amplitude, 0, 0));
	CONTAINER->SetPos_dt(ChVector<>(cos(t) * amplitude * 2 * PI * frequency, 0, 0));
	CONTAINER->SetWvel_loc(ChVector<>(0, 0, 0));
	CONTAINER->SetRot(ChQuaternion<>(1, 0, 0, 0));

	real cont_vol = .4 * (BLOCK->GetPos().y + container_size.y - 2 * container_thickness) * .4;
	cout << layer_gen->total_volume<<" "<<layer_gen->total_mass << " " << cont_vol << " " << layer_gen->total_mass / cont_vol << endl;

}

int main(int argc, char* argv[]) {

	if (argc > 1) {
		particle_slide_friction = atof(argv[1]);
		particle_roll_friction = atof(argv[2]);
		particle_std_dev = atof(argv[3]);
		amplitude = atof(argv[4]);
	}

//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	system_gpu->SetTol(particle_radius);
	system_gpu->SetTolSpeeds(particle_radius);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(particle_radius);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(20);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(30, 30, 30));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);
	material->SetRollingFriction(container_friction);
	material->SetSpinningFriction(container_friction);
	material->SetCompliance(0);
	material->SetCohesion(-100);

	CONTAINER = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(CONTAINER, 100000, Vector(0, 0, 0), Quaternion(1, 0, 0, 0), material, true, false, -20, -20);
	AddCollisionGeometry(CONTAINER, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(-container_size.x + container_thickness, container_height - container_thickness, 0),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(CONTAINER, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(container_size.x - container_thickness, container_height - container_thickness, 0),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(CONTAINER, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, container_height - container_thickness, -container_size.z + container_thickness),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(CONTAINER, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, container_height - container_thickness, container_size.z - container_thickness),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(CONTAINER, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(CONTAINER, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0));

	CONTAINER->GetMaterialSurface()->SetCohesion(container_cohesion);
	FinalizeObject(CONTAINER, (ChSystemParallel *) system_gpu);

	BLOCK = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(BLOCK, 100, Vector(0, container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, false, -1, -20);
	AddCollisionGeometry(BLOCK, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(BLOCK, (ChSystemParallel *) system_gpu);

	layer_gen = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
	layer_gen->SetDensity(particle_density);
	layer_gen->SetRadius(R3(particle_radius, particle_radius * .5, particle_radius));
	layer_gen->SetNormalDistribution(particle_radius, particle_std_dev);
	layer_gen->material->SetFriction(particle_slide_friction);
	layer_gen->material->SetCohesion(particle_cohesion);
	layer_gen->material->SetRollingFriction(0);
	layer_gen->material->SetSpinningFriction(0);
	layer_gen->AddMixtureType(MIX_SPHERE);
	layer_gen->AddMixtureType(MIX_ELLIPSOID);
	//layer_gen->AddMixtureType(MIX_DOUBLESPHERE);

	layer_gen->addPerturbedVolumeMixture(R3(0, 0, 0), I3(85, 20, 85), R3(0,0,0), R3(0, 0, 0));

//=========================================================================================================
//Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	//openGLView.render_camera->camera_pos = Vector(0, -5, -10);
//	//openGLView.render_camera->look_at = Vector(0, -5, 0);
//	openGLView.render_camera->mScale = 20;
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
//		if (i % save_every == 0) {
//			stringstream ss;
//			cout << "Frame: " << file << endl;
//			ss << "data/foam/" << "/" << file << ".txt";
//			//DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
//			file++;
//		}
		RunTimeStep(system_gpu, i);
	}
}
