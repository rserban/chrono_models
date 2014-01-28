#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"

real scale_factor = 1000;
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 30;
real tolerance = .0005;

int max_iter = 30;
int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(.4, .2, .4);
real container_thickness = .1;
real container_height = 0;
real container_friction = 0;
real container_cohesion = -1000;

real particle_radius = .003;
real particle_density = 2650;
real particle_slide_friction = .3;
real particle_roll_friction = .3;
real particle_cohesion = 0;

int particle_grid_x = 2;
int particle_grid_z = 2;
real start_height = 1;

ChSharedBodyPtr BLOCK;
ParticleGenerator<ChSystemParallel>* layer_gen;

template<class T>
void RunTimeStep(T* mSys, const int frame) {

	BLOCK->SetRot(ChQuaternion<>(1, 0, 0, 0));
	BLOCK->SetWvel_loc(ChVector<>(0, 0, 0));
	BLOCK->SetPos(ChVector<>(0, BLOCK->GetPos().y, 0));
	BLOCK->SetPos_dt(ChVector<>(0, BLOCK->GetPos_dt().y, 0));

	real cont_vol = .4*(BLOCK->GetPos().y+container_size.y-2*container_thickness)*.4;
	cout << layer_gen->total_volume << " " << cont_vol<<" "<<layer_gen->total_mass/cont_vol << endl;

}

int main(int argc, char* argv[]) {
	//omp_set_num_threads(8);
//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(max_iter);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(max_iter / 2);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(max_iter / 2);
	system_gpu->SetTol(particle_radius);
	system_gpu->SetTolSpeeds(particle_radius);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(particle_radius);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(.5);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(30, 30, 30));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);
	material->SetRollingFriction(container_friction);
	material->SetSpinningFriction(container_friction);
	material->SetCompliance(0);
	material->SetCohesion(-100);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -19, -20);
	InitObject(Top, 100000, Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	L->GetMaterialSurface()->SetCohesion(container_cohesion);
	R->GetMaterialSurface()->SetCohesion(container_cohesion);
	F->GetMaterialSurface()->SetCohesion(container_cohesion);
	B->GetMaterialSurface()->SetCohesion(container_cohesion);
	Bottom->GetMaterialSurface()->SetCohesion(container_cohesion);
	Top->GetMaterialSurface()->SetCohesion(container_cohesion);

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);
	FinalizeObject(Top, (ChSystemParallel *) system_gpu);

	BLOCK = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(BLOCK, 1, Vector(0, container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, false, -1, -20);
	AddCollisionGeometry(BLOCK, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(BLOCK, (ChSystemParallel *) system_gpu);

	layer_gen = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
	layer_gen->SetDensity(particle_density);
	layer_gen->SetRadius(R3(particle_radius));
	layer_gen->SetNormalDistribution(particle_radius, particle_radius / 3.0);
	layer_gen->material->SetFriction(particle_slide_friction);
	layer_gen->material->SetCohesion(particle_cohesion);
	layer_gen->material->SetRollingFriction(particle_roll_friction);
	layer_gen->material->SetSpinningFriction(particle_roll_friction);
	layer_gen->AddMixtureType(MIX_SPHERE);
	layer_gen->AddMixtureType(MIX_ELLIPSOID);
	//layer_gen->AddMixtureType(MIX_DOUBLESPHERE);

	layer_gen->addPerturbedVolumeMixture(R3(0, 0, 0), I3(40, 12, 40), R3(.01, .01, .01), R3(0, 0, 0));

//=========================================================================================================
//Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	//openGLView.render_camera->camera_pos = Vector(0, -5, -10);
	//openGLView.render_camera->look_at = Vector(0, -5, 0);
	openGLView.render_camera->mScale = .4;
	openGLView.SetCustomCallback(RunTimeStep);
	openGLView.StartSpinning(window_manager);
	window_manager->CallGlutMainLoop();
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
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << "data/foam/" << "/" << file << ".txt";
			//DumpAllObjects(system_gpu, ss.str(), ",", true);
			DumpAllObjectsWithGeometry(system_gpu, ss.str(), ",");
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
