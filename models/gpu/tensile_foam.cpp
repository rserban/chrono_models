#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 30;

int max_iter = 30;
real tolerance = .0;
int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(6, 6, 6);
real container_thickness = .4;
real container_height = 0;
real container_friction = 1;

real particle_radius = .1;

template<class T>
void RunTimeStep(T* mSys, const int frame) {

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(1);
	material->SetSpinningFriction(1);
	material->SetRollingFriction(1);
	material->SetCohesion(0);
	material->SetCompliance(0);

//	if (mSys->GetNbodies() < 1000) {
//		ChSharedBodyPtr sphere;
//		real3 rad = R3(particle_radius * 2, particle_radius * 3, particle_radius * 2);
//		real3 size = container_size;
//		size.y = container_size.y / 3.0;
//
//		int3 num_per_dir = I3(4, 1, 4);
//		real mu = 1;
//		real mass = 1;
//		real3 vel = R3(0, -5, 0);
//		if (frame % 80 == 0) {
//			ParticleGenerator layer_gen(((ChSystemParallel*) mSys));
//			layer_gen.SetDensity(1000);
//			layer_gen.SetRadius(R3(particle_radius));
//			layer_gen.SetNormalDistribution(particle_radius, particle_radius*.1);
//			layer_gen.material->SetFriction(1);
//			layer_gen.material->SetRollingFriction(1);
//			layer_gen.material->SetSpinningFriction(1);
//			layer_gen.material->SetCohesion(0);
//			layer_gen.AddMixtureType(MIX_SPHERE);
//			layer_gen.AddMixtureType(MIX_ELLIPSOID);
//			//layer_gen.AddMixtureType(MIX_DOUBLESPHERE);
//
//			layer_gen.addPerturbedVolumeMixture(
//					R3(0, 0, 0),
//					I3(4, 1, 4),
//					R3(.1, .1, .1),
//					R3(0,-5,0));
//
////addPerturbedLayer(R3(-2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.x, friction.x, cohesion.x, R3(0, 5, 0), (ChSystemParallel*) mSys);
////addPerturbedLayer(R3(5, 0, 0), CONE, rad, num_per_dir, R3(0, 0, 0), 1, 1, 0, R3(-5, 0, 0), (ChSystemParallel*) mSys, 0);
////addPerturbedLayer(R3(2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.z, friction.z, cohesion.z, R3(0, 5, 0), (ChSystemParallel*) mSys);
//		}
//	}

}

int main(int argc, char* argv[]) {
	omp_set_num_threads(4);
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
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(40);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(R3(30, 30, 30));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);

	material->SetRollingFriction(.01);
	material->SetSpinningFriction(.01);
	material->SetCompliance(0);

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Top, 100000, Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	L->GetMaterialSurface()->SetCohesion(0);
	R->GetMaterialSurface()->SetCohesion(0);
	F->GetMaterialSurface()->SetCohesion(0);
	B->GetMaterialSurface()->SetCohesion(0);
	Bottom->GetMaterialSurface()->SetCohesion(0);
	Top->GetMaterialSurface()->SetCohesion(0);

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);
	FinalizeObject(Top, (ChSystemParallel *) system_gpu);

	ChSharedPtr<ChMaterialSurface> material2;
	material2 = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material2->SetFriction(1);

	material2->SetRollingFriction(.1);
	material2->SetSpinningFriction(.1);
	material2->SetCohesion(100);
	material2->SetCompliance(0);

//
	 ChSharedBodyPtr body;
	 body = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	 Quaternion q(1,0,0,0);
	// q.Q_from_AngX(PI);
	 InitObject(body, 1, Vector(0, 1, 0), q, material2, true, false, -1, 1);

	 AddCollisionGeometry(body, SPHERE, ChVector<>(1,2,1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	 body->SetWvel_loc(Vector(0,20,0));
	 //body->SetPos_dt(Vector(5,0,0));
	 FinalizeObject(body, (ChSystemParallel *) system_gpu);

//=========================================================================================================
//Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
	openGLView.render_camera->look_at = Vector(0, -5, 0);
	openGLView.render_camera->mScale = .5;
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
