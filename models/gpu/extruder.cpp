#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .0005;
real seconds_to_simulate = 15;

int max_iter = 30;

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
ParticleGenerator *layer_gen;
string data_folder = "data/extruder";
ChSharedPtr<ChMaterialSurface> material_fiber;
template<class T>
void CreateFiber(T* mSys, ChVector<> position) {

	real length = .1;
	real thickness = .04;

	ChSharedBodyPtr Fiber1 = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	ChSharedBodyPtr Fiber2;
	ChSharedBodyPtr fibers[10];
	fibers[0] = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	Quaternion q;
	q.Q_from_AngZ(PI / 2.0);
	InitObject(fibers[0], .1, Vector(0, -2, 0) + position, q, material_fiber, true, false, -1, -2);
	AddCollisionGeometry(fibers[0], CYLINDER, Vector(thickness, length, length), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(fibers[0], SPHERE, Vector(thickness, length, length), Vector(0, length, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(fibers[0], SPHERE, Vector(thickness, length, length), Vector(0, -length, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(fibers[0], (ChSystemParallel *) mSys);
	fibers[0]->SetPos_dt(Vector(0, -1, 0));

	for (int i = 1; i < 10; i++) {
		fibers[i] = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
		InitObject(fibers[i], .1, Vector(i * (length + thickness) * 2, -2, 0) + position, q, material_fiber, true, false, -1, -2);
		AddCollisionGeometry(fibers[i], CYLINDER, Vector(thickness, length, length), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
		AddCollisionGeometry(fibers[i], SPHERE, Vector(thickness, length, length), Vector(0, length, 0), Quaternion(1, 0, 0, 0));
		AddCollisionGeometry(fibers[i], SPHERE, Vector(thickness, length, length), Vector(0, -length, 0), Quaternion(1, 0, 0, 0));
		FinalizeObject(fibers[i], (ChSystemParallel *) mSys);
		fibers[i]->SetPos_dt(Vector(0, -1, 0));
		ChCoordsys<> pos;
		pos.pos = Vector((i - 1) * (length + thickness) * 2 + length, -2, 0) + position;
		ChSharedPtr<ChLinkLockRevolute> joint(new ChLinkLockRevolute);
		joint->Initialize(fibers[i - 1], fibers[i], pos);
		//		joint->Initialize(fibers[i - 1], fibers[i], true, Vector(length,0,0),Vector(-length,0,0),true);
		//		joint->Set_SpringF(100);
		//		joint->Set_SpringK(100);
		//		joint->Set_SpringR(100);

		mSys->AddLink(joint);

	}

}

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	if (mSys->GetNbodies() < 1071630) {
		ChSharedBodyPtr sphere;
		real3 rad = R3(particle_radius, particle_radius, particle_radius);
		real3 size = container_size;
		size.y = container_size.y / 3.0;

		//int3 num_per_dir = I3(1, 10, 10);

		if (frame % 50 == 0 && frame * timestep < 2.0) {

			//addPerturbedLayer(R3(-2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.x, friction.x, cohesion.x, R3(0, 5, 0), (ChSystemParallel*) mSys);

			layer_gen->addPerturbedVolumeMixture(R3(0, 0, 0), I3(100, 1, 100), R3(0, 0, 0), R3(0, -5, 0));

			//layer_gen->addPerturbedVolume(R3(0, 0, 0), SPHERE, I3(40, 1, 40), R3(0, 0, 0), R3(0, -5, 0));

			//layer_gen->addPerturbedVolume(R3(-2, 0, 0), SPHERE, I3(15, 1, 15), R3(0, 0, 0), R3(0, -5, 0));
			//layer_gen->SetRadius(R3(particle_radius,particle_radius*2.0,particle_radius));
			//layer_gen->addPerturbedVolume(R3(0, 0, 2), SPHERE, I3(15, 1, 15), R3(0, 0, 0), R3(0, -5, 0));
			///layer_gen->SetRadius(R3(particle_radius,particle_radius/2.0,particle_radius));
			//layer_gen->addPerturbedVolume(R3(0, 0, -2), SPHERE, I3(15, 1, 15), R3(0, 0, 0), R3(0, -5, 0));
			//addPerturbedLayer(R3(2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.z, friction.z, cohesion.z, R3(0, 5, 0), (ChSystemParallel*) mSys);
		}
		if (frame % 150 == 0 && frame * timestep > 3 && frame * timestep < 5) {

			for (int i = 0; i < 10; i++) {
				//CreateFiber(mSys, Vector(0, 3, i / 8.0));

			}

		}
	}
//
//	slicer1->SetPos(Vector(3.5,(sin(frame*timestep*6)-1),0));
//	slicer1->SetPos_dt(Vector(0,(cos(frame*timestep*6)*6),0));
//	slicer1->SetRot(Quaternion(1,0,0,0));
//
//
//	slicer2->SetPos(Vector(3.5,(sin(frame*timestep*6+PI)-1),0));
//	slicer2->SetPos_dt(Vector(0,(cos(frame*timestep*6)*6),0));
//	slicer2->SetRot(Quaternion(1,0,0,0));
	ang += CH_C_PI * timestep / 2.0;
	if (ang >= 2 * CH_C_PI) {
		ang = 0;
	}
	Quaternion q1;
	q1.Q_from_AngY(ang);
	spinner->SetPos(Vector(0, container_height - container_size.y + 1.5, 0));
	spinner->SetPos_dt(Vector(0, 0, 0));
	spinner->SetRot(q1);

	spinner->SetWvel_loc(Vector(0, CH_C_PI / 2.0, 0));
}

int main(int argc, char* argv[]) {
	//int threads = 8;

	if (argc > 1) {
		cohesion = atof(argv[1]);
		data_folder = argv[2];
		//threads = (atoi(argv[1]));
	}
	//omp_set_num_threads(threads);
//=========================================================================================================
	ChSystemParallel * system_gpu = new ChSystemParallel;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

//=========================================================================================================
	//system_gpu->SetParallelThreadNumber(threads);
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	//((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iteration);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(30);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(15);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(0);
	system_gpu->SetTol(.001);
	system_gpu->SetTolSpeeds(.001);
	system_gpu->SetMaxPenetrationRecoverySpeed(25);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(.001);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(1e-4, 1e-4, .2);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(30);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .01);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(30, 30, 30));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
	((ChSystemParallel*) system_gpu)->SetAABB(R3(-6, -6, -6), R3(6, 6, 6));
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Tube = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(.4);
	material->SetRollingFriction(0);
	material->SetSpinningFriction(0);
	material->SetCompliance(0);
	material->SetCohesion(-100);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Top, 100000, Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Tube, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	AddCollisionGeometry(Tube, BOX, Vector(2, container_thickness / 6.0, 1), Vector(0, container_size.y / 2.0 + .6 + .4, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Tube, BOX, Vector(2, container_thickness / 6.0, 1), Vector(0, container_size.y / 2.0 - .6 + .4, 0), Quaternion(1, 0, 0, 0));

	AddCollisionGeometry(Tube, BOX, Vector(2, .6, container_thickness / 6.0), Vector(0, container_size.y / 2.0 + .4, -1), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Tube, BOX, Vector(2, .6, container_thickness / 6.0), Vector(0, container_size.y / 2.0 + .4, 1), Quaternion(1, 0, 0, 0));

	//FinalizeObject(L, (ChSystemParallel *) system_gpu);
	//FinalizeObject(R, (ChSystemParallel *) system_gpu);
	//FinalizeObject(F, (ChSystemParallel *) system_gpu);
	//FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);
	//FinalizeObject(Top, (ChSystemParallel *) system_gpu);
	//FinalizeObject(Tube, (ChSystemParallel *) system_gpu);

	ChSharedBodyPtr Ring;

	for (int i = 0; i < 360; i += 5) {

		real angle = i * PI / 180.0;
		real x = cos(angle) * 5;
		real z = sin(angle) * 5;
		Quaternion q;
		q.Q_from_AngAxis(-angle, Vector(0, 1, 0));

		Ring = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
		InitObject(Ring, 100000, Vector(x, container_height - 3, z), q, material, true, true, -20, -20);

		AddCollisionGeometry(Ring, BOX, Vector(.25, 3, .25), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

		FinalizeObject(Ring, (ChSystemParallel *) system_gpu);
	}

//	slicer1= ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
//	InitObject(slicer1, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, false, -20, -20);
//	AddCollisionGeometry(slicer1, BOX, Vector(container_thickness/15.0, .1, 2), Vector(0,  container_size.y/2.0+.6+.4, 0), Quaternion(1, 0, 0, 0));
//	FinalizeObject(slicer1, (ChSystemParallel *) system_gpu);
//
//	slicer2= ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
//		InitObject(slicer2, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, false, -20, -20);
//		AddCollisionGeometry(slicer2, BOX, Vector(container_thickness/15.0, .1, 2), Vector(0,  container_size.y/2.0+.6+.4, 0), Quaternion(1, 0, 0, 0));
//		FinalizeObject(slicer2, (ChSystemParallel *) system_gpu);
	spinner = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(spinner, 100000, Vector(0, container_height - container_size.y + 1.5, 0), Quaternion(1, 0, 0, 0), material, true, false, -20, -20);
	real spinner_h = .5;
	AddCollisionGeometry(spinner, CYLINDER, Vector(.5, .4, .5), Vector(0, 0, 0), chrono::Q_from_AngAxis(0, ChVector<>(0, 0, 1)));
	AddCollisionGeometry(spinner, BOX, Vector(container_thickness / 15.0, spinner_h, 1.5), Vector(0, 0, 1.75), chrono::Q_from_AngAxis(CH_C_PI / 4.0, ChVector<>(0, 0, 1)));
	AddCollisionGeometry(spinner, BOX, Vector(container_thickness / 15.0, spinner_h, 1.5), Vector(0, 0, -1.75), chrono::Q_from_AngAxis(-CH_C_PI / 4.0, ChVector<>(0, 0, 1)));

	AddCollisionGeometry(spinner, BOX, Vector(1.5, spinner_h, container_thickness / 15.0), Vector(1.75, 0, 0), chrono::Q_from_AngAxis(CH_C_PI / 4.0, ChVector<>(1, 0, 0)));
	AddCollisionGeometry(spinner, BOX, Vector(1.5, spinner_h, container_thickness / 15.0), Vector(-1.75, 0, 0), chrono::Q_from_AngAxis(-CH_C_PI / 4.0, ChVector<>(1, 0, 0)));

	FinalizeObject(spinner, (ChSystemParallel *) system_gpu);

	layer_gen = new ParticleGenerator((ChSystemParallel *) system_gpu);
	layer_gen->SetDensity(1000);
	layer_gen->SetRadius(R3(particle_radius));

	layer_gen->material->SetFriction(.5);
	layer_gen->material->SetCohesion(cohesion);
	layer_gen->material->SetCompliance(1e-4);
	layer_gen->material->SetSpinningFriction(0);
	layer_gen->material->SetRollingFriction(0);
	layer_gen->SetRadius(R3(particle_radius));
	layer_gen->SetCylinderRadius(4.5);
	layer_gen->SetNormalDistribution(particle_radius, .01);
	layer_gen->AddMixtureType(MIX_SPHERE);
	layer_gen->AddMixtureType(MIX_ELLIPSOID);
	layer_gen->AddMixtureType(MIX_CUBE);
	layer_gen->AddMixtureType(MIX_CYLINDER);
	//layer_gen->AddMixtureType(MIX_CONE);

	material_fiber = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_fiber->SetFriction(.4);
	material_fiber->SetRollingFriction(1);
	material_fiber->SetSpinningFriction(1);
	material_fiber->SetCompliance(0);
	material_fiber->SetCohesion(0);

//=========================================================================================================
//Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//
//	//openGLView.render_camera->camera_pos = Vector(0, -5, -10);
//	//	openGLView.render_camera->look_at = Vector(0, -5, 0);
//	//	openGLView.render_camera->mScale = .1;
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
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";
			//DumpAllObjects(system_gpu, ss.str(), ",", true);
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
