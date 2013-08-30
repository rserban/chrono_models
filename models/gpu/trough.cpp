#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 15;

int max_iter = 20;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(1.75, 1, 4.7);
real container_thickness = .04;
real container_height = -1.1;
Vector container_pos = Vector(0, container_height, 2.8);
real container_friction = 1;

real particle_radius = .02;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = 0;
Vector particle_initial_vel = Vector(0, -5.5, 0);     //initial velocity

int particle_grid_x = 2;
int particle_grid_z = 2;
real start_height = 1;

double chassisL(2.5);
double chassisW(.4);
double legW(1.0);
double legL(1.0);
double footH(0.5);
double footW(1.5);
double footL(2.0);
double axleL(.8);

bool stream = false;

real3 mass = R3(.034, .034, .034);
real3 friction = R3(0, 0, 0);
real3 cohesion = R3(0, 0, 0);

ChSharedBodyPtr chassis;
ChSharedBodyPtr axle_F;
ChSharedBodyPtr axle_R;
ChSharedBodyPtr leg_FR;
ChSharedBodyPtr leg_FL;
ChSharedBodyPtr leg_RR;
ChSharedBodyPtr leg_RL;
int read_file = 0;

ChSharedPtr<ChLinkEngine> eng_F, eng_R;

void createWheel(ChSharedBodyPtr &body) {

	AddCollisionGeometry(body, ELLIPSOID, Vector(.35, .1, .35), Vector(0, .05, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(body, ELLIPSOID, Vector(.35, .1, .35), Vector(0, -.05, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometryTriangleMesh(body, "wheel_low_scaled.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
}

template<class T>
void RunTimeStep(T* mSys, const int frame) {

	double TIME = mSys->GetChTime();
	double STEP = mSys->GetTimerStep();
	double BROD = mSys->GetTimerCollisionBroad();
	double NARR = mSys->GetTimerCollisionNarrow();
	double LCP = mSys->GetTimerLcp();
	double UPDT = mSys->GetTimerUpdate();
	double RESID = ((ChLcpSolverGPU *) (mSys->GetLcpSolverSpeed()))->GetResidual();
	int BODS = ((ChSystemGPU*) mSys)->GetNbodies();
	int CNTC = ((ChSystemGPU*) mSys)->GetNcontacts();
	int REQ_ITS = ((ChLcpSolverGPU*) (mSys->GetLcpSolverSpeed()))->GetTotalIterations();

	printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

	if (frame > 500) {
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_F->Get_spe_funct())) {
			mfun->Set_yconst(1);     // rad/s  angular speed
		}
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_R->Get_spe_funct())) {
			mfun->Set_yconst(1);     // rad/s  angular speed
		}
	}

}

int main(int argc, char* argv[]) {
	omp_set_num_threads(7);

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
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(1e-5);
	system_gpu->SetTolSpeeds(1e-5);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(8);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->DoStabilization(true);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .04);
	mcollisionengine->setBinsPerAxis(R3(175, 100, 400));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedBodyPtr PF = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr PB = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);
	material->SetCompliance(0);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, -container_thickness, 0) + container_pos, Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	InitObject(R, 100000, Vector(container_size.x - container_thickness, -container_thickness, 0) + container_pos, Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, -container_thickness, -container_size.z + container_thickness) + container_pos, Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, -container_thickness, container_size.z - container_thickness) + container_pos, Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height+container_size.y/1.75, 0) + container_pos, Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//AddCollisionGeometryTriangleMesh(Bottom, "ground.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	//ChSharedPtr<ChAsset> asset = Bottom->GetAssets().at(0);
	ChTriangleMeshConnected trimesh;
	trimesh.LoadWavefrontMesh("ground.obj", false, false);
	std::vector<ChVector<double> > verts = trimesh.m_vertices;
	ChSharedBodyPtr sphere;
	for (int i = 0; i < verts.size(); i++) {
		AddCollisionGeometry(Bottom, ELLIPSOID, ChVector<>(.08, .03, .08), verts[i], Quaternion(1, 0, 0, 0));
	}

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);

	real wheel_mass = 60;

	//wheel = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	//InitObject(wheel, wheel_mass, ChVector<>(0, .3, 0), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, true, false, 2, 2);
	//createWheel(wheel);
	//AddCollisionGeometry(wheel, CYLINDER, Vector(.7,.2,.7), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//FinalizeObject(wheel, (ChSystemGPU *) system_gpu);

	//Vector inertia = Vector(1 / 12.0 * wheel_mass * (3 * .7 * .7 + .2 * .2), 1 / 2.0 * wheel_mass * (.7 * .7), 1 / 12.0 * wheel_mass * (3 * .7 * .7 + .2 * .2));

	//wheel->SetInertiaXX(inertia);
	//wheel->GetMaterialSurface()->SetCohesion(-.01);
	real3 rad = R3(particle_radius, particle_radius, particle_radius);
	real3 size = container_size;
	size.y = container_size.y / 3.0;

	int3 num_per_dir;

	real density = 1250;
	real v = 4 / 3.0 * CH_C_PI * pow(particle_radius, 3);
	real mass = density * v;

	cout << "Density " << density << " mass " << mass << " volume " << v << endl;

	real offsety = -.3;

	chassis = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	axle_F = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	axle_R = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	leg_FR = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	leg_FL = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	leg_RR = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	leg_RL = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

	ChSharedPtr<ChMaterialSurface> material_chassis;
	material_chassis = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_chassis->SetFriction(0);
	material_chassis->SetCompliance(0);
	material_chassis->SetCohesion(-100);

	ChSharedPtr<ChMaterialSurface> material_wheel;
	material_wheel = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_wheel->SetFriction(1);
	material_wheel->SetCompliance(0);
	material_wheel->SetCohesion(-100);

	InitObject(chassis, 2500 / 1.0, ChVector<>(0, 0, 0), Quaternion(1, 0, 0, 0), material_chassis, false, false, -2, -2);
	InitObject(axle_F, 250 / 1.0, ChVector<>(0, offsety, chassisL / 2.0 + .2), Q_from_AngZ(CH_C_PI / 2.0), material_chassis, false, false, -2, -2);
	InitObject(axle_R, 250 / 1.0, ChVector<>(0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), material_chassis, false, false, -2, -2);
	InitObject(leg_FR, 60 / 1.0, ChVector<>((axleL + legW) / 2.0, offsety, chassisL / 2.0 + .2), Q_from_AngZ(CH_C_PI / 2.0), material_wheel, true, false, -2, -2);
	InitObject(leg_FL, 60 / 1.0, ChVector<>(-(axleL + legW) / 2.0, offsety, chassisL / 2.0 + .2), Q_from_AngZ(CH_C_PI / 2.0), material_wheel, true, false, -2, -2);
	InitObject(leg_RR, 60 / 1.0, ChVector<>((axleL + legW) / 2.0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), material_wheel, true, false, -2, -2);
	InitObject(leg_RL, 60 / 1.0, ChVector<>(-(axleL + legW) / 2.0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), material_wheel, true, false, -2, -2);

	AddCollisionGeometry(chassis, BOX, ChVector<>(.5, .1, chassisL / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//AddCollisionGeometryTriangleMesh(chassis, "humvee.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	AddCollisionGeometry(axle_F, ELLIPSOID, ChVector<>(0.5 / 2.0, 0.5 / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(axle_R, ELLIPSOID, ChVector<>(0.5 / 2.0, 0.5 / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	createWheel(leg_FR);
	createWheel(leg_FL);
	createWheel(leg_RR);
	createWheel(leg_RL);

	FinalizeObject(chassis, (ChSystemGPU *) system_gpu);
	FinalizeObject(axle_F, (ChSystemGPU *) system_gpu);
	FinalizeObject(axle_R, (ChSystemGPU *) system_gpu);
	FinalizeObject(leg_FR, (ChSystemGPU *) system_gpu);
	FinalizeObject(leg_FL, (ChSystemGPU *) system_gpu);

	FinalizeObject(leg_RR, (ChSystemGPU *) system_gpu);
	FinalizeObject(leg_RL, (ChSystemGPU *) system_gpu);


	ChSharedBodyPtr chassis_ptr = ChSharedBodyPtr(chassis);
	ChSharedBodyPtr axle_F_ptr = ChSharedBodyPtr(axle_F);
	ChSharedBodyPtr axle_R_ptr = ChSharedBodyPtr(axle_R);
	ChSharedBodyPtr leg_FR_ptr = ChSharedBodyPtr(leg_FR);
	ChSharedBodyPtr leg_FL_ptr = ChSharedBodyPtr(leg_FL);
	ChSharedBodyPtr leg_RR_ptr = ChSharedBodyPtr(leg_RR);
	ChSharedBodyPtr leg_RL_ptr = ChSharedBodyPtr(leg_RL);

	//	// attach legs to axles
	ChSharedPtr<ChLinkLockLock> axle_FR(new ChLinkLockLock);
	axle_FR->Initialize(leg_FR_ptr, axle_F_ptr, ChCoordsys<>(VNULL));
	system_gpu->AddLink(axle_FR);

	ChSharedPtr<ChLinkLockLock> axle_FL(new ChLinkLockLock);
	axle_FL->Initialize(leg_FL_ptr, axle_F_ptr, ChCoordsys<>(VNULL));
	system_gpu->AddLink(axle_FL);

	ChSharedPtr<ChLinkLockLock> axle_RR(new ChLinkLockLock);
	axle_RR->Initialize(leg_RR_ptr, axle_R_ptr, ChCoordsys<>(VNULL));
	system_gpu->AddLink(axle_RR);

	ChSharedPtr<ChLinkLockLock> axle_RL(new ChLinkLockLock);
	axle_RL->Initialize(leg_RL_ptr, axle_R_ptr, ChCoordsys<>(VNULL));
	system_gpu->AddLink(axle_RL);

	eng_F = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	eng_F->Initialize(axle_F_ptr, chassis_ptr, ChCoordsys<>(ChVector<>(0, offsety, chassisL / 2.0 + .2), Q_from_AngY(CH_C_PI / 2.0)));
	eng_F->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);     // also works as revolute support
	eng_F->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);

	system_gpu->AddLink(eng_F);

	eng_R = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	eng_R->Initialize(axle_R_ptr, chassis_ptr, ChCoordsys<>(ChVector<>(0, offsety, -chassisL / 2.0), Q_from_AngY(CH_C_PI / 2.0)));
	eng_R->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);     // also works as revolute support
	eng_R->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);

	system_gpu->AddLink(eng_R);

	if (!stream) {
		real3 rad = R3(particle_radius, particle_radius, particle_radius);
		real3 size = container_size - R3(container_thickness * 2);
		size.y = container_size.y / 3.0;

		int3 num_per_dir = I3(size.x / rad.x * .9, size.y / rad.y * .85, size.z / rad.z * .85);

		//num_per_dir = I3(1, size.y / rad.y * .85, 1);
		num_per_dir = I3(66, 16, 215);
		//num_per_dir = I3(10, 12, 10);
		cout << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
		//addPerturbedLayer(R3(0, -2, 0), SPHERE, rad, num_per_dir, R3(.1, .1, .1), .333, 0, 0, R3(0, 0, 0), system_gpu);

		ParticleGenerator layer_gen(system_gpu);
		layer_gen.SetDensity(1500);
		layer_gen.SetRadius(R3(particle_radius));

		layer_gen.material->SetFriction(1);
		layer_gen.material->SetCohesion(50);
		//layer_gen.addPerturbedVolume(R3(0 + container_pos.x, container_pos.y, 0 + container_pos.z), SPHERE, I3(num_per_dir.x, 1, num_per_dir.z), R3(.1, .1, .1), R3(0, 0, 0), false);
		layer_gen.SetNormalDistribution(particle_radius - particle_radius / 6.0, particle_radius / 6.0);
		layer_gen.addPerturbedVolume(R3(0 + container_pos.x, container_pos.y, 0 + container_pos.z), SPHERE, num_per_dir, R3(.1, .1, .1), R3(0, 0, 0), false);



		//layer_gen.SetRadius(R3(particle_radius * 4));
		//.SetNormalDistribution(particle_radius * 2, particle_radius / 2.0);
		//layer_gen.addPerturbedVolume(R3(0 + container_pos.x, container_pos.y + particle_radius * 25, 0 + container_pos.z + 2), BOX, I3(13, 1, 20), R3(1, 1, 1), R3(0, 0, 0), false);
		//addPerturbedLayer(R3(0 + container_pos.x, container_pos.y, 0 + container_pos.z), SPHERE, rad, num_per_dir, R3(.1, .1, .1), mass, 1, 25, 0, R3(0, 0, 0), system_gpu);
		//addPerturbedLayer(R3(0, 2, 0), SPHERE, rad, num_per_dir, R3(.1, .1, .1), .999, 0, 0, R3(0, 0, 0), system_gpu);
	}
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
			ss << "data/trough/" << "/" << file << ".txt";
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
