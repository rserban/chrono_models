#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 6;

int max_iter = 20;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(1.75, 1, 4.7);
real container_thickness = .04;
real container_height = -1;
Vector container_pos = Vector(0, container_height, 2.8);
real container_friction = 1;

real particle_radius = .02;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = 0;
Vector particle_initial_vel = Vector(0, -5.5, 0); //initial velocity

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

	//AddCollisionGeometry(body, ELLIPSOID, Vector(.35, .1, .35), Vector(0, .05, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(body, ELLIPSOID, Vector(.35, .1, .35), Vector(0, -.05, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometryTriangleMesh(body, "wheel_low_scaled.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	body->GetMaterialSurface()->SetCohesion(-50);
}

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	if (stream) {
		if (((ChSystemGPU*) mSys)->GetNbodies() < 100000) {

			ChSharedBodyPtr sphere;
			real3 rad = R3(particle_radius, particle_radius, particle_radius);
			real3 size = container_size;
			size.y = container_size.y / 3.0;

			int3 num_per_dir = I3(30, 20, 1);

			if (frame % 80 == 0) {
				//addPerturbedLayer(R3(-2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.x, friction.x, cohesion.x, 1e-2, R3(0, 5, 0), (ChSystemGPU*) mSys);
				addPerturbedLayer(
						R3(0, -1, 10),
						SPHERE,
						rad,
						num_per_dir,
						R3(1, 0, 1),
						mass.y,
						0,
						cohesion.y,
						1e-3,
						R3(0, 0, -.5),
						(ChSystemGPU*) mSys);

				addPerturbedLayer(
						R3(0, -1, 6),
						SPHERE,
						rad,
						num_per_dir,
						R3(1, 0, 1),
						mass.y,
						0,
						cohesion.y,
						1e-3,
						R3(0, 0, .5),
						(ChSystemGPU*) mSys);
				//addPerturbedLayer(R3(2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.z, friction.z, cohesion.z, 1e-2, R3(0, 5, 0), (ChSystemGPU*) mSys);
			}
		}
	}

	if (frame > 500) {
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_F->Get_spe_funct())) {
			mfun->Set_yconst(4); // rad/s  angular speed
		}
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_R->Get_spe_funct())) {
			mfun->Set_yconst(4); // rad/s  angular speed
		}
	}

//	if (frame > 5000) {
//		stringstream ss;
//		ss << "data/fording/" << "/" << read_file << ".txt";
//
//		ifstream ifile(ss.str().c_str());
//		string temp;
//		for (int i = 0; i < 14; i++) {
//			getline(ifile, temp);
//
//			std::replace(temp.begin(), temp.end(), ',', '\t');
//			//cout<<temp<<endl;;
//			stringstream st(temp);
//			Vector pos, vel, omg;
//			Quaternion rot;
//
//			st >> pos.x >> pos.y >> pos.z >> rot.e0 >> rot.e1 >> rot.e2 >> rot.e3 >> vel.x >> vel.y >> vel.z >> omg.x >> omg.y >> omg.z;
//			if (i == 7) {
//				chassis->SetPos(pos);
//				chassis->SetRot(rot);
//				chassis->SetPos_dt(vel);
//				chassis->SetWvel_loc(omg);
//			}
//			if (i == 8) {
//				axle_F->SetPos(pos);
//				axle_F->SetRot(rot);
//				axle_F->SetPos_dt(vel);
//				axle_F->SetWvel_loc(omg);
//			}
//			if (i == 9) {
//				axle_R->SetPos(pos);
//				axle_R->SetRot(rot);
//				axle_R->SetPos_dt(vel);
//				axle_R->SetWvel_loc(omg);
//			}
//			if (i == 10) {
//				leg_FR->SetPos(pos);
//				leg_FR->SetRot(rot);
//				leg_FR->SetPos_dt(vel);
//				leg_FR->SetWvel_loc(omg);
//			}
//			if (i == 11) {
//				leg_FL->SetPos(pos);
//				leg_FL->SetRot(rot);
//				leg_FL->SetPos_dt(vel);
//				leg_FL->SetWvel_loc(omg);
//			}
//			if (i == 12) {
//				leg_RR->SetPos(pos);
//				leg_RR->SetRot(rot);
//				leg_RR->SetPos_dt(vel);
//				leg_RR->SetWvel_loc(omg);
//			}
//			if (i == 13) {
//				leg_RL->SetPos(pos);
//				leg_RL->SetRot(rot);
//				leg_RL->SetPos_dt(vel);
//				leg_RL->SetWvel_loc(omg);
//			}
//		}
//		read_file++;
//	}

}

int main(int argc, char* argv[]) {
	omp_set_num_threads(4);
	if (argc == 2) {
		stream = atoi(argv[1]);
	}
	if (argc == 3) {
		stream = atoi(argv[1]);
		int sim = atoi(argv[2]);
		if (sim == 0) {
			mass = R3(.333, .666, .999);
		}
		if (sim == 1) {
			friction = R3(0, .5, 1);
		}
		if (sim == 3) {
			cohesion = R3(0, .5, 1.5);
		}
	}
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
	system_gpu->SetTol(0);
	system_gpu->SetTolSpeeds(0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, .2);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(10);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	mcollisionengine->setBinsPerAxis(R3(100, 100, 100));
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

	InitObject(
			L,
			100000,
			Vector(-container_size.x + container_thickness, -container_thickness, 0) + container_pos,
			Quaternion(1, 0, 0, 0),
			container_friction,
			container_friction,
			0,
			true,
			true,
			-20,
			-20);

	InitObject(
			R,
			100000,
			Vector(container_size.x - container_thickness, -container_thickness, 0) + container_pos,
			Quaternion(1, 0, 0, 0),
			container_friction,
			container_friction,
			0,
			true,
			true,
			-20,
			-20);
	InitObject(
			F,
			100000,
			Vector(0, -container_thickness, -container_size.z + container_thickness) + container_pos,
			Quaternion(1, 0, 0, 0),
			container_friction,
			container_friction,
			0,
			true,
			true,
			-20,
			-20);
	InitObject(
			B,
			100000,
			Vector(0, -container_thickness, container_size.z - container_thickness) + container_pos,
			Quaternion(1, 0, 0, 0),
			container_friction,
			container_friction,
			0,
			true,
			true,
			-20,
			-20);
	InitObject(
			Bottom,
			100000,
			Vector(0, -container_size.y / 3.0, 0) + container_pos,
			Quaternion(1, 0, 0, 0),
			container_friction,
			container_friction,
			0,
			true,
			true,
			-20,
			-20);
	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//AddCollisionGeometryTriangleMesh(Bottom, "ground.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);

	L->GetMaterialSurface()->SetCompliance(0);
	R->GetMaterialSurface()->SetCompliance(0);
	F->GetMaterialSurface()->SetCompliance(0);
	B->GetMaterialSurface()->SetCompliance(0);
	Bottom->GetMaterialSurface()->SetCompliance(0);
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
	num_per_dir.x = (size.x - container_thickness * 2 - rad.x * 2) / rad.x;
	num_per_dir.y = 30;
	num_per_dir.z = (size.z - container_thickness * 2 - rad.z * 2) / rad.z;
	cout << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
	//num_per_dir = I3(10, 10,(size.z-container_thickness*3-rad.z*2) / rad.z);

	real density = 1250;
	real v = 4 / 3.0 * CH_C_PI * pow(particle_radius, 3);
	real mass = density * v;

	cout << "Density " << density << " mass " << mass << " volume " << v << endl;

	//addPerturbedLayer(R3(0, -num_per_dir.y*rad.z+rad.z*2-.5, 0), SPHERE, rad, num_per_dir, R3(0,0,0), mass, .1, .01, R3(0, 0, 0), system_gpu);

	//	impactor = ChSharedBodyGPUPtr(new ChBody(new ChCollisionModelGPU));
	//	InitObject(impactor, 1500, Vector(-container_size.x,container_height + container_size.y*2,0), Quaternion(1, 0, 0, 0), 1, 1, 0, true, false, -1, -2);
	//	AddCollisionGeometry(impactor, SPHERE, ChVector<>(.5,0,0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//	FinalizeObject(impactor, (ChSystemGPU *) system_gpu);
	//	impactor->SetPos_dt(Vector(2.5,0,0));

	//	ChSharedBodyGPUPtr Bunny = ChSharedBodyGPUPtr(new ChBody(new ChCollisionModelGPU));
	////
	//	InitObject(Bunny, 1, Vector(0, 5, 50), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	////	//AddCollisionGeometry(Bunny, BOX, Vector(1, 1, 1),  Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//	AddCollisionGeometryTriangleMesh(Bunny, "trough.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//	FinalizeObject(Bunny, (ChSystemGPU *) system_gpu);
	//
	//	real mass = 1;
	//	Vector r = ChVector<>(1, 1, 1);
	//	Vector inertia = Vector((1 / 5.0 * mass * (r.y * r.y + r.z * r.z), 1 / 5.0 * mass * (r.x * r.x + r.z * r.z), 1 / 5.0 * mass * (r.x * r.x + r.y * r.y)));
	//
	//	Bunny->SetInertiaXX(inertia);

	real offsety = -.4;

	chassis = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	axle_F = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	axle_R = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	leg_FR = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	leg_FL = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	leg_RR = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	leg_RL = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

	InitObject(chassis, 2500 / 4.0, ChVector<>(0, 0, 0), Quaternion(1, 0, 0, 0), 0, 0, 0, false, false, 0, 1);
	InitObject(axle_F, 250 / 4.0, ChVector<>(0, offsety, chassisL / 2.0 + .2), Q_from_AngZ(CH_C_PI / 2.0), 0, 0, 0, false, false, -2, -2);
	InitObject(axle_R, 250 / 4.0, ChVector<>(0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), 0, 0, 0, false, false, -2, -2);
	InitObject(
			leg_FR,
			60 / 4.0,
			ChVector<>((axleL + legW) / 2.0, offsety, chassisL / 2.0 + .2),
			Q_from_AngZ(CH_C_PI / 2.0),
			1,
			1,
			0,
			true,
			false,
			2,
			2);
	InitObject(
			leg_FL,
			60 / 4.0,
			ChVector<>(-(axleL + legW) / 2.0, offsety, chassisL / 2.0 + .2),
			Q_from_AngZ(CH_C_PI / 2.0),
			1,
			1,
			0,
			true,
			false,
			2,
			2);
	InitObject(leg_RR, 60 / 4.0, ChVector<>((axleL + legW) / 2.0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_RL, 60 / 4.0, ChVector<>(-(axleL + legW) / 2.0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, true, false, 2, 2);

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

	//	 floor->SetInertiaXX(Vector(1, 1, 1));
	//	 chassis->SetInertiaXX(Vector(1, 1, 1));
	//	 axle_F->SetInertiaXX(Vector(1, 1, 1));
	//	 axle_R->SetInertiaXX(Vector(1, 1, 1));
	//	 leg_FR->SetInertiaXX(Vector(1, 1, 1));
	//	 leg_FL->SetInertiaXX(Vector(1, 1, 1));
	//	 leg_RR->SetInertiaXX(Vector(1, 1, 1));
	//	 leg_RL->SetInertiaXX(Vector(1, 1, 1));
//
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
	eng_F->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
	eng_F->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);

	system_gpu->AddLink(eng_F);

	eng_R = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	eng_R->Initialize(axle_R_ptr, chassis_ptr, ChCoordsys<>(ChVector<>(0, offsety, -chassisL / 2.0), Q_from_AngY(CH_C_PI / 2.0)));
	eng_R->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
	eng_R->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);

	system_gpu->AddLink(eng_R);

//
	if (!stream) {
		real3 rad = R3(particle_radius, particle_radius, particle_radius);
		real3 size = container_size - R3(container_thickness * 2);
		size.y = container_size.y / 3.0;

		int3 num_per_dir = I3(size.x / rad.x * .9, size.y / rad.y * .85, size.z / rad.z * .85);
		cout << num_per_dir.x * num_per_dir.y * num_per_dir.z * 3 << endl;
		//num_per_dir = I3(1, size.y / rad.y * .85, 1);
		num_per_dir = I3(66, 6, 215);
		//addPerturbedLayer(R3(0, -2, 0), SPHERE, rad, num_per_dir, R3(.1, .1, .1), .333, 0, 0, R3(0, 0, 0), system_gpu);

		addPerturbedLayer(
				R3(0 + container_pos.x, container_pos.y, 0 + container_pos.z),
				SPHERE,
				rad,
				num_per_dir,
				R3(.1, .1, .1),
				mass,
				1,
				50,
				1e-3,
				R3(0, 0, 0),
				system_gpu);
		//addPerturbedLayer(R3(0, 2, 0), SPHERE, rad, num_per_dir, R3(.1, .1, .1), .999, 0, 0, R3(0, 0, 0), system_gpu);
	}
//=========================================================================================================
//Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
//	openGLView.render_camera->look_at = Vector(0, -5, 0);
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
		double RESID = ((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverGPU*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

		int save_every = 1.0 / timestep / 60.0; //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << "data/trough/" << "/" << file << ".txt";
			DumpAllObjects(system_gpu, ss.str(), ",", false);
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
