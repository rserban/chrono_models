#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 9;

int max_iter = 40;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(1.5, 1, 3);
real container_thickness = .04;
real container_height = -1.5;
Vector container_pos = Vector(0, container_height, 8);
real container_friction = 1;

real particle_radius = .02;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = .1;
Vector particle_initial_vel = Vector(0, -5.5, 0); //initial velocity

int particle_grid_x = 14;
int particle_grid_z = 14;
real start_height = 1;

double chassisL(2.5);
double chassisW(.4);
double legW(1.0);
double legL(1.0);
double footH(0.5);
double footW(1.5);
double footL(2.0);
double axleL(.8);

ChSharedBodyPtr wheel;
real ang = 0;
template<class T>
void RunTimeStep(T* mSys, const int frame) {

	double TIME = mSys->GetChTime();
	double STEP = mSys->GetTimerStep();
	double BROD = mSys->GetTimerCollisionBroad();
	double NARR = mSys->GetTimerCollisionNarrow();
	double LCP = mSys->GetTimerLcp();
	double UPDT = mSys->GetTimerUpdate();
	int BODS = mSys->GetNbodies();
	int CNTC = mSys->GetNcontacts();
	int REQ_ITS = ((ChLcpSolverGPU*) (mSys->GetLcpSolverSpeed()))->GetTotalIterations();
	real RESID = ((ChLcpSolverGPU*) (mSys->GetLcpSolverSpeed()))->GetResidual();

	printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

//	ChSharedBodyGPUPtr sphere;
//	if (frame % 50 == 0) {
//		for (int i = 0; i < particle_grid_x; i++) {
//			for (int j = 0; j < particle_grid_z; j++) {
//				sphere = ChSharedBodyGPUPtr(new ChBody(new ChCollisionModelGPU));
//				Quaternion q;
//				q.Q_from_NasaAngles(Vector(rand() % 1000 / 1000.0, rand() % 1000 / 1000.0, rand() % 1000 / 1000.0));
//				q.Normalize();
//
//				ChVector<> position(i * particle_radius * 3 - particle_grid_x * .5, start_height, j * particle_radius * 3 - particle_grid_z * .5);
//
//				//position.x += rand() % 1000 / 100000.0;
//				//position.y += rand() % 1000 / 1000.0;
//				//position.z += rand() % 1000 / 100000.0;
//
//				real percent = 0; //(rand() % 10000 / 10000.0 - .5) / 2.0;
//				real radius = particle_radius + particle_radius * percent;
//				InitObject(sphere, particle_mass, position, Quaternion(1, 0, 0, 0), particle_friction, particle_friction, 0, true, false, -1, i);
//				AddCollisionGeometry(sphere, SPHERE, ChVector<>(radius, radius, radius), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
//				sphere->SetPos_dt(particle_initial_vel);
//				FinalizeObject(sphere, (ChSystemGPU *) mSys);
//			}
//		}
//	}

//	ang += CH_C_PI * timestep*.5;
//	if (ang >= 2 * CH_C_PI) {
//		ang = 0;
//	}
//
//	Quaternion q1;
//	q1 = Q_from_AngZ(CH_C_PI / 2.0);
//	q1 = q1 % Q_from_AngY(ang);
//
//	q1.Normalize();
//	Vector pos = wheel->GetPos();
//	wheel->SetPos(Vector(0, pos.y, pos.z));
//	wheel->SetRot(q1);
//	wheel->SetWvel_loc(Vector(0, CH_C_PI*.5, 0));

}

void createWheel(ChSharedBodyPtr &body) {

//	AddCollisionGeometry(body, ELLIPSOID, ChVector<>(legW / 2.0, .1, legW / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
//
//	AddCollisionGeometry(body, BOX, ChVector<>(.1, legL / 2.0, .1), Vector(0, 0, legW / 2.0), Quaternion(1, 0, 0, 0));
//	AddCollisionGeometry(body, BOX, ChVector<>(.1, legL / 2.0, .1), Vector(0, 0, -legW / 2.0), Quaternion(1, 0, 0, 0));
//	AddCollisionGeometry(body, BOX, ChVector<>(.1, legL / 2.0, .1), Vector(legW / 2.0, 0, 0), Quaternion(1, 0, 0, 0));
//	AddCollisionGeometry(body, BOX, ChVector<>(.1, legL / 2.0, .1), Vector(-legW / 2.0, 0, 0), Quaternion(1, 0, 0, 0));
//
//	AddCollisionGeometry(body, BOX, ChVector<>(.1, legL / 2.0, .1), Vector(legW / 2.0/sqrt(3), 0, legW / 2.0/sqrt(3)), Quaternion(1, 0, 0, 0));
//	AddCollisionGeometry(body, BOX, ChVector<>(.1, legL / 2.0, .1), Vector(legW / 2.0/sqrt(3), 0, -legW / 2.0/sqrt(3)), Quaternion(1, 0, 0, 0));
//
//	AddCollisionGeometry(body, BOX, ChVector<>(.1, legL / 2.0, .1), Vector(-legW / 2.0/sqrt(3), 0, legW / 2.0/sqrt(3)), Quaternion(1, 0, 0, 0));
//	AddCollisionGeometry(body, BOX, ChVector<>(.1, legL / 2.0, .1), Vector(-legW / 2.0/sqrt(3), 0, -legW / 2.0/sqrt(3)), Quaternion(1, 0, 0, 0));

	//	ChSharedBodyGPUPtr Bunny = ChSharedBodyGPUPtr(new ChBody(new ChCollisionModelGPU));
	//
	//	InitObject(Bunny, 1, Vector(-3, 0, 50), Quaternion(1, 0, 0, 0), container_friction, container_friction, 0, true, true, -20, -20);
	AddCollisionGeometry(body, ELLIPSOID, Vector(.4, .1, .4), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//AddCollisionGeometryTriangleMesh(body, "wheel_low.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//	FinalizeObject(Bunny, (ChSystemGPU *) system_gpu);

}

int main(int argc, char* argv[]) {
	omp_set_num_threads(1);
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
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(.5);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	//((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(CONJUGATE_GRADIENT);
	//((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(MINIMUM_RESIDUAL);
	//((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(STEEPEST_DESCENT);
	//BLOCK_JACOBI
	//ACCELERATED_PROJECTED_GRADIENT_DESCENT
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	mcollisionengine->setBinsPerAxis(R3(40, 40, 200));
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

	Quaternion AA, BB;
	AA.Q_from_AngX(-PI / 3.0);
	BB.Q_from_AngX(PI / 3.0);

	InitObject(
			PF,
			100000,
			Vector(0, -container_thickness + container_size.y * .5, -container_size.z * 2.3 + container_thickness) + container_pos,
			Quaternion(1, 0, 0, 0),
			container_friction,
			container_friction,
			0,
			true,
			true,
			-20,
			-20);

	InitObject(
			PB,
			100000,
			Vector(0, -container_thickness + container_size.y * .5, container_size.z * 2.3 + container_thickness) + container_pos,
			Quaternion(1, 0, 0, 0),
			container_friction,
			container_friction,
			0,
			true,
			true,
			-20,
			-20);

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
			AA,
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
			BB,
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
	AddCollisionGeometry(PF, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(PB, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x * 1.5, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x * 1.5, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(PF, (ChSystemGPU *) system_gpu);
	FinalizeObject(PB, (ChSystemGPU *) system_gpu);
	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);

	PF->GetMaterialSurface()->SetCompliance(0);
	PB->GetMaterialSurface()->SetCompliance(0);
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

	ChSharedBodyPtr floor(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr chassis(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr axle_F(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr axle_R(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr leg_FR(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr leg_FL(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr leg_RR(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr leg_RL(new ChBody(new ChCollisionModelGPU));

	//InitObject(floor, 1.0, ChVector<>(0, -3.5, 0), Quaternion(1, 0, 0, 0), 1, 1, 0, true, true, -20, -20);
	InitObject(chassis, 10, ChVector<>(0, 0, 0), Quaternion(1, 0, 0, 0), 0, 0, 0, false, false, 0, 1);
	InitObject(axle_F, 1, ChVector<>(0, offsety, chassisL / 2.0 + .2), Q_from_AngZ(CH_C_PI / 2.0), 0, 0, 0, false, false, -2, -2);
	InitObject(axle_R, 1, ChVector<>(0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), 0, 0, 0, false, false, -2, -2);
	InitObject(leg_FR, 1, ChVector<>((axleL + legW) / 2.0, offsety, chassisL / 2.0 + .2), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_FL, 1, ChVector<>(-(axleL + legW) / 2.0, offsety, chassisL / 2.0 + .2), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_RR, 1, ChVector<>((axleL + legW) / 2.0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, true, false, 2, 2);
	InitObject(leg_RL, 1, ChVector<>(-(axleL + legW) / 2.0, offsety, -chassisL / 2.0), Q_from_AngZ(CH_C_PI / 2.0), 1, 1, 0, true, false, 2, 2);

	//AddCollisionGeometry(floor, BOX, ChVector<>(10, 1 / 2.0, 50), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(chassis, BOX, ChVector<>(.5, .1, chassisL / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometryTriangleMesh(chassis, "humvee.obj", Vector(0, -1, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(axle_F, ELLIPSOID, ChVector<>(0.5 / 2.0, 0.5 / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(axle_R, ELLIPSOID, ChVector<>(0.5 / 2.0, 0.5 / 2.0, 0.5 / 2.0), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	createWheel(leg_FR);
	createWheel(leg_FL);
	createWheel(leg_RR);
	createWheel(leg_RL);

	//FinalizeObject(floor, (ChSystemGPU *) system_gpu);
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

	ChSharedPtr<ChLinkEngine> eng_F(new ChLinkEngine);
	eng_F->Initialize(axle_F_ptr, chassis_ptr, ChCoordsys<>(ChVector<>(0, offsety, chassisL / 2.0 + .2), Q_from_AngY(CH_C_PI / 2.0)));
	eng_F->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
	eng_F->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_F->Get_spe_funct()))
		mfun->Set_yconst(4); // rad/s  angular speed
	system_gpu->AddLink(eng_F);

	ChSharedPtr<ChLinkEngine> eng_R(new ChLinkEngine);
	eng_R->Initialize(axle_R_ptr, chassis_ptr, ChCoordsys<>(ChVector<>(0, offsety, -chassisL / 2.0), Q_from_AngY(CH_C_PI / 2.0)));
	eng_R->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
	eng_R->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_R->Get_spe_funct()))
		mfun->Set_yconst(4); // rad/s  angular speed
	system_gpu->AddLink(eng_R);

	//system_gpu->DoStepDynamics(timestep);
	//system_gpu->DoStepDynamics(timestep);
	//exit(0);
	//=========================================================================================================
	////Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
	openGLView.render_camera->look_at = Vector(0, -5, 0);
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
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverGPU*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS);
		int save_every = 1.0 / timestep / 60.0; //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << "data/trough/" << "/" << file << ".txt";
			DumpAllObjects(system_gpu, ss.str(), ",", true);
			//DumpAllObjectsWithGeometry(system_gpu, ss.str(), ",");
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
