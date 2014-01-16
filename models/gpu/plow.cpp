#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"

real gravity = -9.80665;
real timestep = .0001;
real seconds_to_simulate = 30;
int num_steps = seconds_to_simulate / timestep;
int max_iter = 100;

real3 container_size = R3(6.25, 3, 3);
real container_thickness = .2;
real container_height = 0;
real container_friction = 1;

string data_folder = "data/plow";

ChSystemParallel * system_gpu;

ChSharedPtr<ChMaterialSurface> material_shoes, material_chassis;
ChSharedBodyPtr chassis, engine1, engine2, idler;
ChSharedBodyPtr rollers[6];
ChSharedPtr<ChLinkEngine> eng_roller1, eng_roller2;

vector<ChSharedPtr<ChLinkLockRevolute> > revolutes;
vector<ChBody*> shoes;

double idlerPos = 0;
int sim_type = 0;
float mass_multiplier = 10;
float mass_shoe = 2*mass_multiplier;
float mass_idler = 25*mass_multiplier;
float mass_chasis = 500*mass_multiplier;
float mass_sprocket = 25*mass_multiplier;
float mass_roller = 25*mass_multiplier;

float scale_tank = 1;

double L = .09075;
double R = .039;
double H1 = .08;
double H2 = .0739;
double D = .1692;
double W = .015;
int roller_sprocker_counter = 0;
real particle_radius = .02;
ParticleGenerator<ChSystemParallel>* layer_gen;
ChSharedBodyPtr createTrackShoeM113(ChVector<> position, ChQuaternion<> rotation) {
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mrigidBody, mass_shoe, position * scale_tank, rotation, material_shoes, true, false, 3, 3);

	AddCollisionGeometry(mrigidBody, SPHERE, ChVector<>(R, D * .2, R) * scale_tank, ChVector<>(L, 0, 0) * scale_tank, chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	AddCollisionGeometry(mrigidBody, BOX, ChVector<>((L + 2 * R) * .5, W, D) * scale_tank, ChVector<>(0, -H2 + W, 0) * scale_tank, Quaternion(1, 0, 0, 0));
	real rx = (L + 2 * R) * .5;
	real ry = W;
	real rz = D;

	AddCollisionGeometry(mrigidBody, SPHERE, Vector(.01, .01, .01), Vector(-rx, ry - H2 - W, -rz), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(mrigidBody, SPHERE, Vector(.01, .01, .01), Vector(-rx, ry - H2 - W, rz), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(mrigidBody, SPHERE, Vector(.01, .01, .01), Vector(rx, ry - H2 - W, rz), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(mrigidBody, SPHERE, Vector(.01, .01, .01), Vector(rx, ry - H2 - W, -rz), Quaternion(1, 0, 0, 0));

	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);
	//mrigidBody->SetInertiaXX(ChVector<>(.00067, .00082, .00017));
	return mrigidBody;
}

ChSharedBodyPtr createChassisM113(ChVector<> &position) {
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mrigidBody, mass_chasis, position * scale_tank, Quaternion(1, 0, 0, 0), material_chassis, true, false, 4, 3);
	AddCollisionGeometry(mrigidBody, BOX, ChVector<>(1.5, .5, .6) * scale_tank, ChVector<>(0, 0, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);

	ChSharedBodyPtr mPlow = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mPlow, 10*mass_multiplier, position + ChVector<>(-2.9, -.15, 0), Quaternion(1, 0, 0, 0), material_chassis, true, false, 4, 3);
	AddCollisionGeometryTriangleMesh(mPlow, "plow.obj", Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	//AddCollisionGeometry(mPlow, BOX, ChVector<>(.1, .7, 2.2)*scale_tank, ChVector<>(0, 0, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(mPlow, (ChSystemParallel *) system_gpu);
	real rx = .1;
	real ry = .7;
	real rz = 1.65;

	mPlow->SetInertiaXX(ChVector<>(1 / 12.0 * mass_chasis * (ry * ry + rz * rz), 1 / 12.0 * mass_chasis * (rx * rx + rz * rz), 1 / 12.0 * mass_chasis * (rx * rx + ry * ry)));

	ChSharedPtr<ChLinkLockLock> fixedLink(new ChLinkLockLock);
	fixedLink->Initialize(mrigidBody, mPlow, ChCoordsys<>(position + ChVector<>(-2.75, -.2, 0), QUNIT));
	system_gpu->AddLink(fixedLink);

	return mrigidBody;
}
ChSharedBodyPtr createRollerM113(ChVector<> position, ChSharedBodyPtr mrigidBody1, double springLength, double springK, double springR, double rad) {
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mrigidBody, mass_roller, position * scale_tank, chrono::Q_from_AngAxis(0, VECT_X),material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(rad, .02, rad) * scale_tank, ChVector<>(0, 0, -R - .06) * scale_tank, chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(rad, .02, rad) * scale_tank, ChVector<>(0, 0, R + .06) * scale_tank, chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);
	ChSharedBodyPtr ptr1 = ChSharedBodyPtr(mrigidBody);
	ChSharedBodyPtr ptr2 = ChSharedBodyPtr(mrigidBody1);
	ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
	revJoint->Initialize(ptr1, ptr2, ChCoordsys<>(position * scale_tank, QUNIT));
	system_gpu->AddLink(revJoint);
	return mrigidBody;
}
ChSharedBodyPtr createRollerSprocketM113(ChVector<> position, ChSharedBodyPtr mrigidBody1, double springLength, double springK, double springR, double rad) {
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mrigidBody, mass_roller, position * scale_tank, chrono::Q_from_AngAxis(0, VECT_X),material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(rad, .02, rad) * scale_tank, ChVector<>(0, 0, -R - .06) * scale_tank, chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	for (int i = 0; i < 5; i++) {
		AddCollisionGeometry(mrigidBody, BOX, ChVector<>(rad, .01, .16 * .5) * scale_tank, ChVector<>(0, 0, 0) * scale_tank, chrono::Q_from_AngAxis(i * 36.0 * CH_C_PI / 180.0, VECT_Z));
	}
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(rad, .02, rad) * scale_tank, ChVector<>(0, 0, R + .06) * scale_tank, chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);
	ChSharedBodyPtr ptr1 = ChSharedBodyPtr(mrigidBody);
	ChSharedBodyPtr ptr2 = ChSharedBodyPtr(mrigidBody1);
	ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
	revJoint->Initialize(ptr1, ptr2, ChCoordsys<>(position * scale_tank, QUNIT));
	system_gpu->AddLink(revJoint);
	return mrigidBody;
}

ChSharedBodyPtr createSprocketM113(ChVector<> &position, ChSharedBodyPtr mrigidBody1) {
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mrigidBody, mass_sprocket, position * scale_tank, chrono::Q_from_AngAxis(0, VECT_X), material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, ELLIPSOID, ChVector<>(.254, .16 * .5, .254) * scale_tank, ChVector<>(0, 0, 0) * scale_tank, chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	for (int i = 0; i < 5; i++) {
		AddCollisionGeometry(mrigidBody, BOX, ChVector<>(.28965 * 1.1, .03, .16 * .5) * scale_tank, ChVector<>(0, 0, 0) * scale_tank, chrono::Q_from_AngAxis(i * 36.0 * CH_C_PI / 180.0, VECT_Z));
	}
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);
	//mrigidBody->SetInertiaXX(ChVector<>(.3717, .3717, .736));
	ChSharedPtr<ChLinkLockRevolute> my_link1(new ChLinkLockRevolute);
	my_link1->Initialize(mrigidBody, mrigidBody1, ChCoordsys<>(position * scale_tank, QUNIT));
	system_gpu->AddLink(my_link1);
	return mrigidBody;
}

ChSharedBodyPtr createIdlerM113(ChVector<> position, ChSharedBodyPtr mrigidBody1, double springLength, double springK, double springR) {

	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mrigidBody, mass_idler, position * scale_tank, chrono::Q_from_AngAxis(0, VECT_X), material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(.305 / 1.5, .02, .305 / 1.5) * scale_tank, ChVector<>(0, 0, -R - .06) * scale_tank, chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(.305 / 1.5, .02, .305 / 1.5) * scale_tank, ChVector<>(0, 0, R + .06) * scale_tank, chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);
	ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
	revJoint->Initialize(mrigidBody, mrigidBody1, ChCoordsys<>(position * scale_tank, QUNIT));
	system_gpu->AddLink(revJoint);
	return mrigidBody;
}

ChSharedPtr<ChLinkLockLock> createRollerM113WithSpring(ChVector<> position, ChSharedBodyPtr mrigidBody1, double springLength, double springK, double springR, double rad) {
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mrigidBody, mass_roller, position, chrono::Q_from_AngAxis(0, VECT_X), material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(rad, .16 * .3, rad), ChVector<>(0, 0, -.16 * .3 * 2), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(rad, .16 * .3, rad), ChVector<>(0, 0, .16 * .3 * 2), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);

	ChSharedPtr<ChLinkLockLock> fixedLink(new ChLinkLockLock);
	fixedLink->Initialize(mrigidBody, mrigidBody1, ChCoordsys<>(position, QUNIT));
	system_gpu->AddLink(fixedLink);
	return fixedLink;
}

ChSharedBodyPtr createIdlerM113WithSpring(ChVector<> position, ChSharedBodyPtr mrigidBody1, double springLength, double springK, double springR) {
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(mrigidBody, mass_idler, position, chrono::Q_from_AngAxis(0, VECT_X), material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(.305 / 1.5, .16 * .3, .305 / 1.5), ChVector<>(0, 0, -.16 * .3 * 2), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(.305 / 1.5, .16 * .3, .305 / 1.5), ChVector<>(0, 0, .16 * .3 * 2), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);

	ChSharedPtr<ChLinkLockPointLine> revJoint(new ChLinkLockPointLine);
	revJoint->Initialize(mrigidBody, mrigidBody1, ChCoordsys<>(position, QUNIT));
	system_gpu->AddLink(revJoint);

	ChSharedPtr<ChLinkSpring> link_spring = ChSharedPtr<ChLinkSpring>(new ChLinkSpring);
	link_spring->Initialize(mrigidBody, mrigidBody1, false, position, position);
	link_spring->Set_SpringRestLenght(springLength);
	link_spring->Set_SpringK(springK);
	link_spring->Set_SpringR(springR);
	system_gpu->AddLink(link_spring);

	return mrigidBody;
}

ChSharedBodyPtr inputTrackModelM113(int nTracks, real3 pinLoc1, ChVector<> &position, ChSharedBodyPtr chassisBody) {
	ChSharedBodyPtr mrigidBody;
	ChSharedBodyPtr mrigidBodyPrev;
	ChSharedBodyPtr mrigidBodyFirst;

	int bodyNum;
	double pos_x = 0, pos_y = 0, pos_z = 0, rot_x = 0, rot_y = 0, rot_z = 0;

	string temp_data;
	ifstream ifile("pos_M113_4.dat");
	getline(ifile, temp_data);
	stringstream ss(temp_data);

	for (int i = 0; i < nTracks; i++) {
		getline(ifile, temp_data);
		stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;

		mrigidBodyPrev = mrigidBody;
		ChVector<> pos_temp = ChVector<>(pos_x, pos_y, 0) + position;
		ChQuaternion<> temp_Q = chrono::Q_from_AngAxis(rot_z, VECT_Z);
		mrigidBody = createTrackShoeM113(pos_temp, temp_Q);
		shoes.push_back(mrigidBody.get_ptr());
		//mrigidBody->GetCollisionModel()->SetFamily(3);
		//mrigidBody->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);

		if (i == 0) {
			mrigidBodyFirst = mrigidBody;
		}
		if (i != 0) {
			ChCoordsys<> pos1, pos2;
			pos1.pos = Vector(0, 0, 0);
			pos2.pos = Vector(0, 0, 0);
			ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
			revJoint->Initialize(mrigidBody, mrigidBodyPrev, ChCoordsys<>(mrigidBody->Point_Body2World(ChVector<>(pinLoc1.x, pinLoc1.y, pinLoc1.z)), QUNIT));
			system_gpu->AddLink(revJoint);
			revolutes.push_back(revJoint);
		}
		if (i == (nTracks - 1)) {
			ChCoordsys<> pos1, pos2;
			pos1.pos = Vector(0, 0, 0);
			pos2.pos = Vector(0, 0, 0);

			ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
			revJoint->Initialize(mrigidBody, mrigidBodyFirst, ChCoordsys<>(mrigidBodyFirst->Point_Body2World(ChVector<>(pinLoc1.x, pinLoc1.y, pinLoc1.z)), QUNIT));
			system_gpu->AddLink(revJoint);
			revolutes.push_back(revJoint);
		}
	}
	{
		getline(ifile, temp_data);
		std::stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;
		ChVector<> pos_temp = ChVector<>(pos_x, pos_y, 0) + position;
		mrigidBody = createRollerM113(pos_temp, chassisBody, 0, 0, 0, .305 / 1.5);

	}
	for (int i = 0; i < 3; i++) {
		getline(ifile, temp_data);
		std::stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;
		ChVector<> pos_temp = ChVector<>(pos_x, pos_y, 0) + position;
		rollers[roller_sprocker_counter] = createRollerM113(pos_temp, chassisBody, 0, 0, 0, .305 * 1.2);
		ChQuaternion<> temp;
		temp.Q_from_NasaAngles(ChVector<>(rot_x, rot_y, rot_z));
		rollers[roller_sprocker_counter]->SetRot(temp);
		roller_sprocker_counter++;
	}

	ChSharedBodyPtr my_link1;
	for (int i = 0; i < 1; i++) {
		getline(ifile, temp_data);
		std::stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;
		ChVector<> pos_temp = ChVector<>(pos_x, pos_y, 0) + position;
		my_link1 = createSprocketM113(pos_temp, chassisBody);
		ChQuaternion<> temp;
		temp.Q_from_NasaAngles(ChVector<>(rot_x, rot_y, rot_z));
		my_link1->SetRot(temp);
	}
	for (int i = 0; i < 1; i++) {
		getline(ifile, temp_data);
		std::stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;
		ChVector<> pos_temp = ChVector<>(pos_x - .01, pos_y, 0) + position;
		ChSharedBodyPtr idler = createIdlerM113(pos_temp, chassisBody, 10, 500, 100);
	}

	return my_link1;
}

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	engine1->Empty_forces_accumulators();
	engine2->Empty_forces_accumulators();
	if (frame * timestep > .2) {
		engine1->SetWvel_loc(ChVector<>(0, 0, 4));
		engine2->SetWvel_loc(ChVector<>(0, 0, 4));
		for (int i = 0; i < 6; i++) {
			rollers[i]->SetWvel_loc(ChVector<>(0, 0, 4));
		}
		//engine1->Accumulate_torque(ChVector<>(0, 0, 400), 1);
		//engine2->Accumulate_torque(ChVector<>(0, 0, 400), 1);
//		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_roller1->Get_spe_funct())) {
//			mfun->Set_yconst(2);     // rad/s  angular speed
//		}
//
//		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(eng_roller2->Get_spe_funct())) {
//			mfun->Set_yconst(2);     // rad/s  angular speed
//		}
	}

	//Vector pos = chassis->GetPos();

	//((ChSystemParallel*) mSys)->SetAABB(R3(pos.x - 4, pos.y - 3, pos.z - 2), R3(pos.x + 2, pos.y + 3, pos.z + 2));

}

int main(int argc, char* argv[]) {
	if (argc > 1) {
		data_folder = argv[1];

	}
	system_gpu = new ChSystemParallel;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(20);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(40);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(40);
	system_gpu->SetTol(1e-5);
	system_gpu->SetTolSpeeds(1e-5);
	system_gpu->SetMaxPenetrationRecoverySpeed(100);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(1e-5);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(.001, .001, .2);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(100);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(80, 20, 20));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(200, 100);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);

	material->SetCompliance(0);
	material->SetCohesion(0);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);

	material_shoes = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_shoes->SetFriction(1);
	material_shoes->SetCompliance(0);
	material_shoes->SetCohesion(-1000);
	material_chassis = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_chassis->SetFriction(0);
	material_chassis->SetCohesion(-1000);

	real height = -2;
	real x_offset = 1.25;
	ChVector<> temporary1 = ChVector<>(2.0 + x_offset, .1 + height, 0);
	ChVector<> temporary2 = ChVector<>(x_offset, height, .8);
	ChVector<> temporary3 = ChVector<>(x_offset, height, -.8);

	chassis = createChassisM113(temporary1);
	engine1 = inputTrackModelM113(50, R3(-.09075, 0, 0) * scale_tank, temporary2, chassis);
	engine2 = inputTrackModelM113(50, R3(-.09075, 0, 0) * scale_tank, temporary3, chassis);

//	eng_roller1 = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
//	eng_roller1->Initialize(engine1, chassis, ChCoordsys<>(engine1->GetPos(), QUNIT));
//	eng_roller1->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);     // also works as revolute support
//	eng_roller1->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
//	system_gpu->AddLink(eng_roller1);
//
//	eng_roller2 = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
//	eng_roller2->Initialize(engine2, chassis, ChCoordsys<>(engine2->GetPos(), QUNIT));
//	eng_roller2->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);     // also works as revolute support
//	eng_roller2->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
//	system_gpu->AddLink(eng_roller2);

	int3 num_per_dir;
	//num_per_dir = I3(200, 1, 80);
	//num_per_dir = I3(150, 1, 50);
	//num_per_dir = I3(150, 8, 50);
	num_per_dir = I3(50 * 2, 16, 50 * 2);
	layer_gen = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
	layer_gen->SetDensity(10*200);
	layer_gen->SetRadius(R3(particle_radius, particle_radius, particle_radius));
	layer_gen->material->SetFriction(.2);
	layer_gen->material->SetCohesion(.001);
	layer_gen->material->SetRollingFriction(0);
	layer_gen->material->SetSpinningFriction(0);
	layer_gen->material->SetCompliance(.0001);
	layer_gen->AddMixtureType(MIX_SPHERE);
//	layer_gen->AddMixtureType(MIX_ELLIPSOID);
//	layer_gen->AddMixtureType(MIX_DOUBLESPHERE);
//	layer_gen->AddMixtureType(MIX_CUBE);
//	layer_gen->AddMixtureType(MIX_CYLINDER);

	//layer_gen.SetNormalDistribution(rad.x, rad.x/4.0);
	//layer_gen->UseNormalCohesion(particle_cohesion, 1);

	layer_gen->addPerturbedVolumeMixture(R3(-3.25, -2.2, 0), I3(num_per_dir.x, num_per_dir.y, num_per_dir.z), R3(.01, .01, .01), R3(0, 0, 0));

//=========================================================================================================
//Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
////openGLView.render_camera->camera_pos = Vector(0, -5, -10);
////	openGLView.render_camera->look_at = Vector(0, -5, 0);
//	openGLView.render_camera->mScale = .075;
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
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}
}
