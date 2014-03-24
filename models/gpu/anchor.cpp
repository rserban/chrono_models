#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"

real gravity = -9806.65;
real timestep = .00025;
real seconds_to_simulate = 5;
real tolerance = 0;

int max_iter = 40;
int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(150, 300, 150);
real container_thickness = 10;
real container_height = 0;
real container_friction = .1;
real container_cohesion = -1000;

real particle_radius = 1.5 * 2;
real particle_density = .00265;
real particle_slide_friction = .3;
real particle_roll_friction = .3;
real particle_cohesion = 0;
real particle_std_dev = .5 * 2;

ChSharedBodyPtr BLOCK, CONTAINER, ANCHOR;
ParticleGenerator<ChSystemParallel>* layer_gen;
real amplitude = particle_radius;
real frequency = 10;

real anchor_density = 0.00785;
real anchor_vel = -20;
real anchor_rot = -30;     //rpm
int layers = 0;

string data_folder = "data/anchor/";

ChFunction_Ramp* motionFunc1, *motionFunc2;
ChSharedPtr<ChLinkLockLock> actuator_anchor;
ChSharedPtr<ChLinkEngine> engine_anchor;
bool once = true;
bool save;
template<class T>
void RunTimeStep(T* mSys, const int frame) {

	Vector force = engine_anchor->Get_react_force();
	Vector torque = engine_anchor->Get_react_torque();
	double motor_torque = engine_anchor->Get_mot_torque();
	cout << force.x << " " << force.y << " " << force.z << " " << torque.x << " " << torque.y << " " << torque.z << " " << motor_torque << " " << ANCHOR->GetPos().y << " "
			<< ANCHOR->GetPos_dt().y << endl;

//	real t = frame * timestep * PI * 2 * frequency;
//
//	BLOCK->SetRot(ChQuaternion<>(1, 0, 0, 0));
//	BLOCK->SetWvel_loc(ChVector<>(0, 0, 0));
//	BLOCK->SetPos(ChVector<>(sin(t) * amplitude, BLOCK->GetPos().y, 0));
//	BLOCK->SetPos_dt(ChVector<>(cos(t) * amplitude * 2 * PI * frequency, BLOCK->GetPos_dt().y, 0));
//
//	CONTAINER->SetPos(ChVector<>(sin(t) * amplitude, 0, 0));
//	CONTAINER->SetPos_dt(ChVector<>(cos(t) * amplitude * 2 * PI * frequency, 0, 0));
//	CONTAINER->SetWvel_loc(ChVector<>(0, 0, 0));
//	CONTAINER->SetRot(ChQuaternion<>(1, 0, 0, 0));
//
//
//	real cont_vol = (container_size.x - container_thickness * 2) * 2 * (BLOCK->GetPos().y + container_size.y - 2 * container_thickness) * (container_size.z - container_thickness * 2) * 2;
//	cout << layer_gen->total_volume << " " << layer_gen->total_mass << " " << cont_vol << " " << layer_gen->total_mass / cont_vol << endl;
//

	if (save) {
		if (layers < 130 && frame % 60 == 0) {
			cout << layers << endl;
			layer_gen->addPerturbedVolumeMixture(R3(0, -container_size.y + container_thickness + particle_radius * 5 + frame / 14.0, 0), I3(32, 1, 32), R3(0, 0, 0), R3(0, 0, 0));
			layers++;
		}
	} else {
		//300 - 457.2
//		double time = actuator_anchor->GetChTime();
//		if (ANCHOR->GetPos().y <= 300 - 457.2 && once) {
//			motionFunc1->Set_y0(time * -anchor_vel);
//			motionFunc1->Set_ang(-2);
////			motionFunc2->Set_y0(time * -anchor_rot * 1 / 60.0 * 2 * CH_C_PI);
////			motionFunc2->Set_ang(0);
//			if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(engine_anchor->Get_spe_funct())) {
//				mfun->Set_yconst(0);     // rad/s  angular speed
//			}
//			once = false;
//		}
	}

}

int main(int argc, char* argv[]) {
	save = atoi(argv[1]);
	cout << save << endl;
	if (save) {

		seconds_to_simulate = 5;
		num_steps = seconds_to_simulate / timestep;
	} else {
		seconds_to_simulate = 300;
		num_steps = seconds_to_simulate / timestep;
	}

	if (argc > 2) {
		particle_slide_friction = atof(argv[2]);
		particle_roll_friction = atof(argv[3]);
		particle_cohesion = atof(argv[4]);
		data_folder = argv[5];
		cout<<particle_slide_friction<<" "<<particle_roll_friction<<" "<<particle_roll_friction<<endl;
	}

//=========================================================================================================
	ChSystemParallelDVI * system_gpu = new ChSystemParallelDVI;

//=========================================================================================================
	system_gpu->SetMinThreads(32);
	//system_gpu->SetMaxiter(max_iter);
	//system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(max_iter);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(max_iter);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(max_iter);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(max_iter);
	system_gpu->SetTol(tolerance);
	system_gpu->SetTolSpeeds(tolerance);
	system_gpu->SetMaxPenetrationRecoverySpeed(500);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(tolerance);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(500);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(30, 60, 30));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
	((ChSystemParallel*) system_gpu)->SetAABB(R3(-250, -600, -250), R3(250, 600, 250));

//=========================================================================================================

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);
	material->SetRollingFriction(container_friction);
	material->SetSpinningFriction(container_friction);
	material->SetCompliance(0);
	material->SetCohesion(-100);

	CONTAINER = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	InitObject(CONTAINER, 100000, Vector(0, 0, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	if (save) {
		AddCollisionGeometry(
				CONTAINER,
				BOX,
				Vector(container_thickness, container_size.y, container_size.z),
				Vector(-container_size.x + container_thickness, container_height - container_thickness, 0),
				Quaternion(1, 0, 0, 0));
		AddCollisionGeometry(
				CONTAINER,
				BOX,
				Vector(container_thickness, container_size.y, container_size.z),
				Vector(container_size.x - container_thickness, container_height - container_thickness, 0),
				Quaternion(1, 0, 0, 0));
		AddCollisionGeometry(
				CONTAINER,
				BOX,
				Vector(container_size.x, container_size.y, container_thickness),
				Vector(0, container_height - container_thickness, -container_size.z + container_thickness),
				Quaternion(1, 0, 0, 0));
		AddCollisionGeometry(
				CONTAINER,
				BOX,
				Vector(container_size.x, container_size.y, container_thickness),
				Vector(0, container_height - container_thickness, container_size.z - container_thickness),
				Quaternion(1, 0, 0, 0));
		AddCollisionGeometry(
				CONTAINER,
				BOX,
				Vector(container_size.x, container_thickness, container_size.z),
				Vector(0, container_height - container_size.y, 0),
				Quaternion(1, 0, 0, 0));
		//AddCollisionGeometry(CONTAINER, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, container_height + container_size.y, 0), Quaternion(1, 0, 0, 0));
	}
	CONTAINER->GetMaterialSurface()->SetCohesion(container_cohesion);
	FinalizeObject(CONTAINER, (ChSystemParallel *) system_gpu);

	if (save == false) {
		real anchor_length = 100;
		real anchor_r = 35 / 2.0;
		real anchor_R = 150 / 4.0;
		real anchor_h = 50;
		real anchor_thickness = 6 / 2.0;
		real anchor_blade_width = 4;
		ChVector<> p1(0, 0, 0);
		ChVector<> p2(0, anchor_length, 0);
		real anchor_mass = 6208;
		real number_sections = 150;

		ANCHOR = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
		InitObject(ANCHOR, anchor_mass, Vector(0, 100, 0), Quaternion(1, 0, 0, 0), material, true, false, -15, -15);
		AddCollisionGeometry(ANCHOR, SPHERE, ChVector<>(anchor_r, 0, 0), p1, Quaternion(1, 0, 0, 0));
		AddCollisionGeometry(ANCHOR, CYLINDER, Vector(anchor_r, anchor_length, anchor_r), p2, Quaternion(1, 0, 0, 0));
//
		for (int i = 0; i < number_sections; i++) {
			ChQuaternion<> quat, quat2;
			quat.Q_from_AngAxis(i / number_sections * 2 * PI, ChVector<>(0, 1, 0));
			quat2.Q_from_AngAxis(6 * 2 * PI / 360.0, ChVector<>(0, 0, 1));
			quat = quat % quat2;
			ChVector<> pos(sin(i / number_sections * 2 * PI) * anchor_R, i / number_sections * anchor_h, cos(i / number_sections * 2 * PI) * anchor_R);
			//ChMatrix33<> mat(quat);
			AddCollisionGeometry(ANCHOR, BOX, ChVector<>(anchor_blade_width, anchor_thickness, anchor_R), pos, quat);
		}

		FinalizeObject(ANCHOR, (ChSystemParallel *) system_gpu);

//	real vol = ((ChCollisionModelParallel *) ANCHOR->GetCollisionModel())->getVolume();
//	real den = .00785;
//	anchor_mass = den*vol;
//	cout<<vol<<" "<<anchor_mass<<endl;
//	ANCHOR->SetMass(anchor_mass);

		ANCHOR->SetInertiaXX(
				ChVector<>(
						1 / 12.0 * anchor_mass * (1 * 1 + anchor_R * anchor_R),
						1 / 12.0 * anchor_mass * (anchor_R * anchor_R + anchor_R * anchor_R),
						1 / 12.0 * anchor_mass * (anchor_R * anchor_R + 1 * 1)));

		//BLOCK = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
		//InitObject(BLOCK, anchor_mass/2, Vector(0, 300, 0), Quaternion(1, 0, 0, 0), material, false, false, -20, -20);
		//AddCollisionGeometry(BLOCK, BOX, ChVector<>(1, 1, 1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
		//FinalizeObject(BLOCK, (ChSystemParallel *) system_gpu);

		//actuator_anchor = ChSharedPtr<ChLinkLockLock>(new ChLinkLockLock());
		//actuator_anchor->Initialize(CONTAINER, BLOCK, ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT));
		//system_gpu->AddLink(actuator_anchor);

		// apply motion
		//motionFunc1 = new ChFunction_Ramp(0, -anchor_vel);
		//actuator_anchor->SetMotion_Y(motionFunc1);
		//actuator_anchor->SetMotion_axis(ChVector<>(0, 1, 0));

		engine_anchor = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
		engine_anchor->Initialize(CONTAINER, ANCHOR, ChCoordsys<>(ANCHOR->GetPos(), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X)));
		engine_anchor->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_PRISM);     // also works as revolute support
		engine_anchor->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);

		system_gpu->AddLink(engine_anchor);
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(engine_anchor->Get_spe_funct())) {
			mfun->Set_yconst(anchor_rot * 1 / 60.0 * 2 * CH_C_PI);     // rad/s  angular speed
		}
		ReadAllObjectsWithGeometryChrono(system_gpu, "data/anchor/anchor.dat");
		ChSharedPtr<ChMaterialSurface> material_read;
		material_read = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
		material_read->SetFriction(particle_slide_friction);
		material_read->SetRollingFriction(particle_roll_friction);
		material_read->SetSpinningFriction(particle_roll_friction);
		material_read->SetCompliance(0);
		material_read->SetCohesion(particle_cohesion*timestep);

		for (int i = 0; i < system_gpu->Get_bodylist()->size(); i++) {
			system_gpu->Get_bodylist()->at(i)->SetMaterialSurface(material_read);
		}

	} else {
		layer_gen = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
		layer_gen->SetDensity(particle_density);
		layer_gen->SetRadius(R3(particle_radius, particle_radius * .5, particle_radius));
		layer_gen->SetNormalDistribution(particle_radius, particle_std_dev, 1);
		layer_gen->material->SetFriction(particle_slide_friction);
		layer_gen->material->SetCohesion(particle_cohesion);
		layer_gen->material->SetRollingFriction(0);
		layer_gen->material->SetSpinningFriction(0);
		layer_gen->AddMixtureType(MIX_SPHERE);
	}
	//layer_gen->AddMixtureType(MIX_ELLIPSOID);
	//layer_gen->AddMixtureType(MIX_DOUBLESPHERE);
	//layer_gen->addPerturbedVolumeMixture(R3(0, 0, 0), I3(64, 1, 64), R3(0, 0, 0), R3(0, 0, 0));

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

	ofstream reactionfile;

	if (save == false) {
		reactionfile.open("data/anchor/reactions.txt");
	}
	for (int i = 0; i < num_steps; i++) {
		system_gpu->DoStepDynamics(timestep);
		double TIME = system_gpu->GetChTime();
		double STEP = system_gpu->GetTimerStep();
		double BROD = system_gpu->GetTimerCollisionBroad();
		double NARR = system_gpu->GetTimerCollisionNarrow();
		double LCP = system_gpu->GetTimerLcp();
		double SHUR = ((ChLcpSolverParallelDVI*) (system_gpu->GetLcpSolverSpeed()))->solver.time_shurcompliment;
		double UPDT = system_gpu->GetTimerUpdate();
		double RESID = ((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverParallelDVI*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, SHUR, UPDT, BODS, CNTC, REQ_ITS, RESID);

		system_gpu->gpu_data_manager->system_timer.PrintReport();

		int save_every = 1.0 / timestep / 60.0;     //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";

			if (save) {
				DumpAllObjectsWithGeometryChrono(system_gpu, "data/anchor/anchor.dat");
			} else {
				DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			}
			file++;
		}
		if (save == false) {

			Vector force = engine_anchor->Get_react_force();
			Vector torque = engine_anchor->Get_react_torque();
			double motor_torque = engine_anchor->Get_mot_torque();
			reactionfile << force.x << " " << force.y << " " << force.z << " " << torque.x << " " << torque.y << " " << torque.z << " " << motor_torque << " " << ANCHOR->GetPos().y
					<< " " << ANCHOR->GetPos_dt().y << endl;

		}
		RunTimeStep(system_gpu, i);
	}
	if (save == false) {
		reactionfile.close();
	}
}
