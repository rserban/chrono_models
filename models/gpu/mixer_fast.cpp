#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 2;
real tolerance = 4;

#define USEGPU

#ifdef USEGPU
#define ch_body ChBody(new ChCollisionModelParallel)
#define ch_system ChSystemParallel
#else
#define ch_body ChBody()
#define ch_system ChSystem
#endif

int max_iter = 10000;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(6, 6, 6);
real container_thickness = .4;
real container_height = 0;
real container_friction = 1;

real particle_radius = .04;
real particle_mass = .05;
real particle_density = .5;
real particle_friction = 0;
real particle_cohesion = .1;
real ang = 2 * CH_C_PI;

ParticleGenerator<ch_system> *layer_gen;

ChSharedBodyPtr spinner;

string data_folder = "data/mixer";
int max_particles = 1000;
template<class T>
void RunTimeStep(T* mSys, const int frame) {
	if (mSys->GetNbodies() < max_particles) {
		if (frame % int(100*particle_radius/.2) == 0) {
			layer_gen->addPerturbedVolumeMixture(R3(0, 0, 0), I3(100, 1, 100), R3(.1, 0, .1), R3(0, -5, 0));
		}
	}

	ang -= CH_C_PI * timestep / 2.0;
	if (ang <= 0) {
		ang = 2 * CH_C_PI;
	}
	Quaternion q1;
	q1.Q_from_AngY(ang);
	spinner->SetPos(Vector(0, container_height - container_size.y + 2, 0));
	spinner->SetPos_dt(Vector(0, 0, 0));
	//spinner->SetRot(q1);

	//spinner->SetWvel_loc(Vector(0, -CH_C_PI / 2.0, 0));
}

int main(int argc, char* argv[]) {
	//int threads = 8;
	string solver = "APGD";

	if (argc > 1) {
		solver = argv[1];
		max_particles = atoi(argv[2]);
		data_folder = argv[3];

		//threads = (atoi(argv[1]));
	}
	if (argc == 5) {
		particle_radius = atof(argv[4]);
	}
	omp_set_num_threads(1);
//=========================================================================================================
	ch_system * system_gpu = new ch_system;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);
#ifndef USEGPU
	if (solver == "APGD") {
		system_gpu->SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
	} else if (solver == "JACOBI") {
		system_gpu->SetLcpSolverType(ChSystem::LCP_ITERATIVE_JACOBI);
		//((ChLcpIterativeSolver*) system_gpu->GetLcpSolverSpeed())->SetOmega(.5);
		//((ChLcpIterativeSolver*) system_gpu->GetLcpSolverSpeed())->SetSharpnessLambda(.5);
	} else if (solver == "SOR") {
		system_gpu->SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR);
	}
#endif
	((ChLcpIterativeSolver*) system_gpu->GetLcpSolverSpeed())->SetRecordViolation(true);
//=========================================================================================================
	system_gpu->SetParallelThreadNumber(1);
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	system_gpu->SetTol(tolerance);
	system_gpu->SetTolSpeeds(tolerance);
	system_gpu->SetMaxPenetrationRecoverySpeed(30);
#ifdef USEGPU
	//((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iteration);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(max_iter/2);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(max_iter/2);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	((ChLcpSolverParallel*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(0);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(tolerance);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(1e-4, 1e-4, .2);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(5);
	((ChLcpSolverParallel *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .01);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(30, 30, 30));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(100, 50);
	((ch_system*) system_gpu)->SetAABB(R3(-6, -6, -6), R3(6, 6, 6));
#endif
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ch_body);
	ChSharedBodyPtr R = ChSharedBodyPtr(new ch_body);
	ChSharedBodyPtr F = ChSharedBodyPtr(new ch_body);
	ChSharedBodyPtr B = ChSharedBodyPtr(new ch_body);
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ch_body);
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ch_body);
	ChSharedBodyPtr Tube = ChSharedBodyPtr(new ch_body);
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(.4);
	material->SetRollingFriction(0);
	material->SetSpinningFriction(0);
	material->SetCompliance(0);
	material->SetCohesion(0);

	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, 1, 1);
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	FinalizeObject(Bottom, (ch_system *) system_gpu);

	ChSharedBodyPtr Ring;

	for (int i = 0; i < 360; i += 5) {

		real angle = i * PI / 180.0;
		real x = cos(angle) * 5;
		real z = sin(angle) * 5;
		Quaternion q;
		q.Q_from_AngAxis(-angle, Vector(0, 1, 0));

		Ring = ChSharedBodyPtr(new ch_body);
		InitObject(Ring, 100000, Vector(x, container_height - 3, z), q, material, true, true, 1, 1);

		AddCollisionGeometry(Ring, BOX, Vector(.25, 3, .25), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

		FinalizeObject(Ring, (ch_system *) system_gpu);
	}

	spinner = ChSharedBodyPtr(new ch_body);
	InitObject(spinner, 100000, Vector(0, container_height - container_size.y + 2, 0), Quaternion(1, 0, 0, 0), material, true, false, 1, 1);
	real spinner_h = .5;
	AddCollisionGeometry(spinner, CYLINDER, Vector(.5, .4, .5), Vector(0, 0, 0), chrono::Q_from_AngAxis(0, ChVector<>(0, 0, 1)));
	AddCollisionGeometry(spinner, BOX, Vector(container_thickness / 15.0, spinner_h, 1.5), Vector(0, 0, 1.75), chrono::Q_from_AngAxis(CH_C_PI / 4.0, ChVector<>(0, 0, 1)));
	AddCollisionGeometry(spinner, BOX, Vector(container_thickness / 15.0, spinner_h, 1.5), Vector(0, 0, -1.75), chrono::Q_from_AngAxis(-CH_C_PI / 4.0, ChVector<>(0, 0, 1)));

	AddCollisionGeometry(spinner, BOX, Vector(1.5, spinner_h, container_thickness / 15.0), Vector(1.75, 0, 0), chrono::Q_from_AngAxis(CH_C_PI / 4.0, ChVector<>(1, 0, 0)));
	AddCollisionGeometry(spinner, BOX, Vector(1.5, spinner_h, container_thickness / 15.0), Vector(-1.75, 0, 0), chrono::Q_from_AngAxis(-CH_C_PI / 4.0, ChVector<>(1, 0, 0)));

	FinalizeObject(spinner, (ch_system *) system_gpu);
#ifdef USEGPU
	layer_gen = new ParticleGenerator<ch_system>((ch_system *) system_gpu, true);
#else
	layer_gen = new ParticleGenerator<ch_system>((ch_system *) system_gpu, false);
#endif
	layer_gen->SetDensity(1000);

	layer_gen->SetRadius(R3(particle_radius));

	layer_gen->material->SetFriction(.5);
	layer_gen->material->SetCohesion(particle_cohesion);

#ifdef USEGPU
	layer_gen->material->SetCompliance(1e-4);
#else
	layer_gen->material->SetCompliance(0);
#endif
	layer_gen->material->SetSpinningFriction(0);
	layer_gen->material->SetRollingFriction(0);
	layer_gen->SetRadius(R3(particle_radius, particle_radius * 1.1, particle_radius));
	layer_gen->SetCylinderRadius(4.5);
	layer_gen->SetNormalDistribution(particle_radius, .01);
	layer_gen->AddMixtureType(MIX_SPHERE);
	layer_gen->AddMixtureType(MIX_ELLIPSOID);
	//layer_gen->AddMixtureType(MIX_CYLINDER);
	layer_gen->AddMixtureType(MIX_CUBE);
	//layer_gen->AddMixtureType(MIX_CONE);

//=========================================================================================================
////Rendering specific stuff:
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
	stringstream s1;
	s1 << data_folder << "/residual.txt";
	CSVGen csv_output;
	csv_output.OpenFile(s1.str().c_str());
	int file = 0;
	for (int i = 0; i < num_steps; i++) {
		system_gpu->DoStepDynamics(timestep);
		double TIME = system_gpu->GetChTime();
		double STEP = system_gpu->GetTimerStep();
		double BROD = system_gpu->GetTimerCollisionBroad();
		double NARR = system_gpu->GetTimerCollisionNarrow();
		double LCP = system_gpu->GetTimerLcp();
		double UPDT = system_gpu->GetTimerUpdate();
		std::vector<double> violation = ((ChLcpIterativeSolver*) system_gpu->GetLcpSolverSpeed())->GetViolationHistory();
		int REQ_ITS = violation.size();
		double RESID = 0;
		if (REQ_ITS != 0) {
			RESID = violation.at(violation.size() - 1);
		}

		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

		int save_every = 1.0 / timestep / 60.0;     //save data every n steps
//		if (i % save_every == 0) {
//			stringstream ss;
//			cout << "Frame: " << file << endl;
//			ss << data_folder << "/" << file << ".txt";
//			//DumpAllObjects(system_gpu, ss.str(), ",", true);
//			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
//			file++;
//
//		}
		csv_output << TIME;
		csv_output << STEP;
		csv_output << BROD;
		csv_output << NARR;
		csv_output << LCP;
		csv_output << UPDT;
		csv_output << BODS;
		csv_output << CNTC;
		csv_output << REQ_ITS;
		csv_output << RESID;

		csv_output.Endline();

		RunTimeStep(system_gpu, i);
	}
	csv_output.CloseFile();

}
