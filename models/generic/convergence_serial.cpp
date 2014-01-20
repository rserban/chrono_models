#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
ChVector<> lpos(0, 0, 0);
ChQuaternion<> quat(1, 0, 0, 0);

//all dimensions are in millimeters, milligrams
real plate_height = 0;
real plate_thickness = 1;
real plate_radius = 7;
real plate_friction = 1;

real particle_radius = .14;
real particle_mass = 1;
real particle_friction = 1.0;

real container_width = 3.0;
real container_thickness = .25;
real container_height = 6.0;
real wscale = 1;

int max_iter = 1000000;

real gravity = -9.810;
real timestep = .001;
real recovery_speed = 50;
real tolerance = .2;

real seconds_to_simulate = timestep;
int num_steps = seconds_to_simulate / timestep;

int save_every = 1.0 / timestep / 600.0;     //save data every n steps
ChSharedBodyPtr BLOCK;

string data_folder = "data/convergence";

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	mSys->Get_bodylist()->at(5)->SetBodyFixed(false);
	mSys->Get_bodylist()->at(5)->SetRot(ChQuaternion<>(1, 0, 0, 0));
	mSys->Get_bodylist()->at(5)->SetWvel_loc(ChVector<>(0, 0, 0));
	mSys->Get_bodylist()->at(5)->SetPos(ChVector<>(0, mSys->Get_bodylist()->at(5)->GetPos().y, 0));
	mSys->Get_bodylist()->at(5)->SetPos_dt(ChVector<>(0, mSys->Get_bodylist()->at(5)->GetPos_dt().y, 0));
	//if (frame % particles_every == 0) {
	//addHCPSheet(10, 10, 0, particle_mass, particle_radius, particle_friction, true, 0, 0, particle_initial_vel, (ChSystemParallel*) mSys);
	//}

//	/std::vector<double> violation = ((ChLcpIterativeSolver*) mSys->GetLcpSolverSpeed())->GetViolationHistory();
//	for (int i = 0; i < violation.size(); i++) {
//		cout << violation[i] << endl;
//	}
	//double residual =violation[violation.size()-1];
	//cout<<residual<<endl;
	//real residual = ((ChLcpSolverParallel *) (mSys->GetLcpSolverSpeed()))->GetResidual() * .5;
	//mSys->SetTol(residual);
	//mSys->SetTolSpeeds(residual);
	//((ChLcpSolverParallel *) (mSys->GetLcpSolverSpeed()))->SetTolerance(residual);

//	if(frame >300){
//
//		BLOCK->SetPos(ChVector<>(0, 5, 0));
//
//	}

}
int main(int argc, char* argv[]) {
	omp_set_num_threads(1);
	//max_iter = atof(argv[2]);
	//tolerance = atof(argv[3]);
	//real recovery_speed = atof(argv[4]);
	//particle_radius = atof(argv[5]);
	bool visual = false;     //atoi(argv[4]);
	string inputfile;

	string solver = argv[1];
	real block_mass = atof(argv[2]);
	num_steps = atoi(argv[3]);
	bool save = atoi(argv[4]);
	real perturb = atof(argv[5]);
	if (save) {
		data_folder = argv[6];
	} else {
		data_folder = argv[6];
		inputfile = argv[7];
	}

	cout << solver << " " << max_iter << " " << tolerance << " " << particle_radius << " " << data_folder << endl;

	//=========================================================================================================
	ChSystem * system = new ChSystem;
	system->SetParallelThreadNumber(1);
	system->SetIntegrationType(ChSystem::INT_ANITESCU);

	if (solver == "APGD") {
		system->SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
		//((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetVerbose(false);
	} else if (solver == "JACOBI") {
		system->SetLcpSolverType(ChSystem::LCP_ITERATIVE_JACOBI);
		//((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetVerbose(false);
		//((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetOmega(.5);
		//((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetSharpnessLambda(.5);
	} else if (solver == "SOR") {
		system->SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR);
		//((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetVerbose(false);
		//NORELAX
		//((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetOmega(1.2);
		//((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetSharpnessLambda(1.5);
	}

	((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetVerbose(true);
	//=========================================================================================================
	//ReadInputFile("convergence_config.txt",system_gpu);
	system->SetMaxiter(max_iter);
	system->SetIterLCPmaxItersSpeed(max_iter);
	system->SetTol(tolerance);
	system->SetTolSpeeds(tolerance);
	system->SetMaxPenetrationRecoverySpeed(recovery_speed);

	//=========================================================================================================
	system->Set_G_acc(ChVector<>(0, gravity, 0));
	system->SetStep(timestep);
	((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->SetRecordViolation(true);
	//=========================================================================================================
	if (save == false) {
		cout << "read data" << endl;
		ReadAllObjectsWithGeometryChrono(system, inputfile, false);
		system->Get_bodylist()->at(5)->SetMass(block_mass);

		real rx = container_width - container_thickness * 2;
		real ry = container_thickness;
		real rz = container_width - container_thickness * 2;

		system->Get_bodylist()->at(5)->SetInertiaXX(ChVector<>(1 / 12.0 * block_mass * (ry * ry + rz * rz), 1 / 12.0 * block_mass * (rx * rx + rz * rz), 1 / 12.0 * block_mass * (rx * rx + ry * ry)));
//
		for(int i=0; i<system->Get_bodylist()->size(); i++){
			system->Get_bodylist()->at(i)->GetMaterialSurface()->SetFriction(.1);
		}


	} else {
		ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody());
		ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody());
		ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody());
		ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody());
		ChSharedBodyPtr BTM = ChSharedBodyPtr(new ChBody());
		BLOCK = ChSharedBodyPtr(new ChBody());

		ChSharedPtr<ChMaterialSurface> material;
		material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
		material->SetFriction(plate_friction);
		material->SetCompliance(0);

		InitObject(L, 100000, Vector(-container_width + container_thickness, plate_height, 0), quat, material, true, true, 2, 4);
		InitObject(R, 100000, Vector(container_width - container_thickness, plate_height, 0), quat, material, true, true, 2, 4);
		InitObject(F, 100000, Vector(0, plate_height, -container_width + container_thickness), quat, material, true, true, 2, 4);
		InitObject(B, 100000, Vector(0, plate_height, container_width - container_thickness), quat, material, true, true, 2, 4);
		InitObject(BTM, 100000, Vector(0, plate_height - container_height, 0), quat, material, true, true, 2, 4);

		ChSharedPtr<ChMaterialSurface> material_block;
		material_block = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
		material_block->SetFriction(0);
		material_block->SetCompliance(0);

		InitObject(BLOCK, block_mass, Vector(0, 5 - .5, 0), quat, material_block, true, false, 2, 4);

		AddCollisionGeometry(L, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
		AddCollisionGeometry(R, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
		AddCollisionGeometry(F, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
		AddCollisionGeometry(B, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
		AddCollisionGeometry(BTM, BOX, Vector(container_width * wscale, container_thickness, container_width), lpos, quat);
		AddCollisionGeometry(BLOCK, BOX, Vector(container_width - container_thickness * 2, container_thickness, container_width - container_thickness * 2), lpos, quat);

		FinalizeObject(L, system);
		FinalizeObject(R, system);
		FinalizeObject(F, system);
		FinalizeObject(B, system);
		FinalizeObject(BTM, system);
		FinalizeObject(BLOCK, system);
		BLOCK->SetPos_dt(ChVector<>(0, -4, 0));

		ParticleGenerator<ChSystem> layer_gen(system, false);
		layer_gen.SetMass(1);
		layer_gen.SetRadius(R3(particle_radius));
		//layer_gen.SetNormalDistribution(particle_radius, .01);
		layer_gen.material->SetFriction(.1);
		layer_gen.material->SetRollingFriction(0);
		layer_gen.material->SetSpinningFriction(0);
		layer_gen.material->SetCohesion(0);
		layer_gen.material->SetCompliance(0);
		int3 num_per_dir = I3(1.7 / particle_radius, 4.0 / particle_radius, 1.7 / particle_radius);

		layer_gen.AddMixtureType(MIX_SPHERE);
		//layer_gen.AddMixtureType(MIX_ELLIPSOID);
		//layer_gen.AddMixtureType(MIX_CONE);
		//layer_gen.AddMixtureType(MIX_CUBE);
		//layer_gen.AddMixtureType(MIX_CYLINDER);

		//layer_gen.SetNormalDistribution(rad.x, rad.x/4.0);
		//layer_gen->UseNormalCohesion(particle_cohesion, 1);

		layer_gen.addPerturbedVolumeMixture(R3(0, -.8, 0), num_per_dir, R3(perturb, 0, perturb), R3(0, -4, 0));
	}
	//layer_gen.addPerturbedVolume(R3(0, -.2, 0), SPHERE, num_per_dir, R3(0.0, 0.0, 0.0), R3(0, -4, 0), false);
	//=========================================================================================================
	//////Rendering specific stuff:
	if (visual) {
		ChOpenGLManager * window_manager = new ChOpenGLManager();
		ChOpenGL openGLView(window_manager, system, 800, 600, 0, 0, "Test_Solvers");

		//openGLView.render_camera->camera_pos = Vector(0, -5, -40);
		//openGLView.render_camera->look_at = Vector(0, -5, 0);
		//openGLView.render_camera->mScale = .1;

		openGLView.SetCustomCallback(RunTimeStep);
		openGLView.StartSpinning(window_manager);
		window_manager->CallGlutMainLoop();
	}
	//=========================================================================================================
//	ofstream ofile("convergence.txt");
	ChTimer<double> timer;
	timer.start();
	int file = 0;
	stringstream s1;
	s1 << data_folder << "/residual.txt";
	CSVGen csv_output;
	csv_output.OpenFile(s1.str().c_str());
	for (int i = 0; i < num_steps; i++) {
//
//		cout << "step " << i << endl;
		system->DoStepDynamics(timestep);
		RunTimeStep(system, i);
		double TIME = system->GetChTime();
		double STEP = system->GetTimerStep();
		double BROD = system->GetTimerCollisionBroad();
		double NARR = system->GetTimerCollisionNarrow();
		double LCP = system->GetTimerLcp();
		double UPDT = system->GetTimerUpdate();
		std::vector<double> violation = ((ChLcpIterativeSolver*) system->GetLcpSolverSpeed())->GetViolationHistory();
		int REQ_ITS = violation.size();
		double RESID = 0;
		if (REQ_ITS != 0) {
			RESID = violation.at(violation.size() - 1);
		}

		int BODS = system->GetNbodies();
		int CNTC = system->GetNcontacts();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);
		//if (i % save_every == 0) {
		stringstream ss, s2;
		cout << "Frame: " << file << endl;
		ss << data_folder << "/reshist" << file << ".txt";
		s2 << data_folder << "/" << file << ".txt";
		DumpResidualHist(system, ss.str());
		//DumpAllObjectsWithGeometryPovray(system, "wakawaka.txt");

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
		file++;
		//}
		timer.stop();

	}
	csv_output.CloseFile();
	cout << "TIME: " << timer() << endl;
//	ofile.close();
	if (save == true) {
		stringstream asd;
		asd << data_folder << "/dump.txt";
		DumpAllObjectsWithGeometryChrono(system, asd.str());
	}
//	DumpAllObjectsWithGeometryChrono(system, "dump_matlab.txt");
//	chrono::ChSparseMatrix mdM;
//	chrono::ChSparseMatrix mdCq;
//	chrono::ChSparseMatrix mdE;
//	chrono::ChMatrixDynamic<double> mdf;
//	chrono::ChMatrixDynamic<double> mdb;
//	chrono::ChMatrixDynamic<double> mdfric;
//	system->GetLcpSystemDescriptor()->ConvertToMatrixForm(&mdCq, &mdM, &mdE, &mdf, &mdb, &mdfric);
//
//	cout<<mdCq.GetRows()<<endl;
//
//	chrono::ChStreamOutAsciiFile file_M("H_dump_M.dat");
//	mdM.StreamOUTsparseMatlabFormat(file_M);
//	chrono::ChStreamOutAsciiFile file_Cq("H_dump_Cq.dat");
//	mdCq.StreamOUTsparseMatlabFormat(file_Cq);
//	chrono::ChStreamOutAsciiFile file_E("H_dump_E.dat");
//	mdE.StreamOUTsparseMatlabFormat(file_E);
//	chrono::ChStreamOutAsciiFile file_f("H_dump_f.dat");
//	mdf.StreamOUTdenseMatlabFormat(file_f);
//	chrono::ChStreamOutAsciiFile file_b("H_dump_b.dat");
//	mdb.StreamOUTdenseMatlabFormat(file_b);
//	chrono::ChStreamOutAsciiFile file_fric("H_dump_fric.dat");
//	mdfric.StreamOUTdenseMatlabFormat(file_fric);

	return 0;
}

