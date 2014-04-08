#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .0002;
real particle_radius = .1;
real particle_friction = .5;
real seconds_to_simulate = 200;
real tolerance = 0;

double fstar = .27;
double Gamma = 3.0;
double phi = .58;

double D = particle_radius * 2;
double L = D * 100;
double P = 60000;

real3 container_size = R3(L / 2, 5, L / 2);
real container_thickness = .25;
real container_height = 4;
real container_friction = .5;
real current_time = 0;
//int3 num_per_dir = I3((container_size.x-particle_radius-container_thickness*2)/particle_radius,10,(container_size.z-particle_radius-container_thickness*2)/particle_radius);
int3 num_per_dir = I3(50, 24, 50);
int save_every = 1.0 / timestep / 60.0;     //save data every n steps

double H = P * PI / 6.0 * (pow(D, 3) / pow(L, 2)) / phi;
double frequency = fstar / sqrtf(H / fabs(gravity));
double amplitude = Gamma * fabs(gravity) * pow(PI, -0.2e1) * pow(frequency, -0.2e1) / 0.4e1;

int max_iter = 80;
string data_folder = "data/shaker2";
int num_steps = seconds_to_simulate / timestep;

ChSharedBodyPtr Bottom;
template<class T>
void RunTimeStep(T* mSys, const int frame) {
	//cout<<"Freq: = "<<frequency<<" Amplitude: = "<<amplitude<<endl;
	Bottom->SetPos(Vector(0, (sin(frame * timestep * PI * 2 * frequency)) * amplitude + container_height - container_size.y, 0));
	Bottom->SetPos_dt(Vector(0, (cos(frame * timestep * PI * 2 * frequency)) * amplitude * 2 * PI * frequency, 0));
}

int main(int argc, char* argv[]) {
	bool warm_start = false;
	if (argc == 7) {
		//omp_set_num_threads(atoi(argv[1]));

		fstar = atof(argv[2]);
		Gamma = atof(argv[3]);

		H = P * PI / 6.0 * (pow(D, 3) / pow(L, 2)) / phi;
		frequency = fstar / sqrtf(H / fabs(gravity));
		amplitude = Gamma * fabs(gravity) * pow(PI, -0.2e1) * pow(frequency, -0.2e1) / 0.4e1;
		cout << "Freq: = " << frequency << " Amplitude: = " << amplitude << endl;
		seconds_to_simulate = atof(argv[4]);
		num_steps = seconds_to_simulate / timestep;
		data_folder = argv[5];

		warm_start = atoi(argv[6]);

	} else {
		omp_set_num_threads(1);
	}

	//=========================================================================================================
	ChSystemGPU * system_gpu = new ChSystemGPU;
	ChCollisionSystemGPU *mcollisionengine = new ChCollisionSystemGPU();
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

	//=========================================================================================================
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(tolerance);
	system_gpu->SetTolSpeeds(tolerance);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(tolerance);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(ceil(2 * PI * frequency * amplitude));     //IMPORTANT: this is the max velocity of the plate!!
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(warm_start);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	mcollisionengine->setBinsPerAxis(R3(num_per_dir.x * 2, num_per_dir.y * 2, num_per_dir.z * 2));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
	//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

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

	AddCollisionGeometry(
			Bottom,
			BOX,
			Vector(container_thickness, container_size.y * 1.5, container_size.z),
			Vector(-container_size.x + container_thickness, container_height - container_thickness, 0),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(
			Bottom,
			BOX,
			Vector(container_thickness, container_size.y * 1.5, container_size.z),
			Vector(container_size.x - container_thickness, container_height - container_thickness, 0),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(
			Bottom,
			BOX,
			Vector(container_size.x, container_size.y * 1.5, container_thickness),
			Vector(0, container_height - container_thickness, -container_size.z + container_thickness),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(
			Bottom,
			BOX,
			Vector(container_size.x, container_size.y * 1.5, container_thickness),
			Vector(0, container_height - container_thickness, container_size.z - container_thickness),
			Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	//FinalizeObject(L, (ChSystemGPU *) system_gpu);
	//FinalizeObject(R, (ChSystemGPU *) system_gpu);
	//FinalizeObject(F, (ChSystemGPU *) system_gpu);
	//FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);
	//=========================================================================================================
	cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;

	ParticleGenerator layer_gen(system_gpu);
	layer_gen.SetDensity(11340);
	layer_gen.SetRadius(R3(particle_radius));
	layer_gen.SetNormalDistribution(particle_radius, 1 / 300.0);
	layer_gen.material->SetFriction(particle_friction);
	layer_gen.material->SetCohesion(0);
	layer_gen.material->SetCompliance(0);

	//layer_gen.addHCPCube(num_per_dir,1,R3(0,0,0),R3(0,0,0));
	//num_per_dir=I3(10,10,10);
	layer_gen.addVolume(R3(0, 2, 0), SPHERE, num_per_dir, R3(0, 0, 0));

	//addPerturbedLayer(R3(0, -5 + particle_radius + container_thickness, 0), ELLIPSOID, R3(particle_radius), num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
	//=========================================================================================================
	//////Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	openGLView.render_camera->camera_position = glm::vec3(0, -5, -10);
//	openGLView.render_camera->camera_look_at = glm::vec3(0, -5, 0);
//	openGLView.render_camera->camera_scale = .1;
//	openGLView.SetCustomCallback(RunTimeStep);
//	openGLView.StartSpinning(window_manager);
//	window_manager->CallGlutMainLoop();
	//=========================================================================================================
	int file = 0;

	stringstream ss_m;

	ss_m << data_folder << "/" << "timing.txt";
	string timing_file_name = ss_m.str();
	ofstream ofile(timing_file_name.c_str());
	ofile.close();

	for (int i = 0; i < num_steps; i++) {
		current_time += timestep;
		system_gpu->DoStepDynamics(timestep);
		cout << "step " << i;
		cout << " Residual: " << ((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		cout << " ITER: " << ((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();
		cout << " OUTPUT STEP: Time= " << current_time << " bodies= " << system_gpu->GetNbodies() << " contacts= " << system_gpu->GetNcontacts() << " step time=" << system_gpu->GetTimerStep()
				<< " lcp time=" << system_gpu->GetTimerLcp() << " CDbroad time=" << system_gpu->GetTimerCollisionBroad() << " CDnarrow time=" << system_gpu->GetTimerCollisionNarrow() << " Iterations="
				<< ((ChLcpSolverGPU*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations() << "\n";

		TimingFile(system_gpu, timing_file_name, current_time);

		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}
}
