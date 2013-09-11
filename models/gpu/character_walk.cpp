#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real particle_radius = .025;
real particle_friction = 1;
real seconds_to_simulate = 14;

real3 container_size = R3(5, 5, 20);
real container_thickness = .25;
real container_height = 4;
real container_friction = 1;
real current_time = 0;
int3 num_per_dir = I3((container_size.x - particle_radius - container_thickness * 2) / particle_radius, 20, (container_size.z - particle_radius - container_thickness * 2) / particle_radius);
//int3 num_per_dir = I3(35, 10, 130);
int save_every = 1.0 / timestep / 30.0;     //save data every n steps

int max_iter = 20;

int num_steps = seconds_to_simulate / timestep;

ChSharedBodyPtr Bottom;

real char_fps = 30;
int char_load_every = 1.0 / timestep / char_fps;     //load data every n steps

int char_init = 0;     //initial data frame
int char_final = 400;     //final data frame
int char_file = 1;
int char_next = 2;
int current_alpha = 0;
real char_scale = .06;
const int char_objects = 66;

string prefix = "cat/cat";
string char_init_file = "cat/cat_init.txt";
ChSharedBodyPtr Character[char_objects];

vector<real3> char_pos[char_objects];
vector<real4> char_rot[char_objects];

static inline quaternion slerp_(const quaternion &a, const quaternion &b, real alpha) {

	real dotp = dot(a, b);

	if (dotp > .9995) {
		return lerp(a, b, alpha);
	}
	//clamp
	dotp = fmin(1.f, fmax(-1.f, dotp));

	real theta_0 = acosf(dotp);

	real theta = theta_0 * alpha;

	quaternion c = b - a * dotp;

	c = normalize(c);

	quaternion ans = a * cos(theta) + c * sin(theta);
	;

	return ans;

}

void LoadChar() {
	string temp;
	for (int i = char_init; i < char_final; i++) {
		stringstream ss;
		ss << prefix << i << ".txt";
		//cout << ss.str() << endl;
		ifstream finit(ss.str().c_str());
		for (int j = 0; j < char_objects; j++) {
			getline(finit, temp);
			if (finit.fail() == false) {
				stringstream ss(temp);
				real4 quat;
				real3 pos;
				ss >> pos.x >> pos.y >> pos.z >> quat.w >> quat.x >> quat.y >> quat.z;
				pos = pos * char_scale + R3(0, 0, 35);
				char_pos[j].push_back(pos);
				char_rot[j].push_back(quat);
			}
		}
	}
}

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	real3 offset = R3(0, -.3, 0);

	if (frame % char_load_every == 0) {
		char_file++;
		char_next++;
		current_alpha = 0;
	}

	current_alpha++;
	if (char_next < char_final) {
		//cout << current_alpha / float(char_load_every) << endl;
		for (int i = 0; i < char_objects; i++) {
			real3 pos_a = char_pos[i].at(char_file) + offset;
			real3 pos_b = char_pos[i].at(char_next) + offset;
			real3 pos = lerp(pos_a, pos_b, current_alpha / float(char_load_every));
			real3 vel = (pos_b - pos_a) / timestep;
			real4 rot_a = char_rot[i].at(char_file);
			real4 rot_b = char_rot[i].at(char_next);
			real4 rot = slerp_(rot_a, rot_b, current_alpha / float(char_load_every));
			real4 rot_dt = (rot_b - rot_a) / timestep;

			ChQuaternion<> q(rot_dt.w, rot_dt.x, rot_dt.y, rot_dt.z);
			ChVector<> wvelloc;
			q.Qdt_to_Wrel(wvelloc, q);

			Character[i]->SetPos(Vector(pos.x, pos.y, pos.z));
			Character[i]->SetPos_dt(Vector(vel.x, vel.y, vel.z));

			Character[i]->SetRot(Quaternion(rot.w, rot.x, rot.y, rot.z));
			Character[i]->SetWvel_loc(wvelloc);
		}
	}
}

int main(int argc, char* argv[]) {

	if (argc == 2) {
		omp_set_num_threads(atoi(argv[1]));

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
	system_gpu->SetTol(1e-3);
	system_gpu->SetTolSpeeds(1e-3);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(1e-3);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(10);     //IMPORTANT: this is the max velocity of the plate!!
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	mcollisionengine->setBinsPerAxis(R3(170, 20, 194 * 4));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
	//=========================================================================================================

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);
	material->SetCompliance(0);
	material->SetCohesion(0);

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 11), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 11), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness + 11), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness + 11), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, +11), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);
	//=========================================================================================================
	ifstream ifile(char_init_file.c_str());
	string temp;
	for (int i = 0; i < char_objects; i++) {
		getline(ifile, temp);
		if (ifile.fail() == false) {
			stringstream ss(temp);
			real3 size, local_pos = R3(0, 0, 0);
			ss >> size.x >> size.y >> size.z;
			//ss >> local_pos.x >> local_pos.y >> local_pos.z;
			size = size * char_scale;
			local_pos = local_pos * char_scale;
			Character[i] = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
			InitObject(Character[i], 10, Vector(0, 0, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
			AddCollisionGeometry(Character[i], BOX, Vector(size.x, size.y, size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
			FinalizeObject(Character[i], (ChSystemGPU *) system_gpu);
		}
	}
	ifstream finit("cat/cat1.txt");
	for (int i = 0; i < char_objects; i++) {
		getline(finit, temp);
		if (finit.fail() == false) {
			stringstream ss(temp);
			real junk;
			real4 quat;
			real3 pos;
			ss >> pos.x >> pos.y >> pos.z >> quat.w >> quat.x >> quat.y >> quat.z;
			pos = pos * char_scale;

			Character[i]->SetPos(Vector(pos.x, pos.y, pos.z));
			Character[i]->SetRot(Quaternion(quat.w, quat.x, quat.y, quat.z));
		}
	}
	LoadChar();
	//170
	num_per_dir = I3(170, 20, 194 * 4);
	cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
	//num_per_dir = I3(1, 1, 194*4);
	ParticleGenerator layer_gen(system_gpu);
	layer_gen.SetDensity(1000);
	layer_gen.SetRadius(R3(particle_radius));
	//layer_gen.SetNormalDistribution(particle_radius , 1/300.0);
	layer_gen.material->SetFriction(particle_friction);
	layer_gen.material->SetCohesion(2);
	layer_gen.material->SetCompliance(0);

	//layer_gen.addHCPCube(num_per_dir,1,R3(0,0,0),R3(0,0,0));
	//layer_gen.addVolume(R3(0, 0, 11), SPHERE, num_per_dir, R3(0, 0, 0));

	//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius, particle_friction, true, 0, 0, 0, 0, system_gpu);
	//addPerturbedLayer(R3(0,1+particle_radius + container_thickness, 0), ELLIPSOID, R3(particle_radius), num_per_dir, R3(.01, .01, .01), 1, 1,.5,R3(0,0,0), system_gpu);
	//=========================================================================================================
	//////Rendering specific stuff:
//	ChOpenGLManager * window_manager = new ChOpenGLManager();
//	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
//	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
//	openGLView.render_camera->look_at = Vector(0, -5, 0);
//	openGLView.render_camera->mScale = 1;
//	openGLView.SetCustomCallback(RunTimeStep);
//	openGLView.StartSpinning(window_manager);
//	window_manager->CallGlutMainLoop();
	//=========================================================================================================
	int file = 0;

	stringstream ss_m;
	string data_folder = "data/snow_walk";

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

//		TimingFile(system_gpu, timing_file_name, current_time);

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

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
