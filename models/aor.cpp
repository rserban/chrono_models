#include "unit_POSTPROCESS/CHmitsubaRender.h"
#include "physics/CHsystem.h"
#include "assets/CHsphereShape.h"
#include "physics/CHapidll.h"
#include "physics/CHsystem.h"
#include "lcp/CHlcpIterativeMINRES.h"
#include "core/CHrealtimeStep.h"
#include "lcp/CHlcpIterativeSOR.h"
#include "lcp/CHlcpIterativeJacobi.h"
#include "collision/CHcModelBullet.h"
#include "collision/CHcCollisionSystemBullet.h"
#include "physics/CHcontactContainer.h"
#include "unit_GPU/CHsystemGPU.h"
#include "unit_GPU/collision/CHcontactContainerGPU.h"
#include "unit_OPENGL/CHopenGL.h"
using namespace chrono;
using namespace postprocess;
using namespace geometry;
using namespace std;

#define PI 3.1415

ChVector<> lpos(0, 0, 0);
ChQuaternion<> quat(1, 0, 0, 0);

//all dimensions are in millimeters, milligrams
real plate_height = -10;
real plate_thickness = 1;
real plate_radius = 7;
real plate_friction = 1;

real particle_radius = .1;
real particle_mass = .5263;
real particle_density = 1.14;
real particle_friction = 1.0;
Vector particle_initial_vel = Vector(0, -2.5, 0); //initial velocity

real container_width = 7.0;
real container_thickness = .25;
real container_height = 3.0;
real wscale = 1;

real gravity = -9810;
real timestep = .00005;
int num_steps = 5000000;

int particle_grid_x = 15;
int particle_grid_z = 15;

int particles_every = 50; //add particles every n steps
int save_every = 100; //save data every n steps

int particle_configuration = 0;
//0: single sphere
//1: two spheres joined together
//2: single ellipsoid

bool create_particle_plate = false;
bool all_three_kinds = false;

string data_folder = "data";
template<class T>
void DumpObjects(T* mSys, string filename) {
	ofstream ofile(filename.c_str());

	for (int i = 0; i < mSys->Get_bodylist()->size(); i++) {
		ChBody* abody = (ChBody*) mSys->Get_bodylist()->at(i);
		if (abody->IsActive() == true) {
			ofile << particle_radius << ",";
			ofile << abody->GetPos().x << "," << abody->GetPos().y << "," << abody->GetPos().z << ",";
			ofile << abody->GetRot().e0 << "," << abody->GetRot().e1 << "," << abody->GetRot().e2 << "," << abody->GetRot().e3 << ",\n";
			// ofile << abody->GetPos_dt().x << "," << abody->GetPos_dt().y << "\t" << abody->GetPos_dt().z << "\t";
			// ofile << abody->GetWvel_loc().x << "," << abody->GetWvel_loc().y << "\t" << abody->GetWvel_loc().z << endl;
		}
	}
}
real start_height = plate_height + plate_thickness + particle_radius * 4;

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	cout << frame << endl;
	if (frame % particles_every == 0) {
		ChSharedBodyGPUPtr sphere;
		for (int i = 0; i < particle_grid_x; i++) {
			for (int j = 0; j < particle_grid_z; j++) {
				sphere = ChSharedBodyGPUPtr(new ChBodyGPU);
				Quaternion q;
				q.Q_from_NasaAngles(Vector(rand() % 1000 / 1000.0, rand() % 1000 / 1000.0, rand() % 1000 / 1000.0));
				q.Normalize();

				ChVector<> position(i * particle_radius * 2 - particle_grid_x * .5 * particle_radius * 2, start_height, j * particle_radius * 2 - particle_grid_z * .5 * particle_radius * 2);

				position.x += rand() % 1000 / 100000.0;
				position.y += rand() % 1000 / 1000.0;
				position.z += rand() % 1000 / 100000.0;

				if (particle_configuration == 0) {
					InitObject(sphere, particle_mass, position, quat, particle_friction, particle_friction, 0, true, false, -1, i);
					AddCollisionGeometry(sphere, SPHERE, ChVector<>(particle_radius, particle_radius, particle_radius), Vector(0, 0, 0), quat);
				}
				if (particle_configuration == 1) {
					position.x *= 1.2;
					position.z *= 1.0;
					InitObject(sphere, particle_mass, position, q, particle_friction, particle_friction, 0, true, false, -1, i);

					AddCollisionGeometry(sphere, SPHERE, ChVector<>(particle_radius * .75, particle_radius * .75, particle_radius * .75), Vector(-.025, 0, 0), quat);
					AddCollisionGeometry(sphere, SPHERE, ChVector<>(particle_radius * .75, particle_radius * .75, particle_radius * .75), Vector(.025, 0, 0), quat);
					real rho = particle_density;

					Vector inertia = Vector(3.47e-6 * rho, 2.62e-7 * rho, 3.47e-6 * rho);
					sphere->SetInertiaXX(inertia);
				}
				if (particle_configuration == 2) {

					InitObject(sphere, particle_mass, position, q, particle_friction, particle_friction, 0, true, false, -1, i);

					AddCollisionGeometry(sphere, ELLIPSOID, ChVector<>(particle_radius * .5, particle_radius, particle_radius * .5), Vector(0, 0, 0), quat);

					real mass = particle_mass;
					Vector r = ChVector<>(particle_radius * .5, particle_radius, particle_radius);
					Vector inertia = Vector((1 / 5.0 * mass * (r.y * r.y + r.z * r.z), 1 / 5.0 * mass * (r.x * r.x + r.z * r.z), 1 / 5.0 * mass * (r.x * r.x + r.y * r.y)));
					sphere->SetInertiaXX(inertia);

				}
				sphere->SetPos_dt(particle_initial_vel);
				FinalizeObject(sphere, (ChSystemGPU *) mSys);
			}
		}
		if (all_three_kinds) {
			particle_configuration++;
			if (particle_configuration > 2) {
				particle_configuration = 0;
			}
		}
		if (start_height < -5) {
			start_height += particle_radius * 2;
		}
	}

}
int main(int argc, char* argv[]) {
	omp_set_num_threads(4);

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
	system_gpu->SetMaxiter(20);
	system_gpu->SetIterLCPmaxItersSpeed(20);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(20);
	system_gpu->SetTol(1e-3);
	system_gpu->SetTolSpeeds(1e-3);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(1e-3);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(.6);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .1);
	mcollisionengine->setBinsPerAxis(R3(40, 15, 40));
	mcollisionengine->setBodyPerBin(100, 50);

	//=========================================================================================================
	ChMitsubaRender output(system_gpu);
	output.SetIntegrator("photonmapper");
	output.SetIntegratorOption("integer", "maxDepth", "32");
	output.SetFilm("ldrfilm");
	output.SetFilmOption("integer", "height", "1200");
	output.SetFilmOption("integer", "width", "1920");
	output.camera_target = ChVector<>(0, -5, 0);
	output.camera_origin = ChVector<>(0, -5, -40);
	output.camera_up = ChVector<>(0, 1, 0);
	output.SetDataFolder(data_folder);
	output.ExportScript("test.xml");
	//=========================================================================================================
	cout << "Mass, Radius, Friction_Sphere, Friction_Plate" << endl;
	if (argc == 9) {
		particle_mass = atof(argv[1]);
		particle_radius = atof(argv[2]);
		particle_friction = atof(argv[3]);
		plate_friction = atof(argv[4]);
		data_folder = argv[5];
		create_particle_plate = atoi(argv[6]);
		all_three_kinds = atoi(argv[7]);
		particle_configuration = atoi(argv[8]);
	}

	cout << particle_mass << " " << particle_radius << " " << particle_friction << " " << plate_friction << " " << data_folder << " " << create_particle_plate << " " << all_three_kinds << " "
			<< particle_configuration << endl;
	//=========================================================================================================
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);

	//=========================================================================================================
	real particle_plate_dim = 8;
	real plate_particle_radius = .02;
	int num_particle = particle_plate_dim / plate_particle_radius;

	if (create_particle_plate) {
		ChSharedBodyGPUPtr sphere;
		sphere = ChSharedBodyGPUPtr(new ChBodyGPU);
		InitObject(sphere, particle_mass, Vector(0, 0, 0), quat, plate_friction, plate_friction, 0, true, true, -1, -1);
		for (int i = 0; i < num_particle; i++) {
			for (int j = 0; j < num_particle; j++) {

				ChVector<> position(i * plate_particle_radius * 1.2 - num_particle * plate_particle_radius * .6, plate_height + plate_thickness + plate_particle_radius,
						j * plate_particle_radius * 1.2 - num_particle * plate_particle_radius * .6);

				position.x += rand() % 10000 / 10000.0 * plate_particle_radius * .25 - plate_particle_radius * .25 * .5;
				position.z += rand() % 10000 / 10000.0 * plate_particle_radius * .25 - plate_particle_radius * .25 * .5;

				AddCollisionGeometry(sphere, SPHERE, ChVector<>(plate_particle_radius, plate_particle_radius, plate_particle_radius), position, quat);

			}
		}
		FinalizeObject(sphere, (ChSystemGPU *) system_gpu);
	}
	ChSharedBodyGPUPtr PLATE = ChSharedBodyGPUPtr(new ChBodyGPU);
	InitObject(PLATE, 1, ChVector<>(0, plate_height, 0), quat, plate_friction, plate_friction, 0, true, true, -1000, -20000);
	AddCollisionGeometry(PLATE, BOX, ChVector<>(plate_radius, plate_thickness, plate_radius), lpos, quat);
	FinalizeObject(PLATE, (ChSystemGPU *) system_gpu);

	ChSharedBodyGPUPtr L = ChSharedBodyGPUPtr(new ChBodyGPU);
	ChSharedBodyGPUPtr R = ChSharedBodyGPUPtr(new ChBodyGPU);
	ChSharedBodyGPUPtr F = ChSharedBodyGPUPtr(new ChBodyGPU);
	ChSharedBodyGPUPtr B = ChSharedBodyGPUPtr(new ChBodyGPU);

	InitObject(L, 100000, Vector(-container_width + container_thickness, plate_height, 0), quat, plate_friction, plate_friction, 0, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_width - container_thickness, plate_height, 0), quat, plate_friction, plate_friction, 0, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, plate_height, -container_width + container_thickness), quat, plate_friction, plate_friction, 0, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, plate_height, container_width - container_thickness), quat, plate_friction, plate_friction, 0, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_height, container_width), lpos, quat);
	AddCollisionGeometry(F, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);
	AddCollisionGeometry(B, BOX, Vector(container_width * wscale, container_height, container_thickness), lpos, quat);

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);

	//=========================================================================================================
	//////Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_pos = Vector(0, -5, -40);
	openGLView.render_camera->look_at = Vector(0, -5, 0);
	openGLView.SetCustomCallback(RunTimeStep);
	openGLView.StartSpinning(window_manager);
	window_manager->CallGlutMainLoop();
	//=========================================================================================================

	int file = 0;
	for (int i = 0; i < num_steps; i++) {
		cout << "step " << i << endl;
		system_gpu->DoStepDynamics(timestep);
		RunTimeStep(system_gpu, i);
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";
			DumpObjects(system_gpu, ss.str());
			//output.ExportData(ss.str());
			file++;
		}

	}
	return 0;
}

