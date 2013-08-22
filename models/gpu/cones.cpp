#include "../../common/common.h"
#include "../../common/generation.h"
#include "../../common/parser.h"
#include "../../common/input_output.h"
real gravity = -9.80665;
real timestep = .001;
real seconds_to_simulate = 30;

int max_iter = 20;

int num_steps = seconds_to_simulate / timestep;

real3 container_size = R3(6, 6, 6);
real container_thickness = .4;
real container_height = 0;
real container_friction = 1;

real particle_radius = .1;

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	if (mSys->GetNbodies() < 10000) {
		ChSharedBodyPtr sphere;
		real3 rad = R3(particle_radius*2, particle_radius*3, particle_radius*2);
		real3 size = container_size;
		size.y = container_size.y / 3.0;

		int3 num_per_dir = I3(1, 10, 10);
		real mu = 1;
		real mass = 1;
		real3 vel =  R3(-5, 0, 0);
		if (frame % 40 == 0) {

			ChSharedBodyPtr body;
			int counter = 0;

			for (int i = 0; i < num_per_dir.x; i++) {
				for (int j = 0; j < num_per_dir.y; j++) {
					for (int k = 0; k < num_per_dir.z; k++) {
						real3 r = rad;

						real3 a = r * R3(0, 0, 0);
						real3 d = a + 2 * r;     //compute cell length
						real3 dp, pos;

						body = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

						dp.x = rand() % 10000 / 10000.0 * a.x - a.x / 2.0;
						dp.y = rand() % 10000 / 10000.0 * a.y - a.y / 2.0;
						dp.z = rand() % 10000 / 10000.0 * a.z - a.z / 2.0;

						pos.x = i * d.x - num_per_dir.x * d.x * .5;
						pos.y = j * d.y - num_per_dir.y * d.y * .5;
						pos.z = k * d.z - num_per_dir.z * d.z * .5;

						pos += dp + R3(5, 0, 0) + r;

						InitObject(body, mass, Vector(pos.x, pos.y, pos.z), Quaternion(1, 0, 0, 0), mu, mu, 0, true, false, -1, counter);

						AddCollisionGeometry(body, CONE, ChVector<>(r.x, r.y, r.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
						AddCollisionGeometry(body, SPHERE, ChVector<>(r.x, r.y, r.z), Vector(0, -r.x, 0), Quaternion(1, 0, 0, 0));

						FinalizeObject(body, (ChSystemGPU *) mSys);
						body->GetMaterialSurface()->SetCohesion(0);
						body->SetPos_dt(Vector(vel.x, vel.y, vel.z));
						counter++;

					}
				}
			}

			//addPerturbedLayer(R3(-2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.x, friction.x, cohesion.x, R3(0, 5, 0), (ChSystemGPU*) mSys);
			//addPerturbedLayer(R3(5, 0, 0), CONE, rad, num_per_dir, R3(0, 0, 0), 1, 1, 0, R3(-5, 0, 0), (ChSystemGPU*) mSys, 0);
			//addPerturbedLayer(R3(2, 0, 0), SPHERE, rad, num_per_dir, R3(1, 0, 1), mass.z, friction.z, cohesion.z, R3(0, 5, 0), (ChSystemGPU*) mSys);
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
	system_gpu->SetMaxiter(max_iter);
	system_gpu->SetIterLCPmaxItersSpeed(max_iter);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetMaxIteration(max_iter);
	system_gpu->SetTol(0);
	system_gpu->SetTolSpeeds(0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(1);
	((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);
	((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .05);
	mcollisionengine->setBinsPerAxis(R3(30, 30, 30));
	mcollisionengine->setBodyPerBin(100, 50);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
//=========================================================================================================
//cout << num_per_dir.x << " " << num_per_dir.y << " " << num_per_dir.z << " " << num_per_dir.x * num_per_dir.y * num_per_dir.z << endl;
//addPerturbedLayer(R3(0, -5 +container_thickness-particle_radius.y, 0), ELLIPSOID, particle_radius, num_per_dir, R3(.01, .01, .01), 10, 1, system_gpu);
//addHCPCube(num_per_dir.x, num_per_dir.y, num_per_dir.z, 1, particle_radius.x, 1, true, 0,  -6 +container_thickness+particle_radius.y, 0, 0, system_gpu);
//=========================================================================================================

	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));
	ChSharedBodyPtr Top = ChSharedBodyPtr(new ChBody(new ChCollisionModelGPU));

	InitObject(
			L,
			100000,
			Vector(-container_size.x + container_thickness, container_height - container_thickness, 0),
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
			Vector(container_size.x - container_thickness, container_height - container_thickness, 0),
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
			Vector(0, container_height - container_thickness, -container_size.z + container_thickness),
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
			Vector(0, container_height - container_thickness, container_size.z - container_thickness),
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
			Vector(0, container_height - container_size.y, 0),
			Quaternion(1, 0, 0, 0),
			container_friction,
			container_friction,
			0,
			true,
			true,
			-20,
			-20);
	InitObject(
			Top,
			100000,
			Vector(0, container_height + container_size.y, 0),
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
	AddCollisionGeometry(Top, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	L->GetMaterialSurface()->SetCohesion(0);
	R->GetMaterialSurface()->SetCohesion(0);
	F->GetMaterialSurface()->SetCohesion(0);
	B->GetMaterialSurface()->SetCohesion(0);
	Bottom->GetMaterialSurface()->SetCohesion(0);
	Top->GetMaterialSurface()->SetCohesion(0);

	FinalizeObject(L, (ChSystemGPU *) system_gpu);
	FinalizeObject(R, (ChSystemGPU *) system_gpu);
	FinalizeObject(F, (ChSystemGPU *) system_gpu);
	FinalizeObject(B, (ChSystemGPU *) system_gpu);
	FinalizeObject(Bottom, (ChSystemGPU *) system_gpu);
	FinalizeObject(Top, (ChSystemGPU *) system_gpu);

//=========================================================================================================
//Rendering specific stuff:
	ChOpenGLManager * window_manager = new ChOpenGLManager();
	ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
	openGLView.render_camera->camera_pos = Vector(0, -5, -10);
	openGLView.render_camera->look_at = Vector(0, -5, 0);
	openGLView.render_camera->mScale = .5;
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
		double RESID = ((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverGPU*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

		int save_every = 1.0 / timestep / 60.0;     //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << "data/foam/" << "/" << file << ".txt";
			//DumpAllObjects(system_gpu, ss.str(), ",", true);
			DumpAllObjectsWithGeometry(system_gpu, ss.str(), ",");
			//output.ExportData(ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}

	//DumpObjects(system_gpu, "diagonal_impact_settled.txt", "\t");

}
