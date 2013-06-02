#include "common.h"

void addHCPSheet(int grid_x, int grid_z, real height, real mass, real radius, real mu, bool active, real global_x, real global_z, Vector vel, ChSystemGPU* mSys) {
	real offset = 0;
	real x = 0, y = height, z = 0;
	ChSharedBodyGPUPtr body;
	for (int i = 0; i < grid_x; i++) {
		for (int k = 0; k < grid_z; k++) {
			body = ChSharedBodyGPUPtr(new ChBodyGPU);

			offset = (k % 2 != 0) ? radius : 0;
			x = i * 2 * radius + offset - grid_x * 2 * radius / 2.0 + global_x;
			z = k * (sqrt(3.0) * radius) - grid_z * sqrt(3.0) * radius / 2.0 + global_z;

			InitObject(body, mass, Vector(x, y, z), Quaternion(1, 0, 0, 0), mu, mu, 0, true, !active, -1, i);
			AddCollisionGeometry(body, SPHERE, ChVector<>(radius, radius, radius), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
			body->SetPos_dt(vel);
			FinalizeObject(body, (ChSystemGPU *) mSys);

		}
	}
}
void addHCPCube(int grid_x, int grid_y, int grid_z, real mass, real radius, real mu, bool active, real global_x, real global_y, real global_z, Vector V, ChSystemGPU* mSys) {
	real offset_x = 0;
	real offset_z = 0;
	real height = 0;

	for (int j = 0; j < grid_y; j++) {
		height = j * (sqrt(3.0) * radius);
		offset_x = offset_z = (j % 2 != 0) ? radius : 0;
		addHCPSheet(grid_x, grid_z, height + global_y, mass, radius, mu, active, offset_x + global_x, offset_z + global_z, V, mSys);
	}
}
void addPerturbedLayer(real3 origin, ShapeType type, real3 r, int3 num_per_dir, real3 percent_perturbation, real mass, real mu, ChSystemGPU* mSys) {
	real3 a = r * percent_perturbation;
	real3 d = a + 2 * r; //compute cell length
	real3 dp, pos;
	srand(1);
	ChSharedBodyGPUPtr body;
	int counter = 0;

	for (int i = 0; i < num_per_dir.x; i++) {
		for (int j = 0; j < num_per_dir.y; j++) {
			for (int k = 0; k < num_per_dir.z; k++) {

				body = ChSharedBodyGPUPtr(new ChBodyGPU);

				dp.x = rand() % 10000 / 10000.0 * a.x - a.x / 2.0;
				dp.y = rand() % 10000 / 10000.0 * a.y - a.y / 2.0;
				dp.z = rand() % 10000 / 10000.0 * a.z - a.z / 2.0;

				pos.x = i * d.x - num_per_dir.x * d.x * .5;
				pos.y = j * d.y - num_per_dir.y * d.y * .5;
				pos.z = k * d.z - num_per_dir.z * d.z * .5;

				pos += dp+origin+r;

				InitObject(body, mass, Vector(pos.x, pos.y, pos.z), Quaternion(1, 0, 0, 0), mu, mu, 0, true, false, -1, counter);
				AddCollisionGeometry(body, type, ChVector<>(r.x, r.y, r.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
				FinalizeObject(body, (ChSystemGPU *) mSys);

				counter++;

			}
		}
	}

}