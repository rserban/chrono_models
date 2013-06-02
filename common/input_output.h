#include "common.h"

template<class T>
void DumpObjects(T* mSys, string filename, string delim = ",", bool dump_vel_rot = false) {
	ofstream ofile(filename.c_str());

	for (int i = 0; i < mSys->Get_bodylist()->size(); i++) {
		ChBody* abody = (ChBody*) mSys->Get_bodylist()->at(i);
		if (abody->IsActive() == true) {
			ofile << abody->GetPos().x << delim << abody->GetPos().y << delim << abody->GetPos().z << delim;
			ofile << abody->GetRot().e0 << delim << abody->GetRot().e1 << delim << abody->GetRot().e2 << delim << abody->GetRot().e3;
			if (!dump_vel_rot) {
				ofile << delim << endl;
			}
			if (dump_vel_rot) {
				ofile <<delim<< abody->GetPos_dt().x << delim << abody->GetPos_dt().y << delim << abody->GetPos_dt().z << delim;
				ofile <<delim<< abody->GetWvel_loc().x << delim << abody->GetWvel_loc().y << delim << abody->GetWvel_loc().z << delim << endl;
			}
		}
	}
}
template<class T>
void DumpAllObjects(T* mSys, string filename, string delim = ",", bool dump_vel_rot = false) {
	ofstream ofile(filename.c_str());

	for (int i = 0; i < mSys->Get_bodylist()->size(); i++) {
		ChBody* abody = (ChBody*) mSys->Get_bodylist()->at(i);
		ofile << abody->GetPos().x << delim << abody->GetPos().y << delim << abody->GetPos().z << delim;
		ofile << abody->GetRot().e0 << delim << abody->GetRot().e1 << delim << abody->GetRot().e2 << delim << abody->GetRot().e3;
		if (!dump_vel_rot) {
			ofile << delim << endl;
		}
		if (dump_vel_rot) {
			ofile <<delim<< abody->GetPos_dt().x << delim << abody->GetPos_dt().y << delim << abody->GetPos_dt().z << delim;
			ofile <<delim<< abody->GetWvel_loc().x << delim << abody->GetWvel_loc().y << delim << abody->GetWvel_loc().z << delim << endl;
		}
	}
}
template<class T>
void TimingFile(T* mSys, string filename, real current_time) {
	ofstream ofile(filename.c_str(), std::ofstream::out | std::ofstream::app);

	ofile << " Residual: " << ((ChLcpSolverGPU *) (mSys->GetLcpSolverSpeed()))->GetResidual();
	ofile << " ITER: " << ((ChLcpSolverGPU *) (mSys->GetLcpSolverSpeed()))->GetTotalIterations();
	ofile << " OUTPUT STEP: Time= " << current_time << " bodies= " << mSys->GetNbodies() << " contacts= " << mSys->GetNcontacts() << " step time=" << mSys->GetTimerStep() << " lcp time="
			<< mSys->GetTimerLcp() << " CDbroad time=" << mSys->GetTimerCollisionBroad() << " CDnarrow time=" << mSys->GetTimerCollisionNarrow() << " Iterations="
			<< ((ChLcpSolverGPU*) (mSys->GetLcpSolverSpeed()))->GetTotalIterations() << "\n";
	ofile.close();
}
