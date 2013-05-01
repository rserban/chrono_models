#include "common.h"

////Example Input File
/*
timestep:				.01
solver:					APGD
solver_max_iterations:	100
solver_tolerance:		1e-3
compliance:				0 0 0
----------------------------
contact_recovery_speed: .6
contact_envelope:		.005
collision_BPA:			40 15 40
collision_BPB:			100 50

*/
void ReadInputFile(string input, ChSystemGPU* m_system){
ifstream ifile(input);
if(ifile.fail()){
	cout<<"CONFIG FILE FAIL"<<endl;
}

string token,value_string;
real value_real;
int value_int;
Vector value_vect;
while(ifile.fail()==false){}
	ifile>>token;
//--------------------------------------------------------------------------------
	if(token=="timestep:"){
		ifile>>value_real;
		m_system->SetStep(value_real);
	}
//--------------------------------------------------------------------------------
	if(token=="solver:"){
		ifile>>value_string
		if(value_string=="APGD")
			{((ChLcpSolverGPU *) (m_system->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);}
		else if(value_string=="JACOBI")
			{((ChLcpSolverGPU *) (m_system->GetLcpSolverSpeed()))->SetSolverType(BLOCK_JACOBI);}
		else if(value_string=="CG")
			{((ChLcpSolverGPU *) (m_system->GetLcpSolverSpeed()))->SetSolverType(CG);}
	}
//--------------------------------------------------------------------------------
	if(token=="solver_max_iterations:"){
		ifile>>value_int;
		((ChLcpSolverGPU *) (m_system->GetLcpSolverSpeed()))->SetMaxIteration(value_int);
		m_system->SetMaxiter(value_int);
		m_system->SetIterLCPmaxItersSpeed(value_int);
	}
//--------------------------------------------------------------------------------
	if(token=="solver_tolerance:"){
		ifile>>value_real;
		m_system->SetTol(value_real);
		m_system->SetTolSpeeds(value_real);
		((ChLcpSolverGPU *) (m_system->GetLcpSolverSpeed()))->SetTolerance(value_real);
	}
//--------------------------------------------------------------------------------
	if(token=="compliance:"){
		real compliance, compliance_T, alpha;
		ifile>>compliance>>compliance_T>>alpha;
		((ChLcpSolverGPU *) (m_system->GetLcpSolverSpeed()))->SetCompliance(compliance, compliance_T, alpha);
	}
//--------------------------------------------------------------------------------
	if(token=="contact_recovery_speed:"){
		ifile>>value_real;
		((ChLcpSolverGPU *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(value_real);
	}
//--------------------------------------------------------------------------------
	if(token=="contact_envelope:"){
		ifile>>value_real;
		((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(value_real);
	}
//--------------------------------------------------------------------------------
	if(token=="collision_BPA:"){
			int x, y, z;
		ifile>>x>>y>>z;
		((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->broadphase.setBinsPerAxis(R3(x, y, z));
	}
//--------------------------------------------------------------------------------
	if(token=="collision_BPB:"){
		int min_b, max_b;
		ifile>>max_b>>min_b;
		((ChCollisionSystemGPU *) (system_gpu->GetCollisionSystem()))->broadphase.setBodyPerBin(max_b, min_b);
	}
//--------------------------------------------------------------------------------
}


}
