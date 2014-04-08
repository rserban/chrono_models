#include "../../common/common.h"

int main() {

	quaternion q1(1, 2, 3, 4);
	q1 = normalize(q1);

// real3 v(1,2,3);

// real3 Pl_1 = MatTMult(AMat(q1), v);

// cout<<Pl_1<<endl;

// real3 Pl_2 = quatRotate(v, ~q1);
// cout<<Pl_2<<endl;

// M33 Ps = XMatrix(v);
// M33 Jtemp = MatMult(AMat(q1), Ps);
// cout<<Jtemp.U;
// cout<<Jtemp.V;
// cout<<Jtemp.W;
	quaternion E1(1, 2, 3, 4);
	E1 = normalize(E1);
	real3 sbar(1, 2, 3);

//real3 vnew1 = quatRotate(sbar,~E1);
//real3 vnew2 = quatRotateMatT(sbar,E1);
//
//
//
////real3 t = 2 * cross(real3(E1.x,E1.y,E1.z), sbar);
////real3 vnew2 = sbar + E1.w * t + cross(real3(E1.x,E1.y,E1.z), t);
//
//cout<<vnew1<<endl;
////
//cout<<vnew2<<endl;;

	real3 U = normalize(real3(1, 7, 3));
	real3 V, W;

	W = cross(U, R3(0, 1, 0));
	real mzlen = length(W);
	if (mzlen < .0001) {     // was near singularity? change singularity reference custom_vector!
		real3 mVsingular = R3(1, 0, 0);
		W = cross(U, mVsingular);
		mzlen = length(W);
	}
	W /= mzlen;
	V = cross(W, U);
//
	M33 contact_plane;
	contact_plane.U = U;
	contact_plane.V = V;
	contact_plane.W = W;
//
	sbar = quatRotate(sbar, ~E1);     //A^T*s
////
////			M33 sbartilde = XMatrix(sbar);
////			M33 Jtemp = MatMult(AMat(E1), sbartilde);
////			M33 Jr = MatTMult(contact_plane, Jtemp);
	real3 col1 = cross(quatRotate(U, ~E1), sbar);
	real3 col2 = cross(quatRotate(V, ~E1), sbar);
	real3 col3 = cross(quatRotate(W, ~E1), sbar);
	cout << col1;
	cout << col2;
	cout << col3;

	real3 A = normalize(quatRotate(U, ~E1));
	real3 B, C;
	C = cross(A, R3(0, 1, 0));
	 mzlen = length(C);
	if (mzlen < .0001) {     // was near singularity? change singularity reference custom_vector!
		real3 mVsingular = R3(1, 0, 0);
		C = cross(A, mVsingular);
		mzlen = length(C);
	}
	C /= mzlen;
	B = cross(C, A);

	 col1 = cross(A, sbar);
	 col2 = cross(B, sbar);
	 col3 = cross(C, sbar);

	 cout<<endl;

	 cout << col1;
	 	cout << col2;
	 	cout << col3;
//
//
////
//	M33 Jr = MatTMult(contact_plane, AMat(E1));
////
//	real3 T1 = R3(Jr.U.x, Jr.V.x, Jr.W.x);
//	real3 T2 = R3(Jr.U.y, Jr.V.y, Jr.W.y);
//	real3 T3 = R3(Jr.U.z, Jr.V.z, Jr.W.z);
////
//	real3 col1 = quatRotate(U, ~E1);
//	real3 col2 = quatRotate(V, ~E1);
//	real3 col3 = quatRotate(W, ~E1);
//	cout << col1;
//	cout << col2;
//	cout << col3;
//
//	cout << endl;
//	cout << T1;
//	cout << T2;
//	cout << T3;
	/*
	 [-1.5185, -2.30493, 1.89737]
	 [-0.347087, -0.52684, -2.9093]
	 [1.34496, 2.04151, 1.3914]
	 */
//  [0.650787, 0.178357, 0.737865]
//  [0.0780945, -0.982338, 0.168655]
//  [0.754913, -0.0521351, -0.653537]
	return 0;

}
