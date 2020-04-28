#include "nr3.h"
#include "svd.h"
#include "utilities.h"
#include <iostream>

using namespace std;

// Parameters
const Doub T1 = 1000;
const Doub T2 = 500;
const Doub epsilon1 = 0.8;
const Doub epsilon2 = 0.6;
const Doub sigma = 1.712*pow(10, -9);
const Doub d = 1;
const Doub w = 1;

// Integral limits
const Doub a = -w/2;
const Doub b = w/2;

// Constant terms
const Doub C1 = epsilon1 * sigma * pow(T1, 4);
const Doub C2 = epsilon2 * sigma * pow(T2, 4);


// F(x,y,d)
Doub F(Doub x, Doub y, Doub d) {
	const Doub denom = pow((pow(d, 2) + pow((x - y), 2)), 3/2);

	return (1/denom);
}

// Set up linear system of equations in matrixform
void MatrixForm(MatDoub& Amat, VecDoub& bvec, int N) {

	// Stepsize
	const Doub h = (b - a)/N;

	// Define A-matrix and b-vector
	MatDoub mAmat((2*(N+1)), (2*(N+1)));
	VecDoub mbvec((2*(N+1)));

	// Set up b-vector
	for (int i = 0; i < 2*(N+1); i++) {
		if (i < N + 1)
			mbvec[i] = C1;
		else
			mbvec[i] = C2;
	}

	/* Set up A-matrix by building the 4 "submatrices", consisting of
	   two identity matrices, and two matrices with F-terms. */

	// 
	for (int i = 0; i < N+1; i++) {		// initialize correct places to 1 and 0
		for (int j = 0; j < N+1; j++) {
			if (j == i) {// Put 1 in diagonal from upper left to lower right corner of matrix
				mAmat[i][j] = 1;
				mAmat[N+1+i][N+1+j] = 1;
			}
			else {
				mAmat[i][j] = 0;
				mAmat[N+1+i][N+1+j] = 0; // Put zero elsewhere in two submatrices
			}
		}
	}

	// Vectors for x- and y-values
	VecDoub x(N + 1), y(N + 1);

	// Calculate x- and y-values
	for (int i = 0; i < N + 1; i++) {
		x[i] = (a + h * i);
		y[i] = (a + h * i);
	}

	// Build submatrices consisting of F-terms 
	MatDoub subMatrix1((N + 1), (N + 1)), subMatrix2((N + 1), (N + 1));
	for (int i = 0; i < N+1; i++) {
		for (int j = 0; j < N+1; j++) {
			if (j != 0 && j != N) {
				subMatrix1[i][j] = -(1 - epsilon1)*h*F(x[i], y[j], d);
				subMatrix2[i][j] = -(1 - epsilon2)*h*F(x[j], y[i], d);

				mAmat[i][N+1+j] = subMatrix1[i][j];
				mAmat[N+1+i][j] = subMatrix2[i][j];
			}
			else {
				subMatrix1[i][j] = -(1 - epsilon1)*h*F(x[i], y[j], d) / 2;
				subMatrix2[i][j] = -(1 - epsilon2)*h*F(x[j], y[i], d) / 2;

				mAmat[i][N+1+j] = subMatrix1[i][j];
				mAmat[N+1+i][j] = subMatrix2[i][j];
			}
		}
	}

	Amat = mAmat;
	bvec = mbvec;
}

void FindQs(VecDoub& u, VecDoub& v, int N, Doub& Q1, Doub& Q2) {
	Doub mQ1 = 0, mQ2 = 0;
	const Doub h = (b - a) / N;

	// Vectors for x- and y-values
	VecDoub x(N + 1), y(N + 1);

	// Calculate x- and y-values
	for (int i = 0; i < N + 1; i++) {
		x[i] = (a + h * i);
		y[i] = (a + h * i);
	}

	VecDoub I1(N+1), I2(N+1);				// Calculate the integral
	for (int i = 0; i < N+1; i++) {
		I1[i] = 0;
		I2[i] = 0;
		for (int j = 0; j < N+1; j++) {
			if (j == 0 || j == N) {	// If endpoint
				I1[i] += h * (F(x[i], y[j], d)*v[j]) / 2;	// v
				I2[i] += h * (F(x[j], y[i], d)*u[j]) / 2;	// u
			}
			else {
				I1[i] += h * (F(x[i], y[j], d)*v[j]);		// v
				I2[i] += h * (F(x[j], y[i], d)*u[j]);		// u
			}
		}
	}

	// Calculate Q's
	for (int i = 0; i < N+1; i++) {		
		if (i != 0 && i != N) {	
			mQ1 += h * (u[i] - I1[i]);
			mQ2 += h * (v[i] - I2[i]);
		}
		else {
			mQ1 += h * (u[i] - I1[i]) / 2;
			mQ2 += h * (v[i] - I2[i]) / 2;
		}
	}

	Q1 = mQ1;
	Q2 = mQ2;
}

void Results(VecDoub u, VecDoub v, Doub Q1, Doub Q2, int N) {
	cout << "N: " << N << endl;
	cout << "u(x): " << endl <<
		"u(-1/2): " << u[0] << endl <<
		"u(-1/4): " << u[(N + 1) / 4] << endl <<
		"u(0): " << u[(N + 1) / 2] << endl <<
		"u(1/4): " << u[(N + 1) / 4] << endl <<
		"u(1/2): " << u[N] << endl << endl;

		cout << "v(y): " << endl <<
		"v(-1/2): " << v[0] << endl <<
		"v(-1/4): " << v[(N+1)/4] << endl <<
		"v(0): " << v[(N+1)/2] << endl <<
		"v(1/4): " << v[(N+1)/4] << endl <<
		"v(1/2): " << v[N] << endl;

	cout << "Q1: " << Q1 << endl << "Q2: " << Q2 << endl << endl;
}

void MakeRichTable(int size, int tablesize, VecDoub Q, MatDoub tableQ) {
	Doub Qt1, Qt2, Qt3;		// Holders for previous Q(N) values 
	Doub k;			// Holder for the order k of Q

	for (int i = 0; i < size; i++) {
		// Initialize N
		tableQ[i][0] = pow(2, 2 + i);

		tableQ[i][1] = Q[i];
		for (int j = 2; j < tablesize; j++) {
			tableQ[i][j] = 0;
		}
	}

	// Estimate the order k
	for (int i = 2; i < size; i++) {
		Qt1 = Q[i - 2];
		Qt2 = Q[i - 1];
		Qt3 = Q[i];
		tableQ[i][2] = log2((Qt1 - Qt2) / (Qt2 - Qt3));
		cout << "Order of Q: " << tableQ[i][2] << " | ";

	}
	cout << "Expected order of Q: " << endl;
	cin >> k;

	// Estimate error
	for (int i = 1; i < size; i++) {
		Qt1 = Q[i - 1];
		Qt2 = Q[i];
		tableQ[i][3] = (Qt1 - Qt2) / (pow(2, k) - 1);
	}

	// Calculate Richardson Extrapolation
	for (int i = 1; i < size; i++) {
		Qt1 = Q[i - 1];
		Qt2 = Q[i];
		tableQ[i][4] = Qt2 + (Qt2 - Qt1) / (pow(2, k) - 1);
	}

	// Calculate the order of S_R
	for (int i = 3; i < size; i++) {
		Qt1 = tableQ[i - 2][4];
		Qt2 = tableQ[i - 1][4];
		Qt3 = tableQ[i][4];
		tableQ[i][5] = log2((Qt1 - Qt2) / (Qt2 - Qt3));
		cout << "Order of S_R: " << "3";
		cout << tableQ[i][5] << " | ";
	}

	cout << "Expected order of S_R: " << endl;
	cin >> k;

	// Estimate error using Richardeeson Exprapolation
	for (int i = 2; i < size; i++) {
		Qt1 = tableQ[i - 1][5];
		Qt2 = tableQ[i][5];
		tableQ[i][6] = (Qt1 - Qt2) / (pow(2, k));
	}

	util::print(tableQ);
}




int main() {

	int size = 8;
	VecDoub Q1(size), Q2(size);

	/* Set up and solve system of linear equations for given values of N */
	// N = 4, 8, 16, 32, 64, 128, 264, 512
	for (int i = 0; i < size; i++) {
		int N = pow(2, 2 + i);

		// Define matrix A and b
		MatDoub Amat(2*(N+1), 2*(N+1));
		VecDoub bvec(2*(N+1));

		MatrixForm(Amat, bvec, N);

		// Solve linear equations using SVD
		VecDoub x(2*(N+1));
		SVD SVDmatrix(Amat);
		SVDmatrix.solve(bvec, x);

		VecDoub u(N+1), v(N+1);  	// Store u and v
		for (int i = 0; i < N+1; i++) {
			u[i] = x[i];
			v[i] = x[N+1+i];
		}


 	/* Use the results for u and v to find Q1 and Q2 */
		// Find Q1 and Q2
		FindQs(u, v, N, Q1[i], Q2[i]);

		// Print results
		Results(u, v, Q1[i], Q2[i], N);
	}


	/* Find the error of Q1 and Q2 */
	int tablesize = 7;
	MatDoub tableQ1(size, tablesize), tableQ2(size, tablesize); // N | S | k1 | Error estimate | S_R | k2

	MakeRichTable(size, tablesize, Q1, tableQ1);

	MakeRichTable(size, tablesize, Q2, tableQ2);

	system("pause");
	return 0;
}