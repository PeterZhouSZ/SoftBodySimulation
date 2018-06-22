//#include <Eigen/Core>
#include "Solver.h"
#include "SoftBody.h"
#include "Node.h"
#include "Tetrahedron.h"

#include <iostream>
#include <iomanip>

# include <cstdlib>
# include <cmath>
# include <ctime>

using namespace std;
using namespace Eigen;


Solver::Solver(vector< shared_ptr<SoftBody> > _softbodies, Integrator _time_integrator, Vector2d _damping, Vector3d _grav, bool _isSparse):
	grav(_grav), damping(_damping), softbodies(_softbodies), time_integrator(_time_integrator), isSparse(_isSparse)
{
	int nVerts = 3 * softbodies[0]->getNumNodes();
	A.resize(nVerts, nVerts);
	A_sparse.resize(nVerts, nVerts);
	K.resize(nVerts, nVerts);
	M.resize(nVerts, nVerts);
	v.resize(nVerts);
	b.resize(nVerts);
	f.resize(nVerts);
	x.resize(nVerts);
	reset();
}

void Solver::step(double h) {
	reset();
	
	auto tets = softbodies[0]->getTets();
	auto nodes = softbodies[0]->getNodes();

	Matrix3d I; 
	I.setIdentity();

	for (int i = 0; i < nodes.size(); ++i) {
		nodes[i]->clearForce();
		
		v.segment<3>(3 * i) = nodes[i]->v;
		M.block<3, 3>(3 * i, 3 * i) = I * nodes[i]->m;
		f.segment<3>(3 * i) = nodes[i]->m * grav;

		//
		if (isSparse) {
			for (int ii = 0; ii < 3; ii++) {
				A_.push_back(T(3*i+ii, 3*i+ii, nodes[i]->m + h * damping(0) * nodes[i]->m));
			}
		}
		
	}

	for (int i = 0; i < tets.size(); ++i) {
		tets[i]->computeElasticForces();

		// Assemble K matrix
		tets[i]->computeForceDifferentials();
		for (int row = 0; row < 4; row++) {
			MatrixXd Kb(12, 3);
			Kb = tets[i]->getStiffness().block<12, 3>(0, 3 * row);
			
			for (int jj = 0; jj < 4; ++jj) {
				int aa = tets[i]->nodes[jj]->i;
				int bb = tets[i]->nodes[row]->i;
				K.block<3, 3>(3 * aa, 3 * bb) += Kb.block<3, 3>(3 * jj, 0);

				if (isSparse) {
					for (int irow = 0; irow < 3; irow++) {
						for (int icol = 0; icol < 3; icol++) {
							A_.push_back(T(3 * aa + irow, 3 * bb + icol, - h * h * damping(1) * Kb(3 * jj + irow, icol)));
						}
					}
				}

			}
		}
	}

	for (int i = 0; i < nodes.size(); ++i) {
		f.segment<3>(3 * i) += nodes[i]->f;
	}

	b = M * v + h * f;


	if (isSparse) {
		A_sparse.setFromTriplets(A_.begin(), A_.end());
		ConjugateGradient< SparseMatrix<double> > cg;
		cg.setMaxIterations(25);
		cg.setTolerance(1e-3);
		cg.compute(A_sparse);
		x = cg.solveWithGuess(b, v);
	}
	else {
		A = M + h * damping(0) * M - h * h * damping(1) * K;
		x = A.ldlt().solve(b);
	}
	
	
	for (int i = 0; i < nodes.size(); ++i) {
		if (nodes[i]->fixed) {
			nodes[i]->v.setZero();
		}
		else {
			nodes[i]->v = x.segment<3>(3 * i);
		}	
	}

	for (int i = 0; i < nodes.size(); ++i) {

		if (nodes[i]->fixed) {
			nodes[i]->x = nodes[i]->x0;
		}
		else {
			nodes[i]->x += nodes[i]->v * h;

			/// TODO: USE MOSEK LATER
			if (nodes[i]->x(1) < -6.0 && nodes[i]->v(1) < -0.00001) {
				nodes[i]->v(1) = 0.0;
				nodes[i]->x(1) = -6.0;
			}
		}
		
	}
}

void Solver::reset() {
	A.setZero();
	M.setZero();
	K.setZero();
	v.setZero();
	f.setZero();
	x.setZero();
	b.setZero();
	A_.clear();
	A_sparse.setZero();
}

Solver::~Solver() {

}