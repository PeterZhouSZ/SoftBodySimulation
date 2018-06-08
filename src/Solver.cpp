
#include "Solver.h"
#include "SoftBody.h"
#include "Node.h"
#include "Tetrahedron.h"

#include <iostream>
#include <iomanip>

# include <cstdlib>

# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;
using namespace Eigen;

Solver::Solver(vector< shared_ptr<SoftBody> > _softbodies, Integrator _time_integrator):
	grav(0.0, -9.8, 0.0)
{
	this->softbodies = _softbodies;
	this->time_integrator = _time_integrator;
	
	int nVerts = 3 * softbodies[0]->getNumNodes();
	A.resize(nVerts, nVerts);
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

	}

	for (int i = 0; i < tets.size(); ++i) {
		tets[i]->computeElasticForces();
	}

	for (int i = 0; i < nodes.size(); ++i) {
		f.segment<3>(3 * i) += nodes[i]->f;
	}

	// TODO: SPARSE; ADD K MATRIX 
	b = M * v + h * f;
	A = M ;
	
	x = A.ldlt().solve(b);
	
	for (int i = 0; i < nodes.size(); ++i) {
		nodes[i]->v = x.segment<3>(3 * i);
	}

	for (int i = 0; i < nodes.size(); ++i) {
		nodes[i]->x += nodes[i]->v * h;

		/// TODO: USE MOSEK LATER
		if (nodes[i]->x(1) < -10.0 && nodes[i]->v(1) < -0.00001) {
			nodes[i]->v(1) = 0.0;
			nodes[i]->x(1) = -10.0;
		}
	}
}

void Solver::reset() {
	A.setZero();
	M.setZero();
	v.setZero();
	f.setZero();
	x.setZero();
	b.setZero();
}

Solver::~Solver() {

}