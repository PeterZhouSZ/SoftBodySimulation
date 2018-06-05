#include "Solver.h"
#include "SoftBody.h"

#include <iostream>
#include <iomanip>

# include <cstdlib>

# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;
using namespace Eigen;

Solver::Solver(vector< shared_ptr<SoftBody> > _softbodies, Integrator _time_integrator)
{
	this->softbodies = _softbodies;
	this->time_integrator = _time_integrator;

	if (isReduced) {
		m = 6 * (int)softbodies.size();
		n = 6 + 1 * ((int)softbodies.size() - 1);
		M.resize(m, m);
		J.resize(m, n);
		f.resize(m);
		M.setZero();
		J.setZero();
		f.setZero();

	}
	else {
		n = 6 * (int)softbodies.size() + 6 + 5 * ((int)softbodies.size() - 1);
	}

	A.resize(n, n);
	x.resize(n);
	b.resize(n);

	A.setZero();
	x.setZero();
	b.setZero();
}

void Solver::step(double h) {

}

Solver::~Solver() {

}