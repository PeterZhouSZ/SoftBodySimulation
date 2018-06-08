#pragma once

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class SoftBody;
enum Integrator { RKF45, SYMPLECTIC, IMPLICIT };

class Solver
{
public:
	Solver(std::vector< std::shared_ptr<SoftBody> > _softbodies, Integrator _time_integrator);
	virtual ~Solver();
	void step(double h);
	void reset();
	
	Integrator time_integrator;

private:
	std::vector< std::shared_ptr<SoftBody> > softbodies;
	Eigen::MatrixXd A;
	Eigen::MatrixXd M;
	Eigen::MatrixXd J;
	Eigen::VectorXd x;
	Eigen::VectorXd b;
	Eigen::VectorXd f;
	Eigen::VectorXd v;
	Eigen::Vector3d grav;
};

