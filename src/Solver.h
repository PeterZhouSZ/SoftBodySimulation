#pragma once

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>


class MatrixReplacement;
class SoftBody;
enum Integrator { RKF45, SYMPLECTIC, IMPLICIT };
typedef Eigen::Triplet<double> T;

class Solver
{
public:
	Solver(std::vector< std::shared_ptr<SoftBody> > _softbodies, Integrator _time_integrator, Eigen::Vector2d _damping, Eigen::Vector3d _grav, bool _isSparse, bool _isMatrixFree);
	virtual ~Solver();
	void step(double h);
	void reset();
	void clearVelocity();
	
	Integrator time_integrator;
	double y_floor;
	bool isSparse;
	bool isMatrixFree;
	bool isGravity;
	bool isElasticForce;
	bool isMosek;
	bool isFloor;
	Eigen::Vector3d nor_floor;

private:
	
	std::vector< std::shared_ptr<SoftBody> > softbodies;
	Eigen::SparseMatrix<double> A_sparse;
	Eigen::SparseMatrix<double> G_sparse;
	Eigen::MatrixXd A;
	Eigen::MatrixXd K;
	Eigen::MatrixXd Dx;
	Eigen::MatrixXd M;
	Eigen::MatrixXd J;
	Eigen::VectorXd x;
	Eigen::VectorXd b;
	Eigen::VectorXd f;
	Eigen::VectorXd v;
	Eigen::Vector3d grav;
	Eigen::Vector2d damping;
	std::vector<T> A_;
	std::vector<T> G_;
	
	int mat_n;
};

