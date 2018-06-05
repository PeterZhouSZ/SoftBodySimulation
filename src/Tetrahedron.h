#pragma once
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Node;
enum Material { LINEAR, NEOHOOKEAN, STVK, COROTATED_LINEAR };


class Tetrahedron
{
public:
	Tetrahedron(double _young, double _poisson, Material _material);
	virtual ~Tetrahedron();

	void step(double h, const Eigen::Vector3d &grav);

	void tare();
	void reset();
	void precomputation();
	Eigen::Matrix3d computeDm(std::vector<std::shared_ptr<Node>> nodes);

	void computeElasticForces();
	void computeForceDifferentials();

private:
	std::vector<std::shared_ptr<Node>> nodes;	// i, j, k, l
	double young;
	double poisson;
	double mass;
	Material material;

	// precomputed
	double W;	// undeformed volume of Te
	Eigen::Matrix3d Dm;		// reference shape matrix ("material-space" shape matrix) is constant 
	Eigen::Matrix3d Bm;		// Dm.inv()
	
	// updated each step
	Eigen::Matrix3d Ds;		// deformed shape matrix
	Eigen::Matrix3d F;		// deformation gradient
	Eigen::Matrix3d P;		// Piola stress
	Eigen::Matrix3d H;		// forces matrix


	// for force differentials
	Eigen::Matrix3d dDs;
	Eigen::Matrix3d dF;		// the differential of the deformation gradient
	Eigen::Matrix3d dP;		// the stress differential 
	Eigen::Matrix3d dH;		// the nodal force differential of the first three vertices


};