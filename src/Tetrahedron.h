#pragma once
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "MLCommon.h"

class Node;
enum Material { LINEAR, NEOHOOKEAN, STVK, COROTATED_LINEAR };


class Tetrahedron
{
public:
	Tetrahedron(double _young, double _poisson, double _density, Material _material, std::vector<std::shared_ptr<Node>> _nodes);
	virtual ~Tetrahedron();

	void step(double h, const Eigen::Vector3d &grav);

	void tare();
	void reset();
	void precomputation();
	void setStiffness(double young);
	void setPoisson(double poisson);
	Eigen::Matrix3d computeDeformationGradient();
	Matrix3x4d computeAreaWeightedVertexNormals();
	Eigen::Matrix3d computePKStress(Eigen::Matrix3d F, Material mt, double mu, double lambda);
	Eigen::Matrix3d computePKStressDerivative(Eigen::Matrix3d F, Eigen::Matrix3d dF, Material mt, double mu, double lambda);
	void computeElasticForces(Eigen::VectorXd f);
	void computeInvertibleElasticForces(Eigen::VectorXd &f);
	Eigen::Matrix3d computeInvertiblePKStress(Eigen::Matrix3d F, double mu, double lambda);
	void computeForceDifferentials(Eigen::VectorXd dx, Eigen::VectorXd& df);
	std::vector<std::shared_ptr<Node>> nodes;	// i, j, k, l
	Eigen::MatrixXd getStiffness() const { return this->K; }

	bool isInvert;
	int clamped;

private:
	Matrix3x4d Nm;		// area-weighted vertex normals [Irving 04]	(Bm in paper)		

	double young;
	double poisson;
	// Lame coefficients
	double mu;
	double lambda;
	double mass;
	double density;
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
	
	// SVD 
	Eigen::Matrix3d U;
	Eigen::Matrix3d V;
	Eigen::Matrix3d Fhat;	// diagonalized deformation gradient
	Eigen::Matrix3d Phat;	
	

	// for force differentials
	Eigen::Matrix3d dDs;
	Eigen::Matrix3d dF;		// the differential of the deformation gradient
	Eigen::Matrix3d dP;		// the stress differential 
	Eigen::Matrix3d dH;		// the nodal force differential of the first three vertices
	Matrix12d K;		// 12x12 stiffness matrix
};