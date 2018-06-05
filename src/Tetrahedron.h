#pragma once
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Particle;

class Tetrahedron
{
public:
	Tetrahedron();
	virtual ~Tetrahedron();

	void step(double h, const Eigen::Vector3d &grav);

	void init();
	void tare();
	void reset();


private:
	std::vector<std::shared_ptr<Particle>> nodes;
	double young;
	double poisson;




};