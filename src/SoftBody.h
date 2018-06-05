#pragma once

#include <vector>
#include <memory>
#include <tetgen.h>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Particle;
class MatrixStack;
class Program;
class Tetrahedron;

class SoftBody
{
public:
	SoftBody();
	virtual ~SoftBody();

	void step(double h, const Eigen::Vector3d &grav);
	void init();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;

private:
	std::vector<std::shared_ptr<Tetrahedron> > tets;
	
	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;

};