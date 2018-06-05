#pragma once

#include <vector>
#include <memory>
#include <tetgen.h>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Node;
class MatrixStack;
class Program;
class Tetrahedron;
class TriFace;

class SoftBody
{
public:
	SoftBody();
	SoftBody(double _young, double _poisson, Material _material);
	virtual ~SoftBody();

	void step(double h, const Eigen::Vector3d &grav);
	void init();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void updatePosNor();
private:
	std::vector<std::shared_ptr<Node> > nodes;
	std::vector<std::shared_ptr<TriFace> > trifaces;
	std::vector<std::shared_ptr<Tetrahedron> > tets;
	
	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;

	int nFacets;
	double young;
	double poisson;
	Material mt;
};