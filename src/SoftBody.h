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

#include "Tetrahedron.h"

class SoftBody
{
public:
	SoftBody();
	SoftBody(double _young, double _poisson, Material _material, std::string mesh_name);
	virtual ~SoftBody() {}

	double getStiffness() { return this->young; }
	double getPoisson() { return this->poisson; }
	void setStiffness(double young);
	void setPoisson(double poisson);
	void step(double h, const Eigen::Vector3d &grav);
	void init();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void updatePosNor();
	void fixPointsByYvalue(double y);
	void flatten(double y);
	void computeStiffness(Eigen::MatrixXd &K);
	void computeInvertibleStiffness(Eigen::MatrixXd &K);

	void computeGravityForce(Eigen::Vector3d grav, Eigen::VectorXd &f);
	void computeElasticForce(Eigen::VectorXd &f);
	void computeInvertibleElasticForce(Eigen::VectorXd &f);
	int getNumNodes() const { return nodes.size(); }
	int getNumTets() const { return tets.size(); }
	int getNumTrifaces() const { return trifaces.size(); }
	std::vector<std::shared_ptr<Tetrahedron> > getTets() const { return tets; };
	std::vector<std::shared_ptr<TriFace> > getTrifaces() const { return trifaces; }
	std::vector<std::shared_ptr<Node> > getNodes() const { return nodes; }

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