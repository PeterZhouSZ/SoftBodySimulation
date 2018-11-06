#pragma once
#ifndef __Scene__
#define __Scene__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Node;
class MatrixStack;
class Program;
class Shape;
class Solver;
class SoftBody;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
	void toggleElasticForce();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
	
	double getTime() const { return t; }
	std::vector<std::shared_ptr<SoftBody> > softbodies;
	std::shared_ptr<Solver> solver;
	
private:
	bool isElasticForce;
	bool isGravity;
	double t;
	double h;
	Eigen::Vector3d grav;
	

};

#endif
