#pragma once
#ifndef __NODE__
#define __NODE__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;

class Node
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Node();
	Node(const std::shared_ptr<Shape> shape);
	virtual ~Node();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void clearForce();
	
	double r; // radius
	double m; // mass
	int i;  // starting index

	bool fixed;
	Eigen::Vector3d x0; // initial position
	Eigen::Vector3d v0; // initial velocity
	Eigen::Vector3d x;  // position
	Eigen::Vector3d v;  // velocity
	Eigen::Vector3d f;
	Eigen::Vector3d df;
private:
	const std::shared_ptr<Shape> sphere;

	
};

#endif
