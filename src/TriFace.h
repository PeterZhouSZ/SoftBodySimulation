#pragma once
#ifndef __TRIFACE__
#define __TRIFACE__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Node;

class TriFace
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	TriFace();
	TriFace(std::vector<std::shared_ptr<Node>> _nodes);
	virtual ~TriFace();
	Eigen::Vector3d computeNormal();
	
	std::vector<std::shared_ptr<Node>> nodes;
	Eigen::Vector3d normal;

};

#endif
