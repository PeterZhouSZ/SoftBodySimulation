#include <iostream>

#include "Node.h"
#include "TriFace.h"

using namespace std;
using namespace Eigen;

TriFace::TriFace() {
}

TriFace::TriFace(std::vector<shared_ptr<Node>> _nodes):
nodes(_nodes)
{
}

Eigen::Vector3d TriFace::computeNormal() {
	Vector3d p0 = nodes[0]->x;
	Vector3d p1 = nodes[1]->x;
	Vector3d p2 = nodes[2]->x;

	Vector3d e0 = p1 - p0;
	Vector3d e1 = p2 - p0;
	Vector3d normal = e0.cross(e1);
	normal.normalize();
	this->normal = normal;

	return normal;
}


TriFace::~TriFace() {
}