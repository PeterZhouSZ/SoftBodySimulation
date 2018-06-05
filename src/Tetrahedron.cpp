#include "Tetrahedron.h"
#include "Node.h"

using namespace Eigen;
using namespace std;

Tetrahedron::Tetrahedron(double _young, double _poisson, Material _material):
young(_young),
poisson(_poisson),
material(_material)
{


}

void Tetrahedron::step(double h, const Eigen::Vector3d &grav) {


}

void Tetrahedron::precomputation() {
	this->Dm = computeDm(this->nodes);
	this->Bm = this->Dm.inverse();
	this->W = 1.0 / 6.0 * Dm.determinant();
}

Matrix3d Tetrahedron::computeDm(vector<shared_ptr<Node>> nodes) {
	Matrix3d Dm;
	assert(nodes.size() == 4);
	for (int i = 0; i < nodes.size(); i++) {
		Dm.col(i) = nodes[i]->x0 - nodes[3]->x0;
	}
	return Dm;
}


void Tetrahedron::computeElasticForces() {


}


void Tetrahedron::computeForceDifferentials() {



}


void Tetrahedron::tare() {


}

void Tetrahedron::reset() {


}


Tetrahedron:: ~Tetrahedron() {

}