#include <iostream>

#include "Scene.h"
#include "Node.h"
#include "Shape.h"
#include "SoftBody.h"
#include "Program.h"
#include "Solver.h"
#include "Tetrahedron.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-3),
	grav(0.0, 0.0, 0.0)
{
	isElasticForce = true;
	isGravity = true;
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 1e-4;
	grav << 0.0, -9.8, 0.0;
	
	double mass = 1;
	double stiffness = 1.0e4;
	double possion = 0.35;
	bool isSparse = false;
	bool isMatrixFree = false;
	Vector2d damping(0.0, -0.1);
	double y_floor = -10.0;
	Vector3d nor_floor;
	nor_floor << 0.0, -1.0, 0.0;
	auto softbody = make_shared<SoftBody>(stiffness, possion, NEOHOOKEAN, "bunny100");
	softbodies.push_back(softbody);
	softbody->fixPointsByYvalue(0.5);
	solver = make_shared<Solver>(softbodies, SYMPLECTIC, damping, grav, isSparse, isMatrixFree);
	solver->y_floor = y_floor;
	solver->nor_floor = nor_floor;
	solver->isMosek = false;
	solver->isElasticForce = true;
	solver->isGravity = true;
}

void Scene::init()
{
	for (int i = 0; i < (int)softbodies.size(); ++i) {
		softbodies[i]->init();
	}
}

void Scene::tare()
{

}

void Scene::reset()
{
}

void Scene::toggleElasticForce() {
	isElasticForce = !isElasticForce;
	if (isElasticForce) {
		cout << "Elastic Force is On" << endl;
	}
	else {
		cout << "Elastic Force is Off " << endl;
	}
}

void Scene::toggleGravForce() {
	isGravity = !isGravity;
	if (isGravity) {
		cout << "Gravity Force is On" << endl;
	}
	else {
		cout << "Gravity Force is Off " << endl;
	}
}

void Scene::flattenSoftBody(double y) {
	softbodies[0]->flatten(y);


}

void Scene::step()
{
	t += h;

	solver->isElasticForce = isElasticForce;
	solver->isGravity = isGravity;

	solver->step(h);

	for (int i = 0; i < (int)softbodies.size(); ++i) {
		softbodies[i]->step(h, grav);
	}
	
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	
	for (int i = 0; i < (int)softbodies.size(); ++i) {
		softbodies[i]->draw(MV, prog);
	}
}
