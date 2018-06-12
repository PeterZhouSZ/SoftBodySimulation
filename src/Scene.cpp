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
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 1e-2;
	grav << 0.0, -9.8, 0.0;
	
	int rows = 2;
	int cols = 2;
	double mass = 0.1;
	double stiffness = 1e4;
	Vector2d damping(0.0, 1.0);

	auto softbody = make_shared<SoftBody>(1e3, 0.4, STVK);
	softbodies.push_back(softbody);
	solver = make_shared<Solver>(softbodies, SYMPLECTIC);
}

void Scene::init()
{
	//sphereShape->init();
	for (int i = 0; i < (int)softbodies.size(); ++i) {
		softbodies[i]->init();
	}
}

void Scene::tare()
{
	/*for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->tare();
	}*/
	
}

void Scene::reset()
{
	/*t = 0.0;
	for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->reset();
	}*/
	
}

void Scene::step()
{
	t += h;
	solver->step(h);

	for (int i = 0; i < (int)softbodies.size(); ++i) {
		softbodies[i]->step(h, grav);
	}
	
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	//for(int i = 0; i < (int)spheres.size(); ++i) {
	//	//spheres[i]->draw(MV, prog);
	//}
	for (int i = 0; i < (int)softbodies.size(); ++i) {
		softbodies[i]->draw(MV, prog);
	}

}
