#include <iostream>

#include "Scene.h"
#include "Node.h"
#include "Shape.h"
#include "Program.h"
#include "Solver.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-2),
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
	double stiffness = 1e2;
	Vector2d damping(0.0, 1.0);


	sphereShape = make_shared<Shape>();
	sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	
	auto sphere = make_shared<Node>(sphereShape);
	spheres.push_back(sphere);
	sphere->r = 0.1;
	sphere->x = Vector3d(0.0, 0.2, 0.0);
}

void Scene::init()
{
	sphereShape->init();
	
}

void Scene::tare()
{
	for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->tare();
	}
	
}

void Scene::reset()
{
	t = 0.0;
	for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->reset();
	}
	
}

void Scene::step()
{
	t += h;
	
	
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	for(int i = 0; i < (int)spheres.size(); ++i) {
		//spheres[i]->draw(MV, prog);
	}
}
