#define TETLIBRARY
#include <iostream>
#include <tetgen.h>
#include <cmath>        // std::abs

#include "SoftBody.h"
#include "Tetrahedron.h"
#include "Particle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "stdlib.h"
#define CRTDBG_MAP_ALLOC

#include <crtdbg.h>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;

SoftBody::SoftBody() {
	tetgenio input_mesh, output_mesh;
	input_mesh.load_ply("dodecahedron");
	tetrahedralize("pqz", &input_mesh, &output_mesh);

}



void SoftBody::step(double h, const Eigen::Vector3d &grav) {


}


void SoftBody::init() {


}


void SoftBody::tare() {


}


void SoftBody::reset() {


}


void SoftBody::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const {



}

SoftBody::~SoftBody() {

}