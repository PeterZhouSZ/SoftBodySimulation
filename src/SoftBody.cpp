#define TETLIBRARY
#include <iostream>
#include <tetgen.h>
#include <cmath>        // std::abs

#include "SoftBody.h"
#include "Tetrahedron.h"
#include "TriFace.h"
#include "Node.h"
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

}

SoftBody::SoftBody(double _young, double _poisson, Material _material):
young(_young),
poisson(_poisson),
mt(_material)
{
	tetgenio input_mesh, output_mesh;
	input_mesh.load_ply("dodecahedron");
	tetrahedralize("pqz", &input_mesh, &output_mesh);
	
	nFacets = output_mesh.numberoffacets;
	double r = 0.02;

	// Create Nodes
	for (int i = 0; i < output_mesh.numberofpoints; i++) {
		auto node = make_shared<Node>();
		nodes.push_back(node);
		node->r = r;
		node->x0 << output_mesh.pointlist[3 * i + 0],
			output_mesh.pointlist[3 * i + 1],
			output_mesh.pointlist[3 * i + 2];

		node->x = node->x0;
		node->v0.setZero();
		node->v = node->v0;
		node->m = 0.0;
		node->i = i;
		node->fixed = false;
	}

	// Create TriFaces
	for (int i = 0; i < output_mesh.numberoftrifaces; i++) {
		auto triface = make_shared<TriFace>();
		trifaces.push_back(triface);
		for (int ii = 0; ii < 3; ii++) {
			triface->nodes.push_back(nodes[output_mesh.trifacelist[3 * i + ii]]);
		}
	}

	// Create Tets
	vector<shared_ptr<Node>> tet_nodes;
	for (int i = 0; i < output_mesh.numberoftetrahedra; i++) {
		tet_nodes.clear();
		for (int ii = 0; ii < 4; ii++) {
			tet_nodes.push_back(nodes[output_mesh.tetrahedronlist[4 * i + ii]]);
		}

		auto tet = make_shared<Tetrahedron>(young, poisson, mt, tet_nodes);
		tets.push_back(tet);	
	}

	// Init Buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(trifaces.size() * 9);
	norBuf.resize(trifaces.size() * 9);
	eleBuf.resize(trifaces.size() * 3);
	updatePosNor();

	for (int i = 0; i < trifaces.size(); i++) {
		eleBuf[3 * i + 0] = 3 * i;
		eleBuf[3 * i + 1] = 3 * i + 1;
		eleBuf[3 * i + 2] = 3 * i + 2;
	}
}

void SoftBody::step(double h, const Eigen::Vector3d &grav) {


}


void SoftBody::init() {
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	/*glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size() * sizeof(float), &texBuf[0], GL_STATIC_DRAW);*/

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	assert(glGetError() == GL_NO_ERROR);
}

void SoftBody::updatePosNor() {
	for (int i = 0; i < trifaces.size(); i++) {
		auto triface = trifaces[i];

		Vector3d p0 = triface->nodes[0]->x;
		Vector3d p1 = triface->nodes[1]->x;
		Vector3d p2 = triface->nodes[2]->x;

		Vector3d normal = triface->computeNormal();

		for (int ii = 0; ii < 3; ii++) {
			posBuf[9 * i + 0 + ii] = p0(ii);
			posBuf[9 * i + 3 + ii] = p1(ii);
			posBuf[9 * i + 6 + ii] = p2(ii);

			norBuf[9 * i + 0 + ii] = normal(ii);
			norBuf[9 * i + 3 + ii] = normal(ii);
			norBuf[9 * i + 6 + ii] = normal(ii);
		}
	}
}

void SoftBody::tare() {

}


void SoftBody::reset() {

}


void SoftBody::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const {
	// Draw mesh
	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(p->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());
	MV->pushMatrix();

	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);


	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);


	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);

	glDrawElements(GL_TRIANGLES, 3 * trifaces.size(), GL_UNSIGNED_INT, (const void *)(0 * sizeof(unsigned int)));

	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();


}

SoftBody::~SoftBody() {

}