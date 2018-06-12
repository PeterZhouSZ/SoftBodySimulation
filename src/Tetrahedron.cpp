#include "Tetrahedron.h"
#include "Node.h"

using namespace Eigen;
using namespace std;

Tetrahedron::Tetrahedron(double _young, double _poisson, double _density, Material _material, vector<shared_ptr<Node>> _nodes):
young(_young),
poisson(_poisson),
density(_density),
material(_material),
nodes(_nodes)
{
	this->mu = young / (2.0 * (1.0 + poisson));
	this->lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));

	precomputation();
	K.resize(12, 12);
	K.setZero();
}

void Tetrahedron::step(double h, const Eigen::Vector3d &grav) {
	computeElasticForces();
	computeForceDifferentials();
}

void Tetrahedron::precomputation() {

	assert(nodes.size() == 4);
	for (int i = 0; i < nodes.size() - 1; i++) {
		this->Dm.col(i) = nodes[i]->x0 - nodes[3]->x0;
	}
	this->Bm = this->Dm.inverse();
	this->W = 1.0 / 6.0 * Dm.determinant();
	this->mass = this->W * this->density;

	// Distribute 1/4 mass to each node
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->m += this->mass * 0.25;
	}
}

void Tetrahedron::computeElasticForces() {
	assert(nodes.size() == 4);
	for (int i = 0; i < nodes.size() - 1; i++) {
		this->Ds.col(i) = nodes[i]->x - nodes[3]->x;
	}

	this->F = Ds * Bm;
	this->P = computePKStress(F, material, mu, lambda);
	this->H = -W * P * (Bm.transpose());

	for (int i = 0; i < nodes.size() - 1; i++) {
		nodes[i]->addForce(H.col(i));
		nodes[3]->addForce(-H.col(i));
	}
}

Matrix3d Tetrahedron::computePKStress(Matrix3d F, Material mt, double mu, double lambda) {
	Matrix3d E, P, I;
	E.setZero();
	P.setZero();
	I.setIdentity();
	double psi;
	
	switch (mt)
	{
		case LINEAR:
		{
			E = 0.5 * (F + F.transpose()) - I;
			//psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
			P = 2.0 * mu * E + lambda * E.trace() * I;
			break;
		}
		
		case NEOHOOKEAN: 
		{
			double I1 = (F.transpose() * F).trace();
			double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
			double I3 = (F.transpose() * F).determinant();
			double J = sqrt(I3);
			//psi = 1.0 / 2.0 * mu *(I1 - 3.0) - mu * log(J) + 1.0 / 2.0 * lambda * log(J)*log(J);
			P = mu * (F - F.inverse().transpose()) + lambda * log(J)*(F.inverse().transpose());
			break;
		}
		
		case STVK: 
		{
			E = 0.5 * (F.transpose() * F - I);
			//psi = mu * E.norm()*E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
			P = F * (2.0 * mu * E + lambda * E.trace() * I);
			break;
		}
		
		case COROTATED_LINEAR: 
		{
			// Polar decomposition
			Matrix3d A = F.adjoint() * F;
			SelfAdjointEigenSolver<Matrix3d> es(A);
			Matrix3d S = es.operatorSqrt();
			Matrix3d R = F * S.inverse();

			E = S - I;
			//psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
			//P = R * (2.0 * mu * E + lambda * E.trace() * I);
			P = 2.0 * mu * (F - R) + lambda * (R.transpose() * F - I).trace() * R;
			break; 
		}

		default: 
		{
			break;
		}	
	}

	return P;
}

Matrix3d Tetrahedron::computePKStressDerivative(Matrix3d F, Matrix3d dF, Material mt, double mu, double lambda) {
	Matrix3d E, dE, P, dP, I;
	E.setZero();
	dE.setZero();
	P.setZero();
	dP.setZero();
	I.setIdentity();

	switch (mt) {
	case COROTATED_LINEAR:
	{
		Matrix3d A = F.adjoint() * F;
		SelfAdjointEigenSolver<Matrix3d> es(A);
		Matrix3d S = es.operatorSqrt();
		Matrix3d R = F * S.inverse();
		E = S - I;
		P = 2.0 * mu *(F - R) + lambda * (R.transpose()*F - I).trace() * R;
		break;
	}

	case STVK:
	{
		E = 1.0 / 2.0 * (F.transpose() * F - I);
		dE = 1.0 / 2.0 * (dF.transpose() * F + F.transpose() * dF);
		P = F * (2.0 * mu * E + lambda * E.trace() * I);
		dP = dF * (2.0 * mu * E + lambda * E.trace() * I) + F * (2.0 * mu * dE + lambda * dE.trace() * I);
		break;
	}

	case NEOHOOKEAN:
	{
		//double I1 = (F.norm()) * (F.norm());
		//double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		MatrixXd FT = F.transpose();
		MatrixXd FIT = F.inverse().transpose();
		//double I3 = (F.transpose() * F).determinant();
		double I3 = (FT * F).determinant();
		double J = sqrt(I3);
		P = mu * (F - FIT) + lambda * log(J) * FIT;
		dP = mu * dF + (mu - lambda * log(J)) * FIT * (dF.transpose()) * FIT + lambda * ((F.inverse() * dF)).trace() * FIT;

		//P = mu * (F - (F.inverse().transpose())) + lambda * log(J) * (F.inverse().transpose());
		//dP = mu * dF + (mu - lambda * log(J)) * (F.inverse().transpose()) * (dF.transpose()) * (F.inverse().transpose()) + lambda * ((F.inverse() * dF)).trace() * (F.inverse().transpose());
		break;
	}
	case LINEAR:
	{
		E = 1.0 / 2.0 * (F + F.transpose()) - I;
		dE = 1.0 / 2.0 * (dF + dF.transpose());
		P = 2.0 * mu * E + lambda * E.trace() * I;
		dP = 2.0 * mu * dE + lambda * dE.trace() * I;
		break;
	}
	default:
		break;
	}

	return dP;
}

void Tetrahedron::computeForceDifferentials() {
	assert(nodes.size() == 4);
	for (int i = 0; i < nodes.size() - 1; i++) {
		this->Ds.col(i) = nodes[i]->x - nodes[3]->x;
	}

	this->F = Ds * Bm;

	MatrixXd dFRow(4, 3);
	for (int i = 0; i < 3; ++i) {
		dFRow.row(i) = Bm.row(i);
		dFRow(3, i) = -Bm(0, i) - Bm(1, i) - Bm(2, i);
	}

	K.setZero();
	for (int row = 0; row < 4; ++row) {
		MatrixXd Kb(12, 3);
		Kb.setZero();
		for (int kk = 0; kk < 3; ++kk) {
			Matrix3d dF;
			dF.setZero();
			dF.row(kk) = dFRow.row(row);
			MatrixXd dP = computePKStressDerivative(F, dF, material, mu, lambda);
			dH = -W * dP * (Bm.transpose());
			for (int ii = 0; ii < 3; ii++) {
				for (int ll = 0; ll < 3; ll++) {
					Kb(ii * 3 + ll, kk) = dH(ll, ii);
				}
				Kb(9 + ii, kk) = -dH(ii, 0) - dH(ii, 1) - dH(ii, 2);
			}
		}
		K.block<12, 3>(0, 3 * row) = Kb;
	}
}


void Tetrahedron::tare() {


}

void Tetrahedron::reset() {


}


Tetrahedron:: ~Tetrahedron() {

}