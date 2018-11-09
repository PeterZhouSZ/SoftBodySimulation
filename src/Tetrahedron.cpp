#include "Tetrahedron.h"
#include "Node.h"
#include <iostream>

using namespace Eigen;
using namespace std;
#define Fthreshold 0.2

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
	K.setZero();
}

void Tetrahedron::setStiffness(double young) {
	this->mu = young / (2.0 * (1.0 + poisson));
	this->lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));

}


void Tetrahedron::setPoisson(double poisson) {
	this->mu = young / (2.0 * (1.0 + poisson));
	this->lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));

}

Eigen::Matrix3d Tetrahedron::computeDeformationGradient() {

	for (int i = 0; i < (int)nodes.size() - 1; i++) {
		this->Ds.col(i) = nodes[i]->x - nodes[3]->x;
	}

	this->F = Ds * Bm;
	return this->F;
}

Matrix3x4d Tetrahedron::computeAreaWeightedVertexNormals() {
	Vector3d va, vb, vc, vd;
	va = nodes[0]->x;
	vb = nodes[1]->x;
	vc = nodes[2]->x;
	vd = nodes[3]->x;

	// Computes normals for the four faces: acb, adc, abd, bcd
	Vector3d acb_normal, adc_normal, abd_normal, bcd_normal;
	acb_normal = (vc - va).cross(vb - va);
	adc_normal = (vd - va).cross(vc - va);
	abd_normal = (vb - va).cross(vd - va);
	bcd_normal = (vc - vb).cross(vd - vb);

	// if the tet vertices abcd form a positive orientation, no need to correct
	// if not, flip them
	double orientation = (vd - va).dot((vb - va).cross(vc - va));
	if (orientation < 0.0) {
		acb_normal *= -1.0;
		adc_normal *= -1.0;
		abd_normal *= -1.0;
		bcd_normal *= -1.0;
	}

	// Computes the area of triangles
	// area = 0.5 * | u x v |
	double acb_area, adc_area, abd_area, bcd_area;
	acb_area = 0.5 * sqrt(acb_normal.dot(acb_normal));
	adc_area = 0.5 * sqrt(adc_normal.dot(adc_normal));
	abd_area = 0.5 * sqrt(abd_normal.dot(abd_normal));
	bcd_area = 0.5 * sqrt(bcd_normal.dot(bcd_normal));

	acb_normal.normalize();
	adc_normal.normalize();
	abd_normal.normalize();
	bcd_normal.normalize();

	/*cout << "acb nor" << endl << acb_normal << endl;
	cout << "adc nor" << endl << adc_normal << endl;
	cout << "abd nor" << endl << abd_normal << endl;
	cout << "bcd nor" << endl << bcd_normal << endl;*/


	this->Nm.col(0) = -(acb_area * acb_normal + adc_area * adc_normal + abd_area * abd_normal) / 3.0;
	this->Nm.col(1) = -(acb_area * acb_normal + abd_area * abd_normal + bcd_area * bcd_normal) / 3.0;
	this->Nm.col(2) = -(acb_area * acb_normal + adc_area * adc_normal + bcd_area * bcd_normal) / 3.0;
	this->Nm.col(3) = -(adc_area * adc_normal + abd_area * abd_normal + bcd_area * bcd_normal) / 3.0;

	return this->Nm;
}

void Tetrahedron::step(double h, const Eigen::Vector3d &grav) {
	//computeElasticForces();
	//computeForceDifferentials();
}

void Tetrahedron::precomputation() {

	for (int i = 0; i < nodes.size() - 1; i++) {
		this->Dm.col(i) = nodes[i]->x0 - nodes[3]->x0;
	}

	this->Bm = this->Dm.inverse();
	this->W = abs(1.0 / 6.0 * Dm.determinant()); // probably negative 

	this->mass = this->W * this->density;
	
	// Distribute 1/4 mass to each node
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->m += this->mass * 0.25;
	}

	computeAreaWeightedVertexNormals();
}

void Tetrahedron::computeElasticForces(Eigen::VectorXd &f) {

	this->F = computeDeformationGradient();
	this->P = computePKStress(F, material, mu, lambda);
	this->H = -W * P * (Bm.transpose());

	for (int i = 0; i < nodes.size() - 1; i++) {
		nodes[i]->addForce(H.col(i));
		nodes[3]->addForce(-H.col(i));
		f.segment<3>(3 * nodes[i]->i) += H.col(i);
		f.segment<3>(3 * nodes[3]->i) -= H.col(i);
	}
}

void Tetrahedron::computeInvertibleElasticForces(VectorXd &f) {

	this->F = computeDeformationGradient();
	// The deformation gradient is available in this->F
	if (this->F.determinant() < 0.0) {
		isInvert = true;
	}
	else {
		isInvert = false;
	}

	// SVD on the deformation gradient
	int modifiedSVD = 1;
	Vector3d Fhat_vec;

	if (!SVD(this->F, this->U, Fhat_vec, this->V, 1e-8, modifiedSVD)) {
		//cout << "error in svd " << endl;
	}
	this->Fhat = Fhat_vec.asDiagonal();
	cout << "Fhat:" << endl << this->Fhat << endl;

	// Test if correct. checked
	//cout << "F" << F << endl;
	//cout << "Fhat" << Fhat << endl;
	//cout << "UFV'" << U * Fhat * V.transpose() << endl;

	// SVD result is available in this->U, this->V, Fhat_vec, this->Fhat

	// clamp if below the principal stretch threshold
	clamped = 0;
	for (int i = 0; i < 3; i++)
	{
		if (this->Fhat(i, i) < Fthreshold)
		{
			this->Fhat(i, i) = Fthreshold;
			clamped |= (1 << i);
		}
	}
	cout << "Fhat:" << endl << this->Fhat << endl;
	//clamped = 0; // disable clamping

	// Computes the internal forces
	// Computes P first and computes the nodal forces G=PBm in section 4 of [Irving 04]

	// Computes the diagonal P tensor

	this->Phat = computeInvertiblePKStress(this->Fhat, mu, lambda);
	cout << "Phat:" << endl << this->Phat << endl;
	// The result is the same as 
	// MatrixXd temp = computePKStress(this->Fhat, material, mu, lambda);


	// P = U * diag(Phat) * V'

	MatrixXd temp = this->U.transpose() * this->P * this->V;
	this->P = this->U * this->Phat * this->V.transpose();
	
	//cout << "P:"<< endl<< this->P << endl;

	// Computes the nodal forces by G=PBm=PNm
	//cout << "forces: " << endl << this->H << endl;
	cout << "instead: " << endl << this->P * this->Nm.block<3, 3>(0, 0) << endl;

	for (int i = 0; i < (int)nodes.size()-1; i++) {
		f.segment<3>(3 * nodes[i]->i) += this->P * this->Nm.col(i);
		f.segment<3>(3 * nodes[3]->i) -= this->P * this->Nm.col(i);

	}

}


Matrix3d Tetrahedron::computeInvertiblePKStress(Matrix3d F, double mu, double lambda) {
	// Neohookean
	Vector3d invariants;
	Vector3d lambda1, lambda2;
	lambda1 << F(0, 0), F(1, 1), F(2, 2);

	lambda2 << lambda1(0) * lambda1(0), lambda1(1) * lambda1(1), lambda1(2) * lambda1(2);
	double IC = lambda2(0) + lambda2(1) + lambda2(2);
	double IIC = lambda2(0) * lambda2(0) + lambda2(1) * lambda2(1) + lambda2(2) * lambda2(2);
	double IIIC = lambda2(0) * lambda2(1) * lambda2(2);

	invariants << IC, IIC, IIIC;

	Vector3d dPsidIV;
	dPsidIV << 0.5 * mu, 0.0, (-0.5 * mu + 0.25 * lambda * log(IIIC)) / IIIC;

	// PDiag = [ dI / dlambda ]^T * dPsidI
	Matrix3d matM;
	matM << 2.0 * lambda1(0), 2.0 * lambda1(1), 2.0 * lambda1(2),
		4.0 * lambda1(0) * lambda1(0) * lambda1(0), 4.0 * lambda1(1) * lambda1(1) * lambda1(1), 4.0 * lambda1(2) * lambda1(2) * lambda1(2),
		2.0 * lambda1(0) * lambda2(1) * lambda2(2), 2.0 * lambda1(1) * lambda2(0) * lambda2(2), 2.0 * lambda1(2) * lambda2(0) * lambda2(1);

	Vector3d result;
	result = matM.transpose() * dPsidIV;
	this->Phat = result.asDiagonal();
	return this->Phat;

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

void Tetrahedron::computeForceDifferentials(VectorXd dx, VectorXd& df) {

	this->F = computeDeformationGradient();
	/*for (int i = 0; i < (int)nodes.size(); i++) {
		
		for (int ii = 0; ii < (int)nodes.size() - 1; i++) {
			Vector12d diff;
			diff.setZero();
			diff(4 * i + ii) = 1.0;

			this->dDs.col(ii) = diff.segment<3>(3 * ii) - dx.segment<3>(9);

	}*/


	for (int i = 0; i < nodes.size() - 1; i++) {
		this->dDs.col(i) = dx.segment<3>(3 * nodes[i]->i) - dx.segment<3>(3 * nodes[3]->i);
	}

	this->dF = dDs * Bm;
	this->dP = computePKStressDerivative(F, dF, material, mu, lambda);
	this->dH = -W * dP * (Bm.transpose());

	for (int i = 0; i < nodes.size() - 1; i++) {
		df.segment<3>(3 * nodes[i]->i) += this->dH.col(i);
		df.segment<3>(3 * nodes[3]->i) -= this->dH.col(i);
	}

	/*MatrixXd dFRow(4, 3);
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
	}*/
}
