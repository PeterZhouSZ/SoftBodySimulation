#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include "Tetrahedron.h"
#include "SoftBody.h"

class MatrixReplacement;
using Eigen::SparseMatrix;

namespace Eigen {
	namespace internal {
		template<>
		struct traits<MatrixReplacement> : public
			Eigen::internal::traits<Eigen::SparseMatrix<double> >
		{};
	}
}

class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
	// Required typedefs, constants, and method:
	typedef double Scalar;
	typedef double RealScalar;
	typedef int StorageIndex;

	enum {
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic,
		IsRowMajor = false
	};

	Index rows() const { return M.rows(); }
	Index cols() const { return M.cols(); }

	template<typename Rhs>
	Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
		return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
	}

	// Custon API:
	MatrixReplacement(std::vector< std::shared_ptr<SoftBody> > _softbodies, Eigen::MatrixXd _M) : softbodies(_softbodies), M(_M) {}
	void getCol(int col) {





	}

private:
	std::vector< std::shared_ptr<SoftBody> > softbodies;
	Eigen::MatrixXd M;

};