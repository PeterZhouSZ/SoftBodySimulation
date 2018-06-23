//#pragma once
//#include <iostream>
//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include <Eigen/IterativeLinearSolvers>
//#include <unsupported/Eigen/IterativeSolvers>
//
//class MatrixReplacement;
//using Eigen::SparseMatrix;
//
//namespace Eigen {
//	namespace internal {
//		template<>
//		struct traits<MatrixReplacement> : public
//			Eigen::internal::traits<Eigen::SparseMatrix<double> >
//		{};
//	}
//}
//
//class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
//public:
//	// Required typedefs, constants, and method:
//	typedef double Scalar;
//	typedef double RealScalar;
//	typedef int StorageIndex;
//
//	enum {
//		ColsAtCompileTime = Eigen::Dynamic,
//		MaxColsAtCompileTime = Eigen::Dynamic,
//		IsRowMajor = false
//	};
//
//	Index rows() const { return mp_mat->rows(); }
//	Index cols() const { return mp_mat->cols(); }
//
//	template<typename Rhs>
//	Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
//		return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
//	}
//
//	// Custon API:
//	MatrixReplacement() : mp_mat(0) {}
//
//private:
//
//};