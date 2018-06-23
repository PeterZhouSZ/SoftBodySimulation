//#include "MatrixFree.h"
//
//namespace Eigen {
//	namespace internal {
//
//		template<typename Rhs>
//		struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
//			: generic_product_impl_base<MatrixReplacement, Rhs, generic_product_impl<MatrixReplacement, Rhs> >
//		{
//			assert(alpha == Scalar(1) && "scaling is not implemented");
//			EIGEN_ONLY_USED_FOR_DEBUG(alpha);
//
//			for (Index i = 0; i < lhs.cols(); ++i)
//				dst += rhs(i) * lhs.my_matrix().col(i);
//
//
//		};
//	}
//}