CC      = g++
CFLAGS  = -I. -Wall -Wextra
DEPS    = headers/header.hpp headers/Init.hpp
SOURCES = main.cpp $(DEPS:.hpp=.cpp)
OBJ = $(SOURCES:.cpp=.o)



%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Sim: $(OBJ)
	$(CC) -o ../$@ $(OBJ) $(CFLAGS)

clean:
	rm -f ∗~ $(OBJ) ∗.bak

depend:
	makedepend $(SOURCES)
# DO NOT DELETE

main.o: headers/Init.hpp headers/header.hpp headers/Eigen/Dense
main.o: headers/Eigen/Core
main.o: headers/Eigen/src/Core/util/DisableStupidWarnings.h
main.o: headers/Eigen/src/Core/util/Macros.h
main.o: headers/Eigen/src/Core/util/MKL_support.h
main.o: headers/Eigen/src/Core/util/Constants.h
main.o: headers/Eigen/src/Core/util/Meta.h
main.o: headers/Eigen/src/Core/util/ForwardDeclarations.h
main.o: headers/Eigen/src/Core/util/StaticAssert.h
main.o: headers/Eigen/src/Core/util/XprHelper.h
main.o: headers/Eigen/src/Core/util/Memory.h
main.o: headers/Eigen/src/Core/NumTraits.h
main.o: headers/Eigen/src/Core/MathFunctions.h
main.o: headers/Eigen/src/Core/GenericPacketMath.h
main.o: headers/Eigen/src/Core/MathFunctionsImpl.h
main.o: headers/Eigen/src/Core/arch/Default/ConjHelper.h
main.o: headers/Eigen/src/Core/arch/CUDA/Half.h
main.o: headers/Eigen/src/Core/arch/CUDA/PacketMathHalf.h
main.o: headers/Eigen/src/Core/arch/CUDA/TypeCasting.h
main.o: headers/Eigen/src/Core/arch/Default/Settings.h
main.o: headers/Eigen/src/Core/functors/TernaryFunctors.h
main.o: headers/Eigen/src/Core/functors/BinaryFunctors.h
main.o: headers/Eigen/src/Core/functors/UnaryFunctors.h
main.o: headers/Eigen/src/Core/functors/NullaryFunctors.h
main.o: headers/Eigen/src/Core/functors/StlFunctors.h
main.o: headers/Eigen/src/Core/functors/AssignmentFunctors.h
main.o: headers/Eigen/src/Core/arch/CUDA/Complex.h
main.o: headers/Eigen/src/Core/IO.h headers/Eigen/src/Core/DenseCoeffsBase.h
main.o: headers/Eigen/src/Core/DenseBase.h
main.o: headers/Eigen/src/plugins/BlockMethods.h
main.o: headers/Eigen/src/Core/MatrixBase.h
main.o: headers/Eigen/src/plugins/CommonCwiseUnaryOps.h
main.o: headers/Eigen/src/plugins/CommonCwiseBinaryOps.h
main.o: headers/Eigen/src/plugins/MatrixCwiseUnaryOps.h
main.o: headers/Eigen/src/plugins/MatrixCwiseBinaryOps.h
main.o: headers/Eigen/src/Core/EigenBase.h headers/Eigen/src/Core/Product.h
main.o: headers/Eigen/src/Core/CoreEvaluators.h
main.o: headers/Eigen/src/Core/AssignEvaluator.h
main.o: headers/Eigen/src/Core/Assign.h headers/Eigen/src/Core/ArrayBase.h
main.o: headers/Eigen/src/plugins/ArrayCwiseUnaryOps.h
main.o: headers/Eigen/src/plugins/ArrayCwiseBinaryOps.h
main.o: headers/Eigen/src/Core/util/BlasUtil.h
main.o: headers/Eigen/src/Core/DenseStorage.h
main.o: headers/Eigen/src/Core/NestByValue.h
main.o: headers/Eigen/src/Core/ReturnByValue.h
main.o: headers/Eigen/src/Core/NoAlias.h
main.o: headers/Eigen/src/Core/PlainObjectBase.h
main.o: headers/Eigen/src/Core/Matrix.h headers/Eigen/src/Core/Array.h
main.o: headers/Eigen/src/Core/CwiseTernaryOp.h
main.o: headers/Eigen/src/Core/CwiseBinaryOp.h
main.o: headers/Eigen/src/Core/CwiseUnaryOp.h
main.o: headers/Eigen/src/Core/CwiseNullaryOp.h
main.o: headers/Eigen/src/Core/CwiseUnaryView.h
main.o: headers/Eigen/src/Core/SelfCwiseBinaryOp.h
main.o: headers/Eigen/src/Core/Dot.h headers/Eigen/src/Core/StableNorm.h
main.o: headers/Eigen/src/Core/Stride.h headers/Eigen/src/Core/MapBase.h
main.o: headers/Eigen/src/Core/Map.h headers/Eigen/src/Core/Ref.h
main.o: headers/Eigen/src/Core/Block.h headers/Eigen/src/Core/VectorBlock.h
main.o: headers/Eigen/src/Core/Transpose.h
main.o: headers/Eigen/src/Core/DiagonalMatrix.h
main.o: headers/Eigen/src/Core/Diagonal.h
main.o: headers/Eigen/src/Core/DiagonalProduct.h
main.o: headers/Eigen/src/Core/Redux.h headers/Eigen/src/Core/Visitor.h
main.o: headers/Eigen/src/Core/Fuzzy.h headers/Eigen/src/Core/Swap.h
main.o: headers/Eigen/src/Core/CommaInitializer.h
main.o: headers/Eigen/src/Core/GeneralProduct.h
main.o: headers/Eigen/src/Core/Solve.h headers/Eigen/src/Core/Inverse.h
main.o: headers/Eigen/src/Core/SolverBase.h
main.o: headers/Eigen/src/Core/PermutationMatrix.h
main.o: headers/Eigen/src/Core/Transpositions.h
main.o: headers/Eigen/src/Core/TriangularMatrix.h
main.o: headers/Eigen/src/Core/SelfAdjointView.h
main.o: headers/Eigen/src/Core/products/GeneralBlockPanelKernel.h
main.o: headers/Eigen/src/Core/products/Parallelizer.h
main.o: headers/Eigen/src/Core/ProductEvaluators.h
main.o: headers/Eigen/src/Core/products/GeneralMatrixVector.h
main.o: headers/Eigen/src/Core/products/GeneralMatrixMatrix.h
main.o: headers/Eigen/src/Core/SolveTriangular.h
main.o: headers/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h
main.o: headers/Eigen/src/Core/products/SelfadjointMatrixVector.h
main.o: headers/Eigen/src/Core/products/SelfadjointMatrixMatrix.h
main.o: headers/Eigen/src/Core/products/SelfadjointProduct.h
main.o: headers/Eigen/src/Core/products/SelfadjointRank2Update.h
main.o: headers/Eigen/src/Core/products/TriangularMatrixVector.h
main.o: headers/Eigen/src/Core/products/TriangularMatrixMatrix.h
main.o: headers/Eigen/src/Core/products/TriangularSolverMatrix.h
main.o: headers/Eigen/src/Core/products/TriangularSolverVector.h
main.o: headers/Eigen/src/Core/BandMatrix.h
main.o: headers/Eigen/src/Core/CoreIterators.h
main.o: headers/Eigen/src/Core/ConditionEstimator.h
main.o: headers/Eigen/src/Core/BooleanRedux.h headers/Eigen/src/Core/Select.h
main.o: headers/Eigen/src/Core/VectorwiseOp.h headers/Eigen/src/Core/Random.h
main.o: headers/Eigen/src/Core/Replicate.h headers/Eigen/src/Core/Reverse.h
main.o: headers/Eigen/src/Core/ArrayWrapper.h
main.o: headers/Eigen/src/Core/GlobalFunctions.h
main.o: headers/Eigen/src/Core/util/ReenableStupidWarnings.h headers/Eigen/LU
main.o: headers/Eigen/src/misc/Kernel.h headers/Eigen/src/misc/Image.h
main.o: headers/Eigen/src/LU/FullPivLU.h headers/Eigen/src/LU/PartialPivLU.h
main.o: headers/Eigen/src/LU/Determinant.h headers/Eigen/src/LU/InverseImpl.h
main.o: headers/Eigen/Cholesky headers/Eigen/Jacobi
main.o: headers/Eigen/src/Jacobi/Jacobi.h headers/Eigen/src/Cholesky/LLT.h
main.o: headers/Eigen/src/Cholesky/LDLT.h headers/Eigen/QR
main.o: headers/Eigen/Householder headers/Eigen/src/Householder/Householder.h
main.o: headers/Eigen/src/Householder/HouseholderSequence.h
main.o: headers/Eigen/src/Householder/BlockHouseholder.h
main.o: headers/Eigen/src/QR/HouseholderQR.h
main.o: headers/Eigen/src/QR/FullPivHouseholderQR.h
main.o: headers/Eigen/src/QR/ColPivHouseholderQR.h
main.o: headers/Eigen/src/QR/CompleteOrthogonalDecomposition.h
main.o: headers/Eigen/SVD headers/Eigen/src/misc/RealSvd2x2.h
main.o: headers/Eigen/src/SVD/UpperBidiagonalization.h
main.o: headers/Eigen/src/SVD/SVDBase.h headers/Eigen/src/SVD/JacobiSVD.h
main.o: headers/Eigen/src/SVD/BDCSVD.h headers/Eigen/Geometry
main.o: headers/Eigen/src/Geometry/OrthoMethods.h
main.o: headers/Eigen/src/Geometry/EulerAngles.h
main.o: headers/Eigen/src/Geometry/Homogeneous.h
main.o: headers/Eigen/src/Geometry/RotationBase.h
main.o: headers/Eigen/src/Geometry/Rotation2D.h
main.o: headers/Eigen/src/Geometry/Quaternion.h
main.o: headers/Eigen/src/Geometry/AngleAxis.h
main.o: headers/Eigen/src/Geometry/Transform.h
main.o: headers/Eigen/src/Geometry/Translation.h
main.o: headers/Eigen/src/Geometry/Scaling.h
main.o: headers/Eigen/src/Geometry/Hyperplane.h
main.o: headers/Eigen/src/Geometry/ParametrizedLine.h
main.o: headers/Eigen/src/Geometry/AlignedBox.h
main.o: headers/Eigen/src/Geometry/Umeyama.h headers/Eigen/Eigenvalues
main.o: headers/Eigen/src/Eigenvalues/Tridiagonalization.h
main.o: headers/Eigen/src/Eigenvalues/RealSchur.h
main.o: headers/Eigen/src/Eigenvalues/./HessenbergDecomposition.h
main.o: headers/Eigen/src/Eigenvalues/EigenSolver.h
main.o: headers/Eigen/src/Eigenvalues/./RealSchur.h
main.o: headers/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h
main.o: headers/Eigen/src/Eigenvalues/./Tridiagonalization.h
main.o: headers/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h
main.o: headers/Eigen/src/Eigenvalues/HessenbergDecomposition.h
main.o: headers/Eigen/src/Eigenvalues/ComplexSchur.h
main.o: headers/Eigen/src/Eigenvalues/ComplexEigenSolver.h
main.o: headers/Eigen/src/Eigenvalues/./ComplexSchur.h
main.o: headers/Eigen/src/Eigenvalues/RealQZ.h
main.o: headers/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h
main.o: headers/Eigen/src/Eigenvalues/./RealQZ.h
main.o: headers/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h
headers/header.o: headers/header.hpp headers/Eigen/Dense headers/Eigen/Core
headers/header.o: headers/Eigen/src/Core/util/DisableStupidWarnings.h
headers/header.o: headers/Eigen/src/Core/util/Macros.h
headers/header.o: headers/Eigen/src/Core/util/MKL_support.h
headers/header.o: headers/Eigen/src/Core/util/Constants.h
headers/header.o: headers/Eigen/src/Core/util/Meta.h
headers/header.o: headers/Eigen/src/Core/util/ForwardDeclarations.h
headers/header.o: headers/Eigen/src/Core/util/StaticAssert.h
headers/header.o: headers/Eigen/src/Core/util/XprHelper.h
headers/header.o: headers/Eigen/src/Core/util/Memory.h
headers/header.o: headers/Eigen/src/Core/NumTraits.h
headers/header.o: headers/Eigen/src/Core/MathFunctions.h
headers/header.o: headers/Eigen/src/Core/GenericPacketMath.h
headers/header.o: headers/Eigen/src/Core/MathFunctionsImpl.h
headers/header.o: headers/Eigen/src/Core/arch/Default/ConjHelper.h
headers/header.o: headers/Eigen/src/Core/arch/CUDA/Half.h
headers/header.o: headers/Eigen/src/Core/arch/CUDA/PacketMathHalf.h
headers/header.o: headers/Eigen/src/Core/arch/CUDA/TypeCasting.h
headers/header.o: headers/Eigen/src/Core/arch/Default/Settings.h
headers/header.o: headers/Eigen/src/Core/functors/TernaryFunctors.h
headers/header.o: headers/Eigen/src/Core/functors/BinaryFunctors.h
headers/header.o: headers/Eigen/src/Core/functors/UnaryFunctors.h
headers/header.o: headers/Eigen/src/Core/functors/NullaryFunctors.h
headers/header.o: headers/Eigen/src/Core/functors/StlFunctors.h
headers/header.o: headers/Eigen/src/Core/functors/AssignmentFunctors.h
headers/header.o: headers/Eigen/src/Core/arch/CUDA/Complex.h
headers/header.o: headers/Eigen/src/Core/IO.h
headers/header.o: headers/Eigen/src/Core/DenseCoeffsBase.h
headers/header.o: headers/Eigen/src/Core/DenseBase.h
headers/header.o: headers/Eigen/src/plugins/BlockMethods.h
headers/header.o: headers/Eigen/src/Core/MatrixBase.h
headers/header.o: headers/Eigen/src/plugins/CommonCwiseUnaryOps.h
headers/header.o: headers/Eigen/src/plugins/CommonCwiseBinaryOps.h
headers/header.o: headers/Eigen/src/plugins/MatrixCwiseUnaryOps.h
headers/header.o: headers/Eigen/src/plugins/MatrixCwiseBinaryOps.h
headers/header.o: headers/Eigen/src/Core/EigenBase.h
headers/header.o: headers/Eigen/src/Core/Product.h
headers/header.o: headers/Eigen/src/Core/CoreEvaluators.h
headers/header.o: headers/Eigen/src/Core/AssignEvaluator.h
headers/header.o: headers/Eigen/src/Core/Assign.h
headers/header.o: headers/Eigen/src/Core/ArrayBase.h
headers/header.o: headers/Eigen/src/plugins/ArrayCwiseUnaryOps.h
headers/header.o: headers/Eigen/src/plugins/ArrayCwiseBinaryOps.h
headers/header.o: headers/Eigen/src/Core/util/BlasUtil.h
headers/header.o: headers/Eigen/src/Core/DenseStorage.h
headers/header.o: headers/Eigen/src/Core/NestByValue.h
headers/header.o: headers/Eigen/src/Core/ReturnByValue.h
headers/header.o: headers/Eigen/src/Core/NoAlias.h
headers/header.o: headers/Eigen/src/Core/PlainObjectBase.h
headers/header.o: headers/Eigen/src/Core/Matrix.h
headers/header.o: headers/Eigen/src/Core/Array.h
headers/header.o: headers/Eigen/src/Core/CwiseTernaryOp.h
headers/header.o: headers/Eigen/src/Core/CwiseBinaryOp.h
headers/header.o: headers/Eigen/src/Core/CwiseUnaryOp.h
headers/header.o: headers/Eigen/src/Core/CwiseNullaryOp.h
headers/header.o: headers/Eigen/src/Core/CwiseUnaryView.h
headers/header.o: headers/Eigen/src/Core/SelfCwiseBinaryOp.h
headers/header.o: headers/Eigen/src/Core/Dot.h
headers/header.o: headers/Eigen/src/Core/StableNorm.h
headers/header.o: headers/Eigen/src/Core/Stride.h
headers/header.o: headers/Eigen/src/Core/MapBase.h
headers/header.o: headers/Eigen/src/Core/Map.h headers/Eigen/src/Core/Ref.h
headers/header.o: headers/Eigen/src/Core/Block.h
headers/header.o: headers/Eigen/src/Core/VectorBlock.h
headers/header.o: headers/Eigen/src/Core/Transpose.h
headers/header.o: headers/Eigen/src/Core/DiagonalMatrix.h
headers/header.o: headers/Eigen/src/Core/Diagonal.h
headers/header.o: headers/Eigen/src/Core/DiagonalProduct.h
headers/header.o: headers/Eigen/src/Core/Redux.h
headers/header.o: headers/Eigen/src/Core/Visitor.h
headers/header.o: headers/Eigen/src/Core/Fuzzy.h
headers/header.o: headers/Eigen/src/Core/Swap.h
headers/header.o: headers/Eigen/src/Core/CommaInitializer.h
headers/header.o: headers/Eigen/src/Core/GeneralProduct.h
headers/header.o: headers/Eigen/src/Core/Solve.h
headers/header.o: headers/Eigen/src/Core/Inverse.h
headers/header.o: headers/Eigen/src/Core/SolverBase.h
headers/header.o: headers/Eigen/src/Core/PermutationMatrix.h
headers/header.o: headers/Eigen/src/Core/Transpositions.h
headers/header.o: headers/Eigen/src/Core/TriangularMatrix.h
headers/header.o: headers/Eigen/src/Core/SelfAdjointView.h
headers/header.o: headers/Eigen/src/Core/products/GeneralBlockPanelKernel.h
headers/header.o: headers/Eigen/src/Core/products/Parallelizer.h
headers/header.o: headers/Eigen/src/Core/ProductEvaluators.h
headers/header.o: headers/Eigen/src/Core/products/GeneralMatrixVector.h
headers/header.o: headers/Eigen/src/Core/products/GeneralMatrixMatrix.h
headers/header.o: headers/Eigen/src/Core/SolveTriangular.h
headers/header.o: headers/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h
headers/header.o: headers/Eigen/src/Core/products/SelfadjointMatrixVector.h
headers/header.o: headers/Eigen/src/Core/products/SelfadjointMatrixMatrix.h
headers/header.o: headers/Eigen/src/Core/products/SelfadjointProduct.h
headers/header.o: headers/Eigen/src/Core/products/SelfadjointRank2Update.h
headers/header.o: headers/Eigen/src/Core/products/TriangularMatrixVector.h
headers/header.o: headers/Eigen/src/Core/products/TriangularMatrixMatrix.h
headers/header.o: headers/Eigen/src/Core/products/TriangularSolverMatrix.h
headers/header.o: headers/Eigen/src/Core/products/TriangularSolverVector.h
headers/header.o: headers/Eigen/src/Core/BandMatrix.h
headers/header.o: headers/Eigen/src/Core/CoreIterators.h
headers/header.o: headers/Eigen/src/Core/ConditionEstimator.h
headers/header.o: headers/Eigen/src/Core/BooleanRedux.h
headers/header.o: headers/Eigen/src/Core/Select.h
headers/header.o: headers/Eigen/src/Core/VectorwiseOp.h
headers/header.o: headers/Eigen/src/Core/Random.h
headers/header.o: headers/Eigen/src/Core/Replicate.h
headers/header.o: headers/Eigen/src/Core/Reverse.h
headers/header.o: headers/Eigen/src/Core/ArrayWrapper.h
headers/header.o: headers/Eigen/src/Core/GlobalFunctions.h
headers/header.o: headers/Eigen/src/Core/util/ReenableStupidWarnings.h
headers/header.o: headers/Eigen/LU headers/Eigen/src/misc/Kernel.h
headers/header.o: headers/Eigen/src/misc/Image.h
headers/header.o: headers/Eigen/src/LU/FullPivLU.h
headers/header.o: headers/Eigen/src/LU/PartialPivLU.h
headers/header.o: headers/Eigen/src/LU/Determinant.h
headers/header.o: headers/Eigen/src/LU/InverseImpl.h headers/Eigen/Cholesky
headers/header.o: headers/Eigen/Jacobi headers/Eigen/src/Jacobi/Jacobi.h
headers/header.o: headers/Eigen/src/Cholesky/LLT.h
headers/header.o: headers/Eigen/src/Cholesky/LDLT.h headers/Eigen/QR
headers/header.o: headers/Eigen/Householder
headers/header.o: headers/Eigen/src/Householder/Householder.h
headers/header.o: headers/Eigen/src/Householder/HouseholderSequence.h
headers/header.o: headers/Eigen/src/Householder/BlockHouseholder.h
headers/header.o: headers/Eigen/src/QR/HouseholderQR.h
headers/header.o: headers/Eigen/src/QR/FullPivHouseholderQR.h
headers/header.o: headers/Eigen/src/QR/ColPivHouseholderQR.h
headers/header.o: headers/Eigen/src/QR/CompleteOrthogonalDecomposition.h
headers/header.o: headers/Eigen/SVD headers/Eigen/src/misc/RealSvd2x2.h
headers/header.o: headers/Eigen/src/SVD/UpperBidiagonalization.h
headers/header.o: headers/Eigen/src/SVD/SVDBase.h
headers/header.o: headers/Eigen/src/SVD/JacobiSVD.h
headers/header.o: headers/Eigen/src/SVD/BDCSVD.h headers/Eigen/Geometry
headers/header.o: headers/Eigen/src/Geometry/OrthoMethods.h
headers/header.o: headers/Eigen/src/Geometry/EulerAngles.h
headers/header.o: headers/Eigen/src/Geometry/Homogeneous.h
headers/header.o: headers/Eigen/src/Geometry/RotationBase.h
headers/header.o: headers/Eigen/src/Geometry/Rotation2D.h
headers/header.o: headers/Eigen/src/Geometry/Quaternion.h
headers/header.o: headers/Eigen/src/Geometry/AngleAxis.h
headers/header.o: headers/Eigen/src/Geometry/Transform.h
headers/header.o: headers/Eigen/src/Geometry/Translation.h
headers/header.o: headers/Eigen/src/Geometry/Scaling.h
headers/header.o: headers/Eigen/src/Geometry/Hyperplane.h
headers/header.o: headers/Eigen/src/Geometry/ParametrizedLine.h
headers/header.o: headers/Eigen/src/Geometry/AlignedBox.h
headers/header.o: headers/Eigen/src/Geometry/Umeyama.h
headers/header.o: headers/Eigen/Eigenvalues
headers/header.o: headers/Eigen/src/Eigenvalues/Tridiagonalization.h
headers/header.o: headers/Eigen/src/Eigenvalues/RealSchur.h
headers/header.o: headers/Eigen/src/Eigenvalues/./HessenbergDecomposition.h
headers/header.o: headers/Eigen/src/Eigenvalues/EigenSolver.h
headers/header.o: headers/Eigen/src/Eigenvalues/./RealSchur.h
headers/header.o: headers/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h
headers/header.o: headers/Eigen/src/Eigenvalues/./Tridiagonalization.h
headers/header.o: headers/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h
headers/header.o: headers/Eigen/src/Eigenvalues/HessenbergDecomposition.h
headers/header.o: headers/Eigen/src/Eigenvalues/ComplexSchur.h
headers/header.o: headers/Eigen/src/Eigenvalues/ComplexEigenSolver.h
headers/header.o: headers/Eigen/src/Eigenvalues/./ComplexSchur.h
headers/header.o: headers/Eigen/src/Eigenvalues/RealQZ.h
headers/header.o: headers/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h
headers/header.o: headers/Eigen/src/Eigenvalues/./RealQZ.h
headers/header.o: headers/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h
headers/Init.o: headers/Init.hpp headers/header.hpp headers/Eigen/Dense
headers/Init.o: headers/Eigen/Core
headers/Init.o: headers/Eigen/src/Core/util/DisableStupidWarnings.h
headers/Init.o: headers/Eigen/src/Core/util/Macros.h
headers/Init.o: headers/Eigen/src/Core/util/MKL_support.h
headers/Init.o: headers/Eigen/src/Core/util/Constants.h
headers/Init.o: headers/Eigen/src/Core/util/Meta.h
headers/Init.o: headers/Eigen/src/Core/util/ForwardDeclarations.h
headers/Init.o: headers/Eigen/src/Core/util/StaticAssert.h
headers/Init.o: headers/Eigen/src/Core/util/XprHelper.h
headers/Init.o: headers/Eigen/src/Core/util/Memory.h
headers/Init.o: headers/Eigen/src/Core/NumTraits.h
headers/Init.o: headers/Eigen/src/Core/MathFunctions.h
headers/Init.o: headers/Eigen/src/Core/GenericPacketMath.h
headers/Init.o: headers/Eigen/src/Core/MathFunctionsImpl.h
headers/Init.o: headers/Eigen/src/Core/arch/Default/ConjHelper.h
headers/Init.o: headers/Eigen/src/Core/arch/CUDA/Half.h
headers/Init.o: headers/Eigen/src/Core/arch/CUDA/PacketMathHalf.h
headers/Init.o: headers/Eigen/src/Core/arch/CUDA/TypeCasting.h
headers/Init.o: headers/Eigen/src/Core/arch/Default/Settings.h
headers/Init.o: headers/Eigen/src/Core/functors/TernaryFunctors.h
headers/Init.o: headers/Eigen/src/Core/functors/BinaryFunctors.h
headers/Init.o: headers/Eigen/src/Core/functors/UnaryFunctors.h
headers/Init.o: headers/Eigen/src/Core/functors/NullaryFunctors.h
headers/Init.o: headers/Eigen/src/Core/functors/StlFunctors.h
headers/Init.o: headers/Eigen/src/Core/functors/AssignmentFunctors.h
headers/Init.o: headers/Eigen/src/Core/arch/CUDA/Complex.h
headers/Init.o: headers/Eigen/src/Core/IO.h
headers/Init.o: headers/Eigen/src/Core/DenseCoeffsBase.h
headers/Init.o: headers/Eigen/src/Core/DenseBase.h
headers/Init.o: headers/Eigen/src/plugins/BlockMethods.h
headers/Init.o: headers/Eigen/src/Core/MatrixBase.h
headers/Init.o: headers/Eigen/src/plugins/CommonCwiseUnaryOps.h
headers/Init.o: headers/Eigen/src/plugins/CommonCwiseBinaryOps.h
headers/Init.o: headers/Eigen/src/plugins/MatrixCwiseUnaryOps.h
headers/Init.o: headers/Eigen/src/plugins/MatrixCwiseBinaryOps.h
headers/Init.o: headers/Eigen/src/Core/EigenBase.h
headers/Init.o: headers/Eigen/src/Core/Product.h
headers/Init.o: headers/Eigen/src/Core/CoreEvaluators.h
headers/Init.o: headers/Eigen/src/Core/AssignEvaluator.h
headers/Init.o: headers/Eigen/src/Core/Assign.h
headers/Init.o: headers/Eigen/src/Core/ArrayBase.h
headers/Init.o: headers/Eigen/src/plugins/ArrayCwiseUnaryOps.h
headers/Init.o: headers/Eigen/src/plugins/ArrayCwiseBinaryOps.h
headers/Init.o: headers/Eigen/src/Core/util/BlasUtil.h
headers/Init.o: headers/Eigen/src/Core/DenseStorage.h
headers/Init.o: headers/Eigen/src/Core/NestByValue.h
headers/Init.o: headers/Eigen/src/Core/ReturnByValue.h
headers/Init.o: headers/Eigen/src/Core/NoAlias.h
headers/Init.o: headers/Eigen/src/Core/PlainObjectBase.h
headers/Init.o: headers/Eigen/src/Core/Matrix.h
headers/Init.o: headers/Eigen/src/Core/Array.h
headers/Init.o: headers/Eigen/src/Core/CwiseTernaryOp.h
headers/Init.o: headers/Eigen/src/Core/CwiseBinaryOp.h
headers/Init.o: headers/Eigen/src/Core/CwiseUnaryOp.h
headers/Init.o: headers/Eigen/src/Core/CwiseNullaryOp.h
headers/Init.o: headers/Eigen/src/Core/CwiseUnaryView.h
headers/Init.o: headers/Eigen/src/Core/SelfCwiseBinaryOp.h
headers/Init.o: headers/Eigen/src/Core/Dot.h
headers/Init.o: headers/Eigen/src/Core/StableNorm.h
headers/Init.o: headers/Eigen/src/Core/Stride.h
headers/Init.o: headers/Eigen/src/Core/MapBase.h headers/Eigen/src/Core/Map.h
headers/Init.o: headers/Eigen/src/Core/Ref.h headers/Eigen/src/Core/Block.h
headers/Init.o: headers/Eigen/src/Core/VectorBlock.h
headers/Init.o: headers/Eigen/src/Core/Transpose.h
headers/Init.o: headers/Eigen/src/Core/DiagonalMatrix.h
headers/Init.o: headers/Eigen/src/Core/Diagonal.h
headers/Init.o: headers/Eigen/src/Core/DiagonalProduct.h
headers/Init.o: headers/Eigen/src/Core/Redux.h
headers/Init.o: headers/Eigen/src/Core/Visitor.h
headers/Init.o: headers/Eigen/src/Core/Fuzzy.h headers/Eigen/src/Core/Swap.h
headers/Init.o: headers/Eigen/src/Core/CommaInitializer.h
headers/Init.o: headers/Eigen/src/Core/GeneralProduct.h
headers/Init.o: headers/Eigen/src/Core/Solve.h
headers/Init.o: headers/Eigen/src/Core/Inverse.h
headers/Init.o: headers/Eigen/src/Core/SolverBase.h
headers/Init.o: headers/Eigen/src/Core/PermutationMatrix.h
headers/Init.o: headers/Eigen/src/Core/Transpositions.h
headers/Init.o: headers/Eigen/src/Core/TriangularMatrix.h
headers/Init.o: headers/Eigen/src/Core/SelfAdjointView.h
headers/Init.o: headers/Eigen/src/Core/products/GeneralBlockPanelKernel.h
headers/Init.o: headers/Eigen/src/Core/products/Parallelizer.h
headers/Init.o: headers/Eigen/src/Core/ProductEvaluators.h
headers/Init.o: headers/Eigen/src/Core/products/GeneralMatrixVector.h
headers/Init.o: headers/Eigen/src/Core/products/GeneralMatrixMatrix.h
headers/Init.o: headers/Eigen/src/Core/SolveTriangular.h
headers/Init.o: headers/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h
headers/Init.o: headers/Eigen/src/Core/products/SelfadjointMatrixVector.h
headers/Init.o: headers/Eigen/src/Core/products/SelfadjointMatrixMatrix.h
headers/Init.o: headers/Eigen/src/Core/products/SelfadjointProduct.h
headers/Init.o: headers/Eigen/src/Core/products/SelfadjointRank2Update.h
headers/Init.o: headers/Eigen/src/Core/products/TriangularMatrixVector.h
headers/Init.o: headers/Eigen/src/Core/products/TriangularMatrixMatrix.h
headers/Init.o: headers/Eigen/src/Core/products/TriangularSolverMatrix.h
headers/Init.o: headers/Eigen/src/Core/products/TriangularSolverVector.h
headers/Init.o: headers/Eigen/src/Core/BandMatrix.h
headers/Init.o: headers/Eigen/src/Core/CoreIterators.h
headers/Init.o: headers/Eigen/src/Core/ConditionEstimator.h
headers/Init.o: headers/Eigen/src/Core/BooleanRedux.h
headers/Init.o: headers/Eigen/src/Core/Select.h
headers/Init.o: headers/Eigen/src/Core/VectorwiseOp.h
headers/Init.o: headers/Eigen/src/Core/Random.h
headers/Init.o: headers/Eigen/src/Core/Replicate.h
headers/Init.o: headers/Eigen/src/Core/Reverse.h
headers/Init.o: headers/Eigen/src/Core/ArrayWrapper.h
headers/Init.o: headers/Eigen/src/Core/GlobalFunctions.h
headers/Init.o: headers/Eigen/src/Core/util/ReenableStupidWarnings.h
headers/Init.o: headers/Eigen/LU headers/Eigen/src/misc/Kernel.h
headers/Init.o: headers/Eigen/src/misc/Image.h
headers/Init.o: headers/Eigen/src/LU/FullPivLU.h
headers/Init.o: headers/Eigen/src/LU/PartialPivLU.h
headers/Init.o: headers/Eigen/src/LU/Determinant.h
headers/Init.o: headers/Eigen/src/LU/InverseImpl.h headers/Eigen/Cholesky
headers/Init.o: headers/Eigen/Jacobi headers/Eigen/src/Jacobi/Jacobi.h
headers/Init.o: headers/Eigen/src/Cholesky/LLT.h
headers/Init.o: headers/Eigen/src/Cholesky/LDLT.h headers/Eigen/QR
headers/Init.o: headers/Eigen/Householder
headers/Init.o: headers/Eigen/src/Householder/Householder.h
headers/Init.o: headers/Eigen/src/Householder/HouseholderSequence.h
headers/Init.o: headers/Eigen/src/Householder/BlockHouseholder.h
headers/Init.o: headers/Eigen/src/QR/HouseholderQR.h
headers/Init.o: headers/Eigen/src/QR/FullPivHouseholderQR.h
headers/Init.o: headers/Eigen/src/QR/ColPivHouseholderQR.h
headers/Init.o: headers/Eigen/src/QR/CompleteOrthogonalDecomposition.h
headers/Init.o: headers/Eigen/SVD headers/Eigen/src/misc/RealSvd2x2.h
headers/Init.o: headers/Eigen/src/SVD/UpperBidiagonalization.h
headers/Init.o: headers/Eigen/src/SVD/SVDBase.h
headers/Init.o: headers/Eigen/src/SVD/JacobiSVD.h
headers/Init.o: headers/Eigen/src/SVD/BDCSVD.h headers/Eigen/Geometry
headers/Init.o: headers/Eigen/src/Geometry/OrthoMethods.h
headers/Init.o: headers/Eigen/src/Geometry/EulerAngles.h
headers/Init.o: headers/Eigen/src/Geometry/Homogeneous.h
headers/Init.o: headers/Eigen/src/Geometry/RotationBase.h
headers/Init.o: headers/Eigen/src/Geometry/Rotation2D.h
headers/Init.o: headers/Eigen/src/Geometry/Quaternion.h
headers/Init.o: headers/Eigen/src/Geometry/AngleAxis.h
headers/Init.o: headers/Eigen/src/Geometry/Transform.h
headers/Init.o: headers/Eigen/src/Geometry/Translation.h
headers/Init.o: headers/Eigen/src/Geometry/Scaling.h
headers/Init.o: headers/Eigen/src/Geometry/Hyperplane.h
headers/Init.o: headers/Eigen/src/Geometry/ParametrizedLine.h
headers/Init.o: headers/Eigen/src/Geometry/AlignedBox.h
headers/Init.o: headers/Eigen/src/Geometry/Umeyama.h
headers/Init.o: headers/Eigen/Eigenvalues
headers/Init.o: headers/Eigen/src/Eigenvalues/Tridiagonalization.h
headers/Init.o: headers/Eigen/src/Eigenvalues/RealSchur.h
headers/Init.o: headers/Eigen/src/Eigenvalues/./HessenbergDecomposition.h
headers/Init.o: headers/Eigen/src/Eigenvalues/EigenSolver.h
headers/Init.o: headers/Eigen/src/Eigenvalues/./RealSchur.h
headers/Init.o: headers/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h
headers/Init.o: headers/Eigen/src/Eigenvalues/./Tridiagonalization.h
headers/Init.o: headers/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h
headers/Init.o: headers/Eigen/src/Eigenvalues/HessenbergDecomposition.h
headers/Init.o: headers/Eigen/src/Eigenvalues/ComplexSchur.h
headers/Init.o: headers/Eigen/src/Eigenvalues/ComplexEigenSolver.h
headers/Init.o: headers/Eigen/src/Eigenvalues/./ComplexSchur.h
headers/Init.o: headers/Eigen/src/Eigenvalues/RealQZ.h
headers/Init.o: headers/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h
headers/Init.o: headers/Eigen/src/Eigenvalues/./RealQZ.h
headers/Init.o: headers/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h
