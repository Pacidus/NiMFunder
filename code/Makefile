CC      = g++
CFLAGS  = -Ofast -W -Wall -Wextra -DNDEBUG

DEPS    = headers/Class.hpp headers/Methods.hpp
SSOURCES = main.cpp $(DEPS:.hpp=.cpp)
SOBJ = $(SSOURCES:.cpp=.o)

MSOURCES = Map.cpp $(DEPS:.hpp=.cpp)
MOBJ = $(MSOURCES:.cpp=.o)



%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Sim: $(SOBJ)
	$(CC) -o ../$@ $(SOBJ) $(CFLAGS)

Map: $(MOBJ)
	$(CC) -o ../$@ $(MOBJ) $(CFLAGS)

clean:
	rm -f ∗~ $(MOBJ) $(SOBJ) ∗.bak

Sdepend:
	makedepend $(SSOURCES)

Mdepend:
	makedepend $(MSOURCES)

# DO NOT DELETE

main.o: headers/Methods.hpp headers/Class.hpp headers/Eigen/Dense
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
headers/Class.o: headers/Class.hpp headers/Eigen/Dense headers/Eigen/Core
headers/Class.o: headers/Eigen/src/Core/util/DisableStupidWarnings.h
headers/Class.o: headers/Eigen/src/Core/util/Macros.h
headers/Class.o: headers/Eigen/src/Core/util/MKL_support.h
headers/Class.o: headers/Eigen/src/Core/util/Constants.h
headers/Class.o: headers/Eigen/src/Core/util/Meta.h
headers/Class.o: headers/Eigen/src/Core/util/ForwardDeclarations.h
headers/Class.o: headers/Eigen/src/Core/util/StaticAssert.h
headers/Class.o: headers/Eigen/src/Core/util/XprHelper.h
headers/Class.o: headers/Eigen/src/Core/util/Memory.h
headers/Class.o: headers/Eigen/src/Core/NumTraits.h
headers/Class.o: headers/Eigen/src/Core/MathFunctions.h
headers/Class.o: headers/Eigen/src/Core/GenericPacketMath.h
headers/Class.o: headers/Eigen/src/Core/MathFunctionsImpl.h
headers/Class.o: headers/Eigen/src/Core/arch/Default/ConjHelper.h
headers/Class.o: headers/Eigen/src/Core/arch/CUDA/Half.h
headers/Class.o: headers/Eigen/src/Core/arch/CUDA/PacketMathHalf.h
headers/Class.o: headers/Eigen/src/Core/arch/CUDA/TypeCasting.h
headers/Class.o: headers/Eigen/src/Core/arch/Default/Settings.h
headers/Class.o: headers/Eigen/src/Core/functors/TernaryFunctors.h
headers/Class.o: headers/Eigen/src/Core/functors/BinaryFunctors.h
headers/Class.o: headers/Eigen/src/Core/functors/UnaryFunctors.h
headers/Class.o: headers/Eigen/src/Core/functors/NullaryFunctors.h
headers/Class.o: headers/Eigen/src/Core/functors/StlFunctors.h
headers/Class.o: headers/Eigen/src/Core/functors/AssignmentFunctors.h
headers/Class.o: headers/Eigen/src/Core/arch/CUDA/Complex.h
headers/Class.o: headers/Eigen/src/Core/IO.h
headers/Class.o: headers/Eigen/src/Core/DenseCoeffsBase.h
headers/Class.o: headers/Eigen/src/Core/DenseBase.h
headers/Class.o: headers/Eigen/src/plugins/BlockMethods.h
headers/Class.o: headers/Eigen/src/Core/MatrixBase.h
headers/Class.o: headers/Eigen/src/plugins/CommonCwiseUnaryOps.h
headers/Class.o: headers/Eigen/src/plugins/CommonCwiseBinaryOps.h
headers/Class.o: headers/Eigen/src/plugins/MatrixCwiseUnaryOps.h
headers/Class.o: headers/Eigen/src/plugins/MatrixCwiseBinaryOps.h
headers/Class.o: headers/Eigen/src/Core/EigenBase.h
headers/Class.o: headers/Eigen/src/Core/Product.h
headers/Class.o: headers/Eigen/src/Core/CoreEvaluators.h
headers/Class.o: headers/Eigen/src/Core/AssignEvaluator.h
headers/Class.o: headers/Eigen/src/Core/Assign.h
headers/Class.o: headers/Eigen/src/Core/ArrayBase.h
headers/Class.o: headers/Eigen/src/plugins/ArrayCwiseUnaryOps.h
headers/Class.o: headers/Eigen/src/plugins/ArrayCwiseBinaryOps.h
headers/Class.o: headers/Eigen/src/Core/util/BlasUtil.h
headers/Class.o: headers/Eigen/src/Core/DenseStorage.h
headers/Class.o: headers/Eigen/src/Core/NestByValue.h
headers/Class.o: headers/Eigen/src/Core/ReturnByValue.h
headers/Class.o: headers/Eigen/src/Core/NoAlias.h
headers/Class.o: headers/Eigen/src/Core/PlainObjectBase.h
headers/Class.o: headers/Eigen/src/Core/Matrix.h
headers/Class.o: headers/Eigen/src/Core/Array.h
headers/Class.o: headers/Eigen/src/Core/CwiseTernaryOp.h
headers/Class.o: headers/Eigen/src/Core/CwiseBinaryOp.h
headers/Class.o: headers/Eigen/src/Core/CwiseUnaryOp.h
headers/Class.o: headers/Eigen/src/Core/CwiseNullaryOp.h
headers/Class.o: headers/Eigen/src/Core/CwiseUnaryView.h
headers/Class.o: headers/Eigen/src/Core/SelfCwiseBinaryOp.h
headers/Class.o: headers/Eigen/src/Core/Dot.h
headers/Class.o: headers/Eigen/src/Core/StableNorm.h
headers/Class.o: headers/Eigen/src/Core/Stride.h
headers/Class.o: headers/Eigen/src/Core/MapBase.h
headers/Class.o: headers/Eigen/src/Core/Map.h headers/Eigen/src/Core/Ref.h
headers/Class.o: headers/Eigen/src/Core/Block.h
headers/Class.o: headers/Eigen/src/Core/VectorBlock.h
headers/Class.o: headers/Eigen/src/Core/Transpose.h
headers/Class.o: headers/Eigen/src/Core/DiagonalMatrix.h
headers/Class.o: headers/Eigen/src/Core/Diagonal.h
headers/Class.o: headers/Eigen/src/Core/DiagonalProduct.h
headers/Class.o: headers/Eigen/src/Core/Redux.h
headers/Class.o: headers/Eigen/src/Core/Visitor.h
headers/Class.o: headers/Eigen/src/Core/Fuzzy.h headers/Eigen/src/Core/Swap.h
headers/Class.o: headers/Eigen/src/Core/CommaInitializer.h
headers/Class.o: headers/Eigen/src/Core/GeneralProduct.h
headers/Class.o: headers/Eigen/src/Core/Solve.h
headers/Class.o: headers/Eigen/src/Core/Inverse.h
headers/Class.o: headers/Eigen/src/Core/SolverBase.h
headers/Class.o: headers/Eigen/src/Core/PermutationMatrix.h
headers/Class.o: headers/Eigen/src/Core/Transpositions.h
headers/Class.o: headers/Eigen/src/Core/TriangularMatrix.h
headers/Class.o: headers/Eigen/src/Core/SelfAdjointView.h
headers/Class.o: headers/Eigen/src/Core/products/GeneralBlockPanelKernel.h
headers/Class.o: headers/Eigen/src/Core/products/Parallelizer.h
headers/Class.o: headers/Eigen/src/Core/ProductEvaluators.h
headers/Class.o: headers/Eigen/src/Core/products/GeneralMatrixVector.h
headers/Class.o: headers/Eigen/src/Core/products/GeneralMatrixMatrix.h
headers/Class.o: headers/Eigen/src/Core/SolveTriangular.h
headers/Class.o: headers/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h
headers/Class.o: headers/Eigen/src/Core/products/SelfadjointMatrixVector.h
headers/Class.o: headers/Eigen/src/Core/products/SelfadjointMatrixMatrix.h
headers/Class.o: headers/Eigen/src/Core/products/SelfadjointProduct.h
headers/Class.o: headers/Eigen/src/Core/products/SelfadjointRank2Update.h
headers/Class.o: headers/Eigen/src/Core/products/TriangularMatrixVector.h
headers/Class.o: headers/Eigen/src/Core/products/TriangularMatrixMatrix.h
headers/Class.o: headers/Eigen/src/Core/products/TriangularSolverMatrix.h
headers/Class.o: headers/Eigen/src/Core/products/TriangularSolverVector.h
headers/Class.o: headers/Eigen/src/Core/BandMatrix.h
headers/Class.o: headers/Eigen/src/Core/CoreIterators.h
headers/Class.o: headers/Eigen/src/Core/ConditionEstimator.h
headers/Class.o: headers/Eigen/src/Core/BooleanRedux.h
headers/Class.o: headers/Eigen/src/Core/Select.h
headers/Class.o: headers/Eigen/src/Core/VectorwiseOp.h
headers/Class.o: headers/Eigen/src/Core/Random.h
headers/Class.o: headers/Eigen/src/Core/Replicate.h
headers/Class.o: headers/Eigen/src/Core/Reverse.h
headers/Class.o: headers/Eigen/src/Core/ArrayWrapper.h
headers/Class.o: headers/Eigen/src/Core/GlobalFunctions.h
headers/Class.o: headers/Eigen/src/Core/util/ReenableStupidWarnings.h
headers/Class.o: headers/Eigen/LU headers/Eigen/src/misc/Kernel.h
headers/Class.o: headers/Eigen/src/misc/Image.h
headers/Class.o: headers/Eigen/src/LU/FullPivLU.h
headers/Class.o: headers/Eigen/src/LU/PartialPivLU.h
headers/Class.o: headers/Eigen/src/LU/Determinant.h
headers/Class.o: headers/Eigen/src/LU/InverseImpl.h headers/Eigen/Cholesky
headers/Class.o: headers/Eigen/Jacobi headers/Eigen/src/Jacobi/Jacobi.h
headers/Class.o: headers/Eigen/src/Cholesky/LLT.h
headers/Class.o: headers/Eigen/src/Cholesky/LDLT.h headers/Eigen/QR
headers/Class.o: headers/Eigen/Householder
headers/Class.o: headers/Eigen/src/Householder/Householder.h
headers/Class.o: headers/Eigen/src/Householder/HouseholderSequence.h
headers/Class.o: headers/Eigen/src/Householder/BlockHouseholder.h
headers/Class.o: headers/Eigen/src/QR/HouseholderQR.h
headers/Class.o: headers/Eigen/src/QR/FullPivHouseholderQR.h
headers/Class.o: headers/Eigen/src/QR/ColPivHouseholderQR.h
headers/Class.o: headers/Eigen/src/QR/CompleteOrthogonalDecomposition.h
headers/Class.o: headers/Eigen/SVD headers/Eigen/src/misc/RealSvd2x2.h
headers/Class.o: headers/Eigen/src/SVD/UpperBidiagonalization.h
headers/Class.o: headers/Eigen/src/SVD/SVDBase.h
headers/Class.o: headers/Eigen/src/SVD/JacobiSVD.h
headers/Class.o: headers/Eigen/src/SVD/BDCSVD.h headers/Eigen/Geometry
headers/Class.o: headers/Eigen/src/Geometry/OrthoMethods.h
headers/Class.o: headers/Eigen/src/Geometry/EulerAngles.h
headers/Class.o: headers/Eigen/src/Geometry/Homogeneous.h
headers/Class.o: headers/Eigen/src/Geometry/RotationBase.h
headers/Class.o: headers/Eigen/src/Geometry/Rotation2D.h
headers/Class.o: headers/Eigen/src/Geometry/Quaternion.h
headers/Class.o: headers/Eigen/src/Geometry/AngleAxis.h
headers/Class.o: headers/Eigen/src/Geometry/Transform.h
headers/Class.o: headers/Eigen/src/Geometry/Translation.h
headers/Class.o: headers/Eigen/src/Geometry/Scaling.h
headers/Class.o: headers/Eigen/src/Geometry/Hyperplane.h
headers/Class.o: headers/Eigen/src/Geometry/ParametrizedLine.h
headers/Class.o: headers/Eigen/src/Geometry/AlignedBox.h
headers/Class.o: headers/Eigen/src/Geometry/Umeyama.h
headers/Class.o: headers/Eigen/Eigenvalues
headers/Class.o: headers/Eigen/src/Eigenvalues/Tridiagonalization.h
headers/Class.o: headers/Eigen/src/Eigenvalues/RealSchur.h
headers/Class.o: headers/Eigen/src/Eigenvalues/./HessenbergDecomposition.h
headers/Class.o: headers/Eigen/src/Eigenvalues/EigenSolver.h
headers/Class.o: headers/Eigen/src/Eigenvalues/./RealSchur.h
headers/Class.o: headers/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h
headers/Class.o: headers/Eigen/src/Eigenvalues/./Tridiagonalization.h
headers/Class.o: headers/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h
headers/Class.o: headers/Eigen/src/Eigenvalues/HessenbergDecomposition.h
headers/Class.o: headers/Eigen/src/Eigenvalues/ComplexSchur.h
headers/Class.o: headers/Eigen/src/Eigenvalues/ComplexEigenSolver.h
headers/Class.o: headers/Eigen/src/Eigenvalues/./ComplexSchur.h
headers/Class.o: headers/Eigen/src/Eigenvalues/RealQZ.h
headers/Class.o: headers/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h
headers/Class.o: headers/Eigen/src/Eigenvalues/./RealQZ.h
headers/Class.o: headers/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h
headers/Methods.o: headers/Methods.hpp headers/Class.hpp headers/Eigen/Dense
headers/Methods.o: headers/Eigen/Core
headers/Methods.o: headers/Eigen/src/Core/util/DisableStupidWarnings.h
headers/Methods.o: headers/Eigen/src/Core/util/Macros.h
headers/Methods.o: headers/Eigen/src/Core/util/MKL_support.h
headers/Methods.o: headers/Eigen/src/Core/util/Constants.h
headers/Methods.o: headers/Eigen/src/Core/util/Meta.h
headers/Methods.o: headers/Eigen/src/Core/util/ForwardDeclarations.h
headers/Methods.o: headers/Eigen/src/Core/util/StaticAssert.h
headers/Methods.o: headers/Eigen/src/Core/util/XprHelper.h
headers/Methods.o: headers/Eigen/src/Core/util/Memory.h
headers/Methods.o: headers/Eigen/src/Core/NumTraits.h
headers/Methods.o: headers/Eigen/src/Core/MathFunctions.h
headers/Methods.o: headers/Eigen/src/Core/GenericPacketMath.h
headers/Methods.o: headers/Eigen/src/Core/MathFunctionsImpl.h
headers/Methods.o: headers/Eigen/src/Core/arch/Default/ConjHelper.h
headers/Methods.o: headers/Eigen/src/Core/arch/CUDA/Half.h
headers/Methods.o: headers/Eigen/src/Core/arch/CUDA/PacketMathHalf.h
headers/Methods.o: headers/Eigen/src/Core/arch/CUDA/TypeCasting.h
headers/Methods.o: headers/Eigen/src/Core/arch/Default/Settings.h
headers/Methods.o: headers/Eigen/src/Core/functors/TernaryFunctors.h
headers/Methods.o: headers/Eigen/src/Core/functors/BinaryFunctors.h
headers/Methods.o: headers/Eigen/src/Core/functors/UnaryFunctors.h
headers/Methods.o: headers/Eigen/src/Core/functors/NullaryFunctors.h
headers/Methods.o: headers/Eigen/src/Core/functors/StlFunctors.h
headers/Methods.o: headers/Eigen/src/Core/functors/AssignmentFunctors.h
headers/Methods.o: headers/Eigen/src/Core/arch/CUDA/Complex.h
headers/Methods.o: headers/Eigen/src/Core/IO.h
headers/Methods.o: headers/Eigen/src/Core/DenseCoeffsBase.h
headers/Methods.o: headers/Eigen/src/Core/DenseBase.h
headers/Methods.o: headers/Eigen/src/plugins/BlockMethods.h
headers/Methods.o: headers/Eigen/src/Core/MatrixBase.h
headers/Methods.o: headers/Eigen/src/plugins/CommonCwiseUnaryOps.h
headers/Methods.o: headers/Eigen/src/plugins/CommonCwiseBinaryOps.h
headers/Methods.o: headers/Eigen/src/plugins/MatrixCwiseUnaryOps.h
headers/Methods.o: headers/Eigen/src/plugins/MatrixCwiseBinaryOps.h
headers/Methods.o: headers/Eigen/src/Core/EigenBase.h
headers/Methods.o: headers/Eigen/src/Core/Product.h
headers/Methods.o: headers/Eigen/src/Core/CoreEvaluators.h
headers/Methods.o: headers/Eigen/src/Core/AssignEvaluator.h
headers/Methods.o: headers/Eigen/src/Core/Assign.h
headers/Methods.o: headers/Eigen/src/Core/ArrayBase.h
headers/Methods.o: headers/Eigen/src/plugins/ArrayCwiseUnaryOps.h
headers/Methods.o: headers/Eigen/src/plugins/ArrayCwiseBinaryOps.h
headers/Methods.o: headers/Eigen/src/Core/util/BlasUtil.h
headers/Methods.o: headers/Eigen/src/Core/DenseStorage.h
headers/Methods.o: headers/Eigen/src/Core/NestByValue.h
headers/Methods.o: headers/Eigen/src/Core/ReturnByValue.h
headers/Methods.o: headers/Eigen/src/Core/NoAlias.h
headers/Methods.o: headers/Eigen/src/Core/PlainObjectBase.h
headers/Methods.o: headers/Eigen/src/Core/Matrix.h
headers/Methods.o: headers/Eigen/src/Core/Array.h
headers/Methods.o: headers/Eigen/src/Core/CwiseTernaryOp.h
headers/Methods.o: headers/Eigen/src/Core/CwiseBinaryOp.h
headers/Methods.o: headers/Eigen/src/Core/CwiseUnaryOp.h
headers/Methods.o: headers/Eigen/src/Core/CwiseNullaryOp.h
headers/Methods.o: headers/Eigen/src/Core/CwiseUnaryView.h
headers/Methods.o: headers/Eigen/src/Core/SelfCwiseBinaryOp.h
headers/Methods.o: headers/Eigen/src/Core/Dot.h
headers/Methods.o: headers/Eigen/src/Core/StableNorm.h
headers/Methods.o: headers/Eigen/src/Core/Stride.h
headers/Methods.o: headers/Eigen/src/Core/MapBase.h
headers/Methods.o: headers/Eigen/src/Core/Map.h headers/Eigen/src/Core/Ref.h
headers/Methods.o: headers/Eigen/src/Core/Block.h
headers/Methods.o: headers/Eigen/src/Core/VectorBlock.h
headers/Methods.o: headers/Eigen/src/Core/Transpose.h
headers/Methods.o: headers/Eigen/src/Core/DiagonalMatrix.h
headers/Methods.o: headers/Eigen/src/Core/Diagonal.h
headers/Methods.o: headers/Eigen/src/Core/DiagonalProduct.h
headers/Methods.o: headers/Eigen/src/Core/Redux.h
headers/Methods.o: headers/Eigen/src/Core/Visitor.h
headers/Methods.o: headers/Eigen/src/Core/Fuzzy.h
headers/Methods.o: headers/Eigen/src/Core/Swap.h
headers/Methods.o: headers/Eigen/src/Core/CommaInitializer.h
headers/Methods.o: headers/Eigen/src/Core/GeneralProduct.h
headers/Methods.o: headers/Eigen/src/Core/Solve.h
headers/Methods.o: headers/Eigen/src/Core/Inverse.h
headers/Methods.o: headers/Eigen/src/Core/SolverBase.h
headers/Methods.o: headers/Eigen/src/Core/PermutationMatrix.h
headers/Methods.o: headers/Eigen/src/Core/Transpositions.h
headers/Methods.o: headers/Eigen/src/Core/TriangularMatrix.h
headers/Methods.o: headers/Eigen/src/Core/SelfAdjointView.h
headers/Methods.o: headers/Eigen/src/Core/products/GeneralBlockPanelKernel.h
headers/Methods.o: headers/Eigen/src/Core/products/Parallelizer.h
headers/Methods.o: headers/Eigen/src/Core/ProductEvaluators.h
headers/Methods.o: headers/Eigen/src/Core/products/GeneralMatrixVector.h
headers/Methods.o: headers/Eigen/src/Core/products/GeneralMatrixMatrix.h
headers/Methods.o: headers/Eigen/src/Core/SolveTriangular.h
headers/Methods.o: headers/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h
headers/Methods.o: headers/Eigen/src/Core/products/SelfadjointMatrixVector.h
headers/Methods.o: headers/Eigen/src/Core/products/SelfadjointMatrixMatrix.h
headers/Methods.o: headers/Eigen/src/Core/products/SelfadjointProduct.h
headers/Methods.o: headers/Eigen/src/Core/products/SelfadjointRank2Update.h
headers/Methods.o: headers/Eigen/src/Core/products/TriangularMatrixVector.h
headers/Methods.o: headers/Eigen/src/Core/products/TriangularMatrixMatrix.h
headers/Methods.o: headers/Eigen/src/Core/products/TriangularSolverMatrix.h
headers/Methods.o: headers/Eigen/src/Core/products/TriangularSolverVector.h
headers/Methods.o: headers/Eigen/src/Core/BandMatrix.h
headers/Methods.o: headers/Eigen/src/Core/CoreIterators.h
headers/Methods.o: headers/Eigen/src/Core/ConditionEstimator.h
headers/Methods.o: headers/Eigen/src/Core/BooleanRedux.h
headers/Methods.o: headers/Eigen/src/Core/Select.h
headers/Methods.o: headers/Eigen/src/Core/VectorwiseOp.h
headers/Methods.o: headers/Eigen/src/Core/Random.h
headers/Methods.o: headers/Eigen/src/Core/Replicate.h
headers/Methods.o: headers/Eigen/src/Core/Reverse.h
headers/Methods.o: headers/Eigen/src/Core/ArrayWrapper.h
headers/Methods.o: headers/Eigen/src/Core/GlobalFunctions.h
headers/Methods.o: headers/Eigen/src/Core/util/ReenableStupidWarnings.h
headers/Methods.o: headers/Eigen/LU headers/Eigen/src/misc/Kernel.h
headers/Methods.o: headers/Eigen/src/misc/Image.h
headers/Methods.o: headers/Eigen/src/LU/FullPivLU.h
headers/Methods.o: headers/Eigen/src/LU/PartialPivLU.h
headers/Methods.o: headers/Eigen/src/LU/Determinant.h
headers/Methods.o: headers/Eigen/src/LU/InverseImpl.h headers/Eigen/Cholesky
headers/Methods.o: headers/Eigen/Jacobi headers/Eigen/src/Jacobi/Jacobi.h
headers/Methods.o: headers/Eigen/src/Cholesky/LLT.h
headers/Methods.o: headers/Eigen/src/Cholesky/LDLT.h headers/Eigen/QR
headers/Methods.o: headers/Eigen/Householder
headers/Methods.o: headers/Eigen/src/Householder/Householder.h
headers/Methods.o: headers/Eigen/src/Householder/HouseholderSequence.h
headers/Methods.o: headers/Eigen/src/Householder/BlockHouseholder.h
headers/Methods.o: headers/Eigen/src/QR/HouseholderQR.h
headers/Methods.o: headers/Eigen/src/QR/FullPivHouseholderQR.h
headers/Methods.o: headers/Eigen/src/QR/ColPivHouseholderQR.h
headers/Methods.o: headers/Eigen/src/QR/CompleteOrthogonalDecomposition.h
headers/Methods.o: headers/Eigen/SVD headers/Eigen/src/misc/RealSvd2x2.h
headers/Methods.o: headers/Eigen/src/SVD/UpperBidiagonalization.h
headers/Methods.o: headers/Eigen/src/SVD/SVDBase.h
headers/Methods.o: headers/Eigen/src/SVD/JacobiSVD.h
headers/Methods.o: headers/Eigen/src/SVD/BDCSVD.h headers/Eigen/Geometry
headers/Methods.o: headers/Eigen/src/Geometry/OrthoMethods.h
headers/Methods.o: headers/Eigen/src/Geometry/EulerAngles.h
headers/Methods.o: headers/Eigen/src/Geometry/Homogeneous.h
headers/Methods.o: headers/Eigen/src/Geometry/RotationBase.h
headers/Methods.o: headers/Eigen/src/Geometry/Rotation2D.h
headers/Methods.o: headers/Eigen/src/Geometry/Quaternion.h
headers/Methods.o: headers/Eigen/src/Geometry/AngleAxis.h
headers/Methods.o: headers/Eigen/src/Geometry/Transform.h
headers/Methods.o: headers/Eigen/src/Geometry/Translation.h
headers/Methods.o: headers/Eigen/src/Geometry/Scaling.h
headers/Methods.o: headers/Eigen/src/Geometry/Hyperplane.h
headers/Methods.o: headers/Eigen/src/Geometry/ParametrizedLine.h
headers/Methods.o: headers/Eigen/src/Geometry/AlignedBox.h
headers/Methods.o: headers/Eigen/src/Geometry/Umeyama.h
headers/Methods.o: headers/Eigen/Eigenvalues
headers/Methods.o: headers/Eigen/src/Eigenvalues/Tridiagonalization.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/RealSchur.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/./HessenbergDecomposition.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/EigenSolver.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/./RealSchur.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/./Tridiagonalization.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/HessenbergDecomposition.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/ComplexSchur.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/ComplexEigenSolver.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/./ComplexSchur.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/RealQZ.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/./RealQZ.h
headers/Methods.o: headers/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h
