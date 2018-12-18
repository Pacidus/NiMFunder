CC      = g++
CFLAGS  = -Ofast -W -Wall -Wextra -DNDEBUG

DEPS    = headers/Class.hpp
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

Map.o: headers/Class.hpp headers/Eigen/Dense headers/Eigen/Core
Map.o: headers/Eigen/src/Core/util/DisableStupidWarnings.h
Map.o: headers/Eigen/src/Core/util/Macros.h
Map.o: headers/Eigen/src/Core/util/MKL_support.h
Map.o: headers/Eigen/src/Core/util/Constants.h
Map.o: headers/Eigen/src/Core/util/Meta.h
Map.o: headers/Eigen/src/Core/util/ForwardDeclarations.h
Map.o: headers/Eigen/src/Core/util/StaticAssert.h
Map.o: headers/Eigen/src/Core/util/XprHelper.h
Map.o: headers/Eigen/src/Core/util/Memory.h
Map.o: headers/Eigen/src/Core/NumTraits.h
Map.o: headers/Eigen/src/Core/MathFunctions.h
Map.o: headers/Eigen/src/Core/GenericPacketMath.h
Map.o: headers/Eigen/src/Core/MathFunctionsImpl.h
Map.o: headers/Eigen/src/Core/arch/Default/ConjHelper.h
Map.o: headers/Eigen/src/Core/arch/CUDA/Half.h
Map.o: headers/Eigen/src/Core/arch/CUDA/PacketMathHalf.h
Map.o: headers/Eigen/src/Core/arch/CUDA/TypeCasting.h
Map.o: headers/Eigen/src/Core/arch/Default/Settings.h
Map.o: headers/Eigen/src/Core/functors/TernaryFunctors.h
Map.o: headers/Eigen/src/Core/functors/BinaryFunctors.h
Map.o: headers/Eigen/src/Core/functors/UnaryFunctors.h
Map.o: headers/Eigen/src/Core/functors/NullaryFunctors.h
Map.o: headers/Eigen/src/Core/functors/StlFunctors.h
Map.o: headers/Eigen/src/Core/functors/AssignmentFunctors.h
Map.o: headers/Eigen/src/Core/arch/CUDA/Complex.h headers/Eigen/src/Core/IO.h
Map.o: headers/Eigen/src/Core/DenseCoeffsBase.h
Map.o: headers/Eigen/src/Core/DenseBase.h
Map.o: headers/Eigen/src/plugins/BlockMethods.h
Map.o: headers/Eigen/src/Core/MatrixBase.h
Map.o: headers/Eigen/src/plugins/CommonCwiseUnaryOps.h
Map.o: headers/Eigen/src/plugins/CommonCwiseBinaryOps.h
Map.o: headers/Eigen/src/plugins/MatrixCwiseUnaryOps.h
Map.o: headers/Eigen/src/plugins/MatrixCwiseBinaryOps.h
Map.o: headers/Eigen/src/Core/EigenBase.h headers/Eigen/src/Core/Product.h
Map.o: headers/Eigen/src/Core/CoreEvaluators.h
Map.o: headers/Eigen/src/Core/AssignEvaluator.h
Map.o: headers/Eigen/src/Core/Assign.h headers/Eigen/src/Core/ArrayBase.h
Map.o: headers/Eigen/src/plugins/ArrayCwiseUnaryOps.h
Map.o: headers/Eigen/src/plugins/ArrayCwiseBinaryOps.h
Map.o: headers/Eigen/src/Core/util/BlasUtil.h
Map.o: headers/Eigen/src/Core/DenseStorage.h
Map.o: headers/Eigen/src/Core/NestByValue.h
Map.o: headers/Eigen/src/Core/ReturnByValue.h
Map.o: headers/Eigen/src/Core/NoAlias.h
Map.o: headers/Eigen/src/Core/PlainObjectBase.h
Map.o: headers/Eigen/src/Core/Matrix.h headers/Eigen/src/Core/Array.h
Map.o: headers/Eigen/src/Core/CwiseTernaryOp.h
Map.o: headers/Eigen/src/Core/CwiseBinaryOp.h
Map.o: headers/Eigen/src/Core/CwiseUnaryOp.h
Map.o: headers/Eigen/src/Core/CwiseNullaryOp.h
Map.o: headers/Eigen/src/Core/CwiseUnaryView.h
Map.o: headers/Eigen/src/Core/SelfCwiseBinaryOp.h
Map.o: headers/Eigen/src/Core/Dot.h headers/Eigen/src/Core/StableNorm.h
Map.o: headers/Eigen/src/Core/Stride.h headers/Eigen/src/Core/MapBase.h
Map.o: headers/Eigen/src/Core/Map.h headers/Eigen/src/Core/Ref.h
Map.o: headers/Eigen/src/Core/Block.h headers/Eigen/src/Core/VectorBlock.h
Map.o: headers/Eigen/src/Core/Transpose.h
Map.o: headers/Eigen/src/Core/DiagonalMatrix.h
Map.o: headers/Eigen/src/Core/Diagonal.h
Map.o: headers/Eigen/src/Core/DiagonalProduct.h
Map.o: headers/Eigen/src/Core/Redux.h headers/Eigen/src/Core/Visitor.h
Map.o: headers/Eigen/src/Core/Fuzzy.h headers/Eigen/src/Core/Swap.h
Map.o: headers/Eigen/src/Core/CommaInitializer.h
Map.o: headers/Eigen/src/Core/GeneralProduct.h headers/Eigen/src/Core/Solve.h
Map.o: headers/Eigen/src/Core/Inverse.h headers/Eigen/src/Core/SolverBase.h
Map.o: headers/Eigen/src/Core/PermutationMatrix.h
Map.o: headers/Eigen/src/Core/Transpositions.h
Map.o: headers/Eigen/src/Core/TriangularMatrix.h
Map.o: headers/Eigen/src/Core/SelfAdjointView.h
Map.o: headers/Eigen/src/Core/products/GeneralBlockPanelKernel.h
Map.o: headers/Eigen/src/Core/products/Parallelizer.h
Map.o: headers/Eigen/src/Core/ProductEvaluators.h
Map.o: headers/Eigen/src/Core/products/GeneralMatrixVector.h
Map.o: headers/Eigen/src/Core/products/GeneralMatrixMatrix.h
Map.o: headers/Eigen/src/Core/SolveTriangular.h
Map.o: headers/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h
Map.o: headers/Eigen/src/Core/products/SelfadjointMatrixVector.h
Map.o: headers/Eigen/src/Core/products/SelfadjointMatrixMatrix.h
Map.o: headers/Eigen/src/Core/products/SelfadjointProduct.h
Map.o: headers/Eigen/src/Core/products/SelfadjointRank2Update.h
Map.o: headers/Eigen/src/Core/products/TriangularMatrixVector.h
Map.o: headers/Eigen/src/Core/products/TriangularMatrixMatrix.h
Map.o: headers/Eigen/src/Core/products/TriangularSolverMatrix.h
Map.o: headers/Eigen/src/Core/products/TriangularSolverVector.h
Map.o: headers/Eigen/src/Core/BandMatrix.h
Map.o: headers/Eigen/src/Core/CoreIterators.h
Map.o: headers/Eigen/src/Core/ConditionEstimator.h
Map.o: headers/Eigen/src/Core/BooleanRedux.h headers/Eigen/src/Core/Select.h
Map.o: headers/Eigen/src/Core/VectorwiseOp.h headers/Eigen/src/Core/Random.h
Map.o: headers/Eigen/src/Core/Replicate.h headers/Eigen/src/Core/Reverse.h
Map.o: headers/Eigen/src/Core/ArrayWrapper.h
Map.o: headers/Eigen/src/Core/GlobalFunctions.h
Map.o: headers/Eigen/src/Core/util/ReenableStupidWarnings.h headers/Eigen/LU
Map.o: headers/Eigen/src/misc/Kernel.h headers/Eigen/src/misc/Image.h
Map.o: headers/Eigen/src/LU/FullPivLU.h headers/Eigen/src/LU/PartialPivLU.h
Map.o: headers/Eigen/src/LU/Determinant.h headers/Eigen/src/LU/InverseImpl.h
Map.o: headers/Eigen/Cholesky headers/Eigen/Jacobi
Map.o: headers/Eigen/src/Jacobi/Jacobi.h headers/Eigen/src/Cholesky/LLT.h
Map.o: headers/Eigen/src/Cholesky/LDLT.h headers/Eigen/QR
Map.o: headers/Eigen/Householder headers/Eigen/src/Householder/Householder.h
Map.o: headers/Eigen/src/Householder/HouseholderSequence.h
Map.o: headers/Eigen/src/Householder/BlockHouseholder.h
Map.o: headers/Eigen/src/QR/HouseholderQR.h
Map.o: headers/Eigen/src/QR/FullPivHouseholderQR.h
Map.o: headers/Eigen/src/QR/ColPivHouseholderQR.h
Map.o: headers/Eigen/src/QR/CompleteOrthogonalDecomposition.h
Map.o: headers/Eigen/SVD headers/Eigen/src/misc/RealSvd2x2.h
Map.o: headers/Eigen/src/SVD/UpperBidiagonalization.h
Map.o: headers/Eigen/src/SVD/SVDBase.h headers/Eigen/src/SVD/JacobiSVD.h
Map.o: headers/Eigen/src/SVD/BDCSVD.h headers/Eigen/Geometry
Map.o: headers/Eigen/src/Geometry/OrthoMethods.h
Map.o: headers/Eigen/src/Geometry/EulerAngles.h
Map.o: headers/Eigen/src/Geometry/Homogeneous.h
Map.o: headers/Eigen/src/Geometry/RotationBase.h
Map.o: headers/Eigen/src/Geometry/Rotation2D.h
Map.o: headers/Eigen/src/Geometry/Quaternion.h
Map.o: headers/Eigen/src/Geometry/AngleAxis.h
Map.o: headers/Eigen/src/Geometry/Transform.h
Map.o: headers/Eigen/src/Geometry/Translation.h
Map.o: headers/Eigen/src/Geometry/Scaling.h
Map.o: headers/Eigen/src/Geometry/Hyperplane.h
Map.o: headers/Eigen/src/Geometry/ParametrizedLine.h
Map.o: headers/Eigen/src/Geometry/AlignedBox.h
Map.o: headers/Eigen/src/Geometry/Umeyama.h headers/Eigen/Eigenvalues
Map.o: headers/Eigen/src/Eigenvalues/Tridiagonalization.h
Map.o: headers/Eigen/src/Eigenvalues/RealSchur.h
Map.o: headers/Eigen/src/Eigenvalues/./HessenbergDecomposition.h
Map.o: headers/Eigen/src/Eigenvalues/EigenSolver.h
Map.o: headers/Eigen/src/Eigenvalues/./RealSchur.h
Map.o: headers/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h
Map.o: headers/Eigen/src/Eigenvalues/./Tridiagonalization.h
Map.o: headers/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h
Map.o: headers/Eigen/src/Eigenvalues/HessenbergDecomposition.h
Map.o: headers/Eigen/src/Eigenvalues/ComplexSchur.h
Map.o: headers/Eigen/src/Eigenvalues/ComplexEigenSolver.h
Map.o: headers/Eigen/src/Eigenvalues/./ComplexSchur.h
Map.o: headers/Eigen/src/Eigenvalues/RealQZ.h
Map.o: headers/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h
Map.o: headers/Eigen/src/Eigenvalues/./RealQZ.h
Map.o: headers/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h
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
