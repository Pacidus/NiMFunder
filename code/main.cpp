/* on importe nos Lib*/
#include "headers/Class.hpp"

int main()
{
	SolNim PayF;
	PayF.Init();

	int N;
	PayF.get_Num(N);
	int D;
	PayF.get_Dim(D);
	MatrixXd Ps(N,D);
	PayF.get_Pos(Ps);
	MatrixXd Sigma(N,D);
	PayF.get_Sig(Ps);
	VectorXd H(N);
	PayF.get_Hau(H);
	VectorXd p(D);
	VectorXi I1(D);
	VectorXi I2(D);
	VectorXi I3(D);
	I1 << 1,0;
	I2 << 0,1;
	I3 << 1,1;
	ostringstream filename;

	filename << "/home/yohan/Bureau/NiMFunder/Results/Min" << N << "Dim" << D << ".res";

	ofstream file(filename.str().c_str());

	int n;
	cin >> n;

	for(int j = 0; j < n; j++)
	{
		PayF.Rpos(2,0);
		PayF.get_pos(p);
		file << p.transpose() << " " << PayF.landscape() << " ";
		file << PayF.dlandscape(I1) << " " << PayF.dlandscape(I2) << " ";
		file << PayF.dlandscape(I3) << endl;
	}
    return 0;
}
