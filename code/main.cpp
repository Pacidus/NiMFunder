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
	VectorXd b(D);
	I1 << 1,0;
	I2 << 0,1;
	I3 << 1,1;
	ostringstream filename;

	filename << "/home/yohan/Bureau/NiMFunder/Results/Min" << N << "Dim" << D << "Solve.res";

	ofstream file(filename.str().c_str());

	int n;
	cin >> n;
	double epsilon = 1e-6;
	int i = 0;

	for(int j = 0; j < n; j++)
	{
		PayF.Rpos(1.5,.25);
		int i = 0;

   		do
    	{
		i++; 				//On implémente de 1
		PayF.SteepDescent();
		PayF.get_b(b);
		PayF.get_pos(p);
		i++;
		}
	    while(b.norm() > epsilon and i < 1000000);

		file << p.transpose() << " " << PayF.landscape() << endl;
	}
    return 0;
}
/*int D = p.size();
    VectorXd b(D);
    double df2;*/

