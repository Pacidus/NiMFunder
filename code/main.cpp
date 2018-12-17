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
	VectorXd b(D);
	ostringstream filename;

	filename << "/home/yohan/Bureau/NiMFunder/Results/Min" << N << "Dim" << D << ".res";

	ofstream file(filename.str().c_str());

	int n;
	cin >> n;
	double epsilon = 1e-6;
	int i = 0;
	double dt = 1;
	for(int j = 0; j < n; j++)
	{
		PayF.Rpos(1.5,.25);
		i = 1;

   		do
    	{
    	dt = 0.0001;
		i++; 				//On implémente de 1

		PayF.SteepDescent(dt);
		PayF.get_b(b);
		PayF.get_pos(p);
		i++;
		}
	    while(b.norm() > epsilon and i < 1000000);
		file << p.transpose() << " " << PayF.landscape() << endl;
	}
    return 0;
}
