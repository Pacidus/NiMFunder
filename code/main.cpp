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
	double epsilon = 1e-5;
	int i = 0;
	double dt = 1;
	for(int j = 0; j < n; j++)
	{
		PayF.Rpos(1.5,.25);
		i = 1;

   		do
    	{
    	dt = 0.1/i;
		i++; 				//On implémente de 1
		PayF.SteepDescent(dt);
		PayF.get_b(b);
		PayF.get_pos(p);
		file << p.transpose() << " " << PayF.landscape() << endl;
		}
	    while(b.norm()*b.norm() > epsilon*epsilon and i < 1000);

	}
    return 0;
}
