/* on importe nos Lib*/
#include "headers/Class.hpp"

int main()
{
	SolNim PayF;
	PayF.Init();
	int D;
	int N;

	PayF.get_Dim(D);
	PayF.get_Num(N);
	VectorXd p(D);
	ostringstream filename;

	filename << "../Results/Min" << N << "Dim" << D << "Map.csv";

	ofstream file(filename.str().c_str());

	int n = 60;
	double a, c;
	a = 100./n;
	c = 0;

	for(int j = 0; j <= n; j++)
	{
		p(0) = a*j+c;
		for(int i = 0; i <= n; i++)
		{
			p(1) = a*i+c;
			PayF.set_pos(p);
			PayF.Rpos(a,c);
			PayF.get_pos(p);
			file << p.transpose() << " " << PayF.landscape() << endl;
		}
		file << endl;
	}
	file.close();
    return 0;
}
