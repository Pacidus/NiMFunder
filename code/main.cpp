/* on importe nos Lib*/
#include "headers/Methods.hpp"

void FindStep(int& NRe, int& NRa, int N_max)
{
	ostringstream filename;
	filename << "/home/yohan/Bureau/NiMFunder/Results/FindStep.res";
	ofstream file(filename.str().c_str());

	if(NRe+NRa == 0)
	{
		ofstream file(filename.str().c_str());
	}
	else
	{
		fstream file;
		file.open(filename.str().c_str(), fstream::out | fstream::app);
	}

	int Ntot = 0;
	double epsilon = 1e-15;
	double l = 0.3;
	double e = 3;
	double L0 = 0;
	double L = 2;
	int Nbr = 100;
	double val[Nbr];
	VectorXd p[Nbr];
	double couts = Nbr;
	couts +=	Gabriele(p, val, Nbr, L, L0, e, l, epsilon);

	SolNim Release;
	Release.Init();

	int j;

	int D;
	Release.get_Dim(D);
	VectorXd b(D);

	double dt = 1;

	for(int i = 1; i < Nbr; i++)
	{
		j = 1;
		Release.Rpos(L,L0);

		do
		{
			dt = 0.2/j;
			Release.SteepDescent(dt);
			Release.get_b(b);
			j++;
		}
		while(b.norm() < epsilon);
		Ntot += j;
	}
	cout << Ntot <<" "<< couts << endl;
}



int main()
{
	int NRe = 0;
	int NRa = 0;
	FindStep(NRe, NRa, 10);
    return 0;
}
