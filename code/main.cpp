/* on importe nos Lib*/
#include "headers/Methods.hpp"

void FindStep(int& CoutR,int& CoutG, double L, double L0)
{
	double epsilon = 1e-5;
	double l = 4;
	double e = 4;
	int Nbr = 10000;
	double val[Nbr];
	VectorXd p[Nbr];
	CoutG = Nbr;
	CoutG +=	Gabriele(p, val, Nbr, L, L0, e, l, epsilon);

	SolNim Release;
	Release.Init();
	for(int i = 1; i < Nbr; i++)
	{
		CoutR += Steepdes(Release, L, L0, epsilon);
	}
}



int main()
{
	ostringstream filename;
	filename << "/home/yohan/Bureau/NiMFunder/Results/FindStep.res";
	ofstream file(filename.str().c_str());

	int CoutR;
	int CoutG;
	int D = 2;
	double L = 100;
	double L0 = 0;
	double sigma = 3.5;
	double H = -4;
	double H0 = -12;
	double sigma0 = 0.5;
	int CTotR = 0;
	int CTotG = 0;
	for(int i = 0; i < 100; i+=10)
	{
		CTotR = 0;
		CTotG = 0;
		GenPays(D, 1+i, L, L0, sigma, sigma0, H, H0);
		for(int j = 0; j < (1000-i*10); j++)
		{
			CoutR = 0;
			CoutG = 0;
			cout << j/10. << endl;
			FindStep(CoutR, CoutG, L, L0);
			CTotR += CoutR;
			CTotG += CoutG;
    	}

    	file << (i+1) << " " << CTotR/(1000.-i*10) << " " << CTotG/(1000.-i*10) << endl;
    }

    return 0;
}
