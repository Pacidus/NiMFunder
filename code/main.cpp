/*On importe nos Lib*/
#include "headers/Methods.hpp"

void FindStep(int& CoutR,int& CoutG, double L, double L0)
{
	double epsilon = 1e-9;
	double l = 6;
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
	filename << "../Results/FindStep.csv";
	ofstream file(filename.str().c_str());

	int CoutR;
	int CoutG;
	int D = 2;
	double L = 100;
	double L0 = 0;
	double sigma = 5.5;
	double sigma0 = 1.;
	double H = -4;
	double H0 = -5;
	double CTotR = 0;
	double CTotG = 0;
	ArrayXd Rval(100);
	ArrayXd Gval(100);

	for(int i = 1; i < 100; i+=1)
	{
		for(int j = 3; j < 20; j+=3)
		{
			CTotR = 0;
			CTotG = 0;

			for(int k = 0; k < 100; k++)
			{
				GenPays(i, j, L, L0, sigma, sigma0, H, H0);
				CoutR = 0;
				CoutG = 0;
				FindStep(CoutR, CoutG, L, L0);
				CTotR += CoutR;
				CTotG += CoutG;
				Rval(k) = CoutR;
				Gval(k) = CoutG;
    		}

    		file << i << " " << j << " " << CTotR / 100 << " ";
    		file << CTotG/100 << " ";
    		file << sqrt(((Rval-(CTotR/100))*(Rval-(CTotR/100))).sum()/100);
    		file << " ";
    		file << sqrt(((Gval-(CTotG/100))*(Gval-(CTotG/100))).sum()/100);
    		file << endl;


    	}
    }
	file.close();
    return 0;
}
