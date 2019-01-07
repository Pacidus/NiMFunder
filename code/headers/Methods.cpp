/* on importe nos Lib*/
#include "Methods.hpp"

int Gabriele(VectorXd p[], double val[], int& Nbr, double L, double L0, double e, double l, double epsilon)
{
	map(p, val, Nbr, L, L0);
	cout << "==== TEST ====\n"<< endl;
	reject(p, val, Nbr, e);
	cout << "==== TEST ====\n"<< endl;
	toclose(p, val, Nbr, l);
	cout << "==== TEST ====\n"<< endl;
	int couts;
	couts = release(p, val, Nbr, epsilon);
	cout << "==== TEST ====\n"<< endl;
	return couts;
}

void map(VectorXd p[], double val[], int Nbr, double L, double L0)
{
	SolNim Map;
	Map.Init();
	int D;
	Map.get_Dim(D);
	VectorXd p1(D);
	for(int i = 0; i < Nbr; i++)
	{
		Map.Rpos(L,L0);
		Map.get_pos(p1);
		p[i] = p1;
		val[i] = Map.landscape();
	}
}

void reject(VectorXd p[], double val[], int& Nbr, double e)
{
	double min = val[0];
	int imin = 0;
	VectorXd palt[Nbr];
	double valt[Nbr];
	int N = 0;

	for(int i = 1; i < Nbr; i++)
	{
		if(min > val[i])
		{
			imin = i;
			min = val[i];
		}
	}

	for(int i = 0; i < Nbr; i++)
	{
		if(val[i] <= min+e)
		{
			palt[N] = p[i];
			valt[N] = val[i];
			N++;
		}
	}

	Nbr = N;
	VectorXd palt2[Nbr];
	double valt2[Nbr];

	for(int i = 0; i < Nbr; i++)
	{
		palt2[i] = palt[i];
		valt2[i] = valt[i];
	}
	p = palt2;
	val = valt2;
}

void toclose(VectorXd p[], double val[], int& Nbr, double l)
{

	VectorXd palt[Nbr];
	double valt[Nbr];
	int N = 0;
	bool k = 1;
	for(int i = 0; i < Nbr; i++)
	{
		k = 1;
		for(int j = 0; j < Nbr; j ++)
		{
			if((p[i]-p[j]).norm() < l)
			{
				if(val[i] > val[j])
				{
					k = 0;
					break;
				}
			}
		}
		if(k)
		{
			palt[N] = p[i];
			valt[N] = val[i];
			N++;
		}
	}

	Nbr = N;
	VectorXd palt2[Nbr];
	double valt2[Nbr];

	for(int i = 0; i < Nbr; i++)
	{
		palt2[i] = palt[i];
		valt2[i] = valt[i];
	}
	p = palt2;
	val = valt2;
}

int release(VectorXd p[], double val[], int Nbr, double epsilon)
{
	int Ntot = 0;
	SolNim Release;
	Release.Init();

	int D;
	Release.get_Dim(D);
	VectorXd b(D);
	double dt = 1;
	int j;

	VectorXd p1(D);

	for(int i = 1; i < Nbr; i++)
	{
		j = 1;
		p1 = p[i];
		Release.set_pos(p1);

		do
		{
			dt = 2./j;
			Release.SteepDescent(dt);
			Release.get_b(b);
			j++;
		}
		while(b.norm() > epsilon and j > 1001);

		Ntot += (j-1);
	}

	return Ntot;
}

void GenPays(int D, int N, double L, double L0, double sigma, double sigma0, double H, double H0)
{
	ostringstream filename;
	filename << "/home/yohan/Bureau/NiMFunder/Paysage.init";
	ofstream file(filename.str().c_str());

	file << N << " " << D << endl;
	MatrixXd P(D,N);
	MatrixXd Sigmas(D,N);
	VectorXd Hs(N);

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < D; j++)
		{
			P(j,i) = L*(rand()/(RAND_MAX+1.0)) + L0;
			Sigmas(j,i) = sigma*(rand()/(RAND_MAX+1.0)) + sigma0;
		}
		Hs(i) = H*(rand()/(RAND_MAX+1.0)) + H0;
	}

	file << P << endl;
	file << Sigmas << endl;
	file << Hs << endl;
}

int Steepdes(SolNim Release, double L, double L0, double epsilon)
{
	int j = 1;
	Release.Rpos(L,L0);

	int D;
	Release.get_Dim(D);
	VectorXd b(D);
	double dt = 1;

	do
	{
		dt = 2./j;
		Release.SteepDescent(dt);
		Release.get_b(b);
		j++;
	}
	while(b.norm() > epsilon and j > 1001);

	return j-1;
}
