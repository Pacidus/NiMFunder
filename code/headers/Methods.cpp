/* on importe nos Lib*/
#include "Methods.hpp"

int Gabriele(VectorXd p[], double val[], int& Nbr, double L, double L0, double e, double l, double epsilon)
{
	map(p, val, Nbr, L, L0);
	cout << p[0] << " " << val[0]<<" "<< Nbr << endl;
	reject(p, val, Nbr, e);
	cout << p[0] << " " << val[0]<<" "<< Nbr << endl;
	toclose(p, val, Nbr, l);
	cout << p[0] << " " << val[0]<<" "<< Nbr << endl;
	int couts;
	couts = release(p, val, Nbr, epsilon);
	cout << p[0] << " " << val[0]<<" "<< Nbr << endl;
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
	int j;

	int D;
	Release.get_Dim(D);
	VectorXd b(D);

	double dt = 1;

	for(int i = 1; i < Nbr; i++)
	{
		j = 1;
		Release.set_pos(p[i]);

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

	return Ntot;
}
