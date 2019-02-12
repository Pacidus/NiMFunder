/*On importe nos Lib*/
#include "headers/Methods.hpp"

int main()
{
	ostringstream filename;
	filename << "../Results/AIRSS.res";
	ofstream file(filename.str().c_str());

	double Cout;	 		//Le cout de la recherche

	int Rep;				//Nombre de répétition

	double Hmin, Hrange;	//Hauteur minimale & range
	double Lmin, Lrange;	//Largeur minimale & range
	double Smin, Srange;	//Sigma   minimale & range

	ArrayXd val(Rep);		//Valeurs stockées

	double epsilon = 1e-9;

	SolNim Release;							//Système

	for(int D = 1; D <= 100; D++)			//Pour les dimension de 1 à 100
	{
		for(int N = 1; N <= 16; N*=2)		//Pour 1 à 16 minimas
		{
			for(int i = 0; i < Rep; i++)	//On répète jusqu'à Rep
			{
				GenPays(D, N, Lmin, Lrange, Smin, Srange, Hmin, Hrange);
				/*On génère le paysage*/
				Cout = 0;							//Le cout initial et nul
				Release.Init();						//On initialise le paysage
				while(1 == i)		//
				{
					Cout += Steepdes(Release, Lmin, Lrange, epsilon);
				}
				val(i) = Cout;
			}
		}
	}

	return 0;
}
