/*On importe nos Lib*/
#include "headers/Methods.hpp"

double Proba(VectorXi T)
/* 	Proba: Renvoie la proba de trouver un nouveau minima en
	fonction du nombre de minimas trouvées 0 <= Proba(T) <= 1
	P(N|T) = (dim(T)/(||T||_n))^r */
{
	cout(T.rows());
	return T.rows()/T.sum();
}

bool Criterion(VectorXi T, double Veps)
/*Criterion: Critère qui définie si il est intérésant de chercher un minima*/
{
	return (Veps < Proba(T));
}

bool Equiv(double x, double y, double Sign)
/*Equiv: Retourne true si x = y à une précision près*/
{
	int X = (int)(x * Sign);
	int Y = (int)(y * Sign);
	return (X == Y);
}

void New(MatrixXd& Min, VectorXi& T, VectorXd p)
/*New: verifie si p n'est pas dans la Matrice Min*/
{
	int D = Min.cols();
	int N = Min.rows();

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < D; j++)
		{
			if(not(Equiv(M(j,i),p(j),100.)))
			{
				break;
			}
			if(j == (D-1))
			{
				T.row(i) += 1;
				return;
			}
		}
	}
	T.conservativeResize(N+1);
	Min.conservativeResize(D,N+1);
	T.row(N) = 1;
	Min.row(N) = p.transpose();
	return;
}

int main()
{
	ostringstream filename;					//Stream
	filename << "../Results/AIRSS.csv";		//Nom du fichier
	ofstream file(filename.str().c_str()); 	//Fichier de stockage

	int Rep = 100;			//Nombre de répétition

	double Cout;	 		//Le cout de la recherche
	ArrayXd Couts(Rep);		//Couts stockées

	double Hmin, Hrange;	//Hauteur minimale & range
	double Lmin, Lrange;	//Largeur minimale & range
	double Smin, Srange;	//Sigma   minimale & range

	double epsilon = 1e-9; 	//Précision de la dérivée

	bool Doubts;			//doutes sur la présences d'autres minimas

	double Veps = 0.2;		//Vraissemblance (plus aucun minimas à trouver)

	SolNim Release;							//Système

	for(int D = 1; D <= 100; D++)			//Pour les dimension de 1 à 100
	{
		VectorXd p(D);			//On initialise le vecteur position
		MatrixXd Min(D,0);		//On initialise la Matrice des positions trouvée
		VectorXi T(0);			//On initialise le Vecteur des minimas trouvée

		for(int N = 1; N <= 16; N*=2)		//Pour 1 à 16 minimas
		{
			for(int i = 0; i < Rep; i++)	//On répète jusqu'à Rep
			{
				GenPays(D, N, Lmin, Lrange, Smin, Srange, Hmin, Hrange);
				/*On génère le paysage*/
				Cout = 0;			//Le cout initial et nul
				Release.Init();		//On initialise le paysage
				Doubts = 1; 		//On as aucun doute il existe des minimas
				while(Doubts)
				/*Tant qu'il y as des doutes sur la présences d'autres minimas*/
				{
					Cout += Steepdes(Release, Lmin, Lrange, epsilon);
					Release.get_pos(p);
					New(Min, T, p);
					Doubts = Criterion(T,Veps);
				}
				Couts(i) = Cout;
			}
		}
	}

	return 0;
}
