/*On importe nos Lib*/
#include "Methods.hpp"

/*=============================Méthode de Gabriele============================*/

int Gabriele(VectorXd p[], double val[], int& Nbr, double L, double L0, double e, double l, double epsilon)
{
	/*Gabriele: Implémentation de la méthode de Gabriele*/
	map(p, val, Nbr, L, L0); //1ere etape
	reject(p, val, Nbr, e); //2eme etape
	toclose(p, val, Nbr, l); //3eme etape
	int couts; //couts de la recherche de minimas
	couts = release(p, val, Nbr, epsilon); //4eme etape
	return couts;
}//Gabriele

void map(VectorXd p[], double val[], int Nbr, double L, double L0)
{
	/*map: 1ere étape de la méthode de Gabriele "mapper" l'espace
	aléatoirement*/
	SolNim Map;		//Map et un objet SolNim
	Map.Init();		//On initialise SolNim avec Le terrain
	int D;			//On initialise D (le nombre de dimentions de l'espace)
	Map.get_Dim(D);	//On obtiens la valeur de D
	VectorXd p1(D);	//On génère un vecteur position
	for(int i = 0; i < Nbr; i++)
	{/*On parcourt les positions qui servirons à mapper l'espace*/
		Map.Rpos(L,L0);		//On génère une position au hasard
		Map.get_pos(p1);	//On récupère la position
		p[i] = p1;			//On stocke la position dans la liste des positions
		val[i] = Map.landscape();
		//On récupère la valeur de la fonction au point p
	}
}//map

void reject(VectorXd p[], double val[], int& Nbr, double e)
{
	/*reject: 2ème étape de la méthode de Gabriele rejeter tous les points qui
	sont trop "haut" par rapport au point le plus bas*/
	double min = val[0];
	//On prend comme valeur la plus faible la première valeur
	int imin = 0; //Index de la valeur la plus faible
	VectorXd palt[Nbr]; //positions apres le tri
	double valt[Nbr]; //valeures apres le tri
	int N = 0; //Valeur de N (le Nombre de positions après tri)

	for(int i = 1; i < Nbr; i++) //Pour toutes les positions
	{
		if(min > val[i])//si la valeur en cette position et inferieur au minima
		{
			imin = i; 		//l'index du minima change
			min = val[i]; 	//la valeur du minima change
		}
	}

	for(int i = 0; i < Nbr; i++)	//pour toutes les positions
	{
		if(val[i] <= min+e)
		//si la valeur en cette position n'est pas trop loin du minima
		{
			palt[N] = p[i];		//On garde ce point
			valt[N] = val[i];	//On garde cette valeur
			N++;				//On incrémente l'indice
		}
	}

	Nbr = N;				//La nouvelle valeur de Nbr
	VectorXd palt2[Nbr];	//Les nouvelles positions
	double valt2[Nbr];		//Les nouvelles valeurs

	for(int i = 0; i < Nbr; i++) // On réaloue les valeurs
	{
		palt2[i] = palt[i];
		valt2[i] = valt[i];
	}
	p = palt2; 		//On change la liste des positions
	val = valt2;	//On change la liste des valeurs
}//reject

void toclose(VectorXd p[], double val[], int& Nbr, double l)
{
	/*toclose: 3ème étape de la méthode de Gabriele rejeter tous les points qui
	sont trop "proches" les un des autres*/
	VectorXd palt[Nbr];		//On crée une nouvelle liste de vecteur position
	double valt[Nbr];		//On crée une nouvelle liste de valeurs
	int N = 0;				//On initie le nombre de valeur dans les listes à 0
	bool k = 1;				//Booléen il est le plus faible de ces voisins
	for(int i = 0; i < Nbr; i++) //Pour tous les éléments de la liste
	{
		k = 1;					 //On considère vrai jusqu'à preuve du contraire
		for(int j = 0; j < Nbr; j ++)	//Pour tous les éléments de la liste
		{
			if(val[i] > val[j]) //si la valeur de i est sup a j
			{
				if((p[i]-p[j]).norm() < l) //Et que leurs distances est faible
				{
					k = 0;		//Il n'est pas un point valide
					break;		//On sort de la boucle
				}
			}
		}
		if(k)	//Si le point verifie les conditions alors
		{
			palt[N] = p[i]; 	//On conserve sa position
			valt[N] = val[i];	//On conserve sa valeur
			N++;				//On as donc un élément en plus
		}
	}

	Nbr = N;						//On met à jour la valeur de Nbr
	VectorXd palt2[Nbr];			//On gen les valeurs de p
	double valt2[Nbr];				//On gen les valeurs de val

	for(int i = 0; i < Nbr; i++)	//On met à jour les valeurs de p & val
	{
		palt2[i] = palt[i];			//On met à jour les valeurs de p
		valt2[i] = valt[i];			//On met à jour les valeurs de val
	}
	p = palt2;						//On met à jour p
	val = valt2;					//On met à jour val
}//toclose

int release(VectorXd p[], double val[], int Nbr, double epsilon)
{
	/*release: 4ème étape de la méthode de Gabriele faire une "release" des points
	restants*/
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
			dt = .2/j;
			Release.SteepDescent(dt);
			Release.get_b(b);
			j++;
		}
		while(b.norm() > epsilon and j > 1001);

		Ntot += (j-1);
	}

	return Ntot;
}//release

void GenPays(int D, int N, double L, double L0, double sigma, double sigma0, double H, double H0)
{
	/*GenPays: Permet de créer un paysage de manière procédurale*/
	ostringstream filename;
	filename << "../Paysage.init";
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
}//GenPays

int Steepdes(SolNim Release, double L, double L0, double epsilon)
{
	/*Steepdes: Permet de faire une "release" (pour la méthode AIRSS)*/
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
}//Steepdes
