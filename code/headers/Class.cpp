/* on importe nos Lib*/
#include "Class.hpp"



/*============================================================================*/
								/*Private func*/
/*============================================================================*/

double SolNim::alea()
{
	/*alea: retourne un double aléatoirement entre [0,1[*/
    return rand()/(RAND_MAX+1.0); // rand() -> [0,RAND_MAX[
}//alea

/*============================================================================*/

double SolNim::DGauss(VectorXd& P, VectorXd& sigma, double A)
{
	/*DGauss: retourne la valeur d'une gaussienne à la position p a D dimensions
    dont le l'extremum se situe à là position P et à la hauteur A*/

    VectorXd r = (p-P).array()/sigma.array(); /* On gen le vecteur de la dist du
    											 point au minima que l'on divise
    											 par le vecteur des ecarts-types
    										  */

    r = r.array()*r.array(); 	/* On génère r^2 */
    double arg = r.sum();		/* On génère l'argument de l'exponentielle
    							   (gaussiènne) */

    return A*exp(-arg/2); /* On retourne la valeur de la gaussiènne */
}//DGauss

double SolNim::dDGauss(VectorXd& P, VectorXd& sigma, double A, int d)
{
	/*
	dDGauss: retourne la valeur de la dérivée d'une gaussienne à la position p a
	D dimensions dont le l'extremum se situe à là position P et à la hauteur A
	dérivée selon l'axe d
	*/

    VectorXd r = (p-P).array()/sigma.array(); /* On gen le vecteur de la dist du
    											 point au minima que l'on divise
    											 par le vecteur des ecarts-types
    										  */
    r = r.array()*r.array();				  /* On génère r^2 */
    double arg = r.sum();	/* On génère l'argument de l'exponentielle
    						   (gaussiènne) */

    return -A*(r(d)/(p-P)(d))*exp(-arg/2);
    /* On retourne la valeur de la dérivée de la gaussiènne */
}//dDGauss

/*============================================================================*/
								/*Public func*/
/*============================================================================*/


/*====================================Init====================================*/
SolNim::SolNim()
{
	/*SolNim: /!\ constructeur de la classe */
}

void SolNim::Init()
{
	/*Init: Initialise les valeurs de SolNim */
	ifstream Initialiser("../Paysage.init");
	/*On importe un fichier text dans lequel on paramètrise le système*/

	if(Initialiser.fail()) cout << "fichier non trouvé" << endl; //test

	char temp[100]; //variable temporaire


	Initialiser >> temp; //on récupaire la première chaine de caractère

	N = atoi(temp); //on la transforme en entier

	Initialiser >> temp; //on récupaire la deuxième chaine de caractère

	D = atoi(temp); //On la transforme en entier

	/*On redimensionne nos variables*/
	Ps.resize(D,N);
	Sigma.resize(D,N);
	A.resize(D,D);
	H.resize(N);
	p.resize(D);
	b.resize(D);

	/*On récupère la position des minimas*/
	for(int i = 0; i < D; i++)
	{
	    for(int j = 0; j < N; j++)
	    {
	        Initialiser >> temp;
	        Ps(i,j) = atof(temp);
	    }
	}

	/*On récupère les Sigmas*/
	for(int i = 0; i < D; i++)
	{
	    for(int j = 0; j < N; j++)
	    {
	        Initialiser >> temp;
	        Sigma(i,j) = atof(temp);
	    }
	}

	/*On récupère la profondeur des minimas*/
	for(int i = 0; i < N; i++)
	{
	    Initialiser >> temp;
	    H(i) = atof(temp);
	}

	return;
}

void SolNim::Rpos(double a, double b)
{
	/*Rpos: tire une position aléatoire dans un hypercube
	taille a et débutant à b*/
	for(int i = 0; i < D; i++)
	{
		p(i) = a*alea() + b;
	}
}

/*=====================================Set====================================*/
void SolNim::set_Num(int N1)
{
	/*set_Num: définis le nombre de minimas*/
	N = N1;
}

void SolNim::set_Dim(int D1)
{
	/*set_Dim: définis la dimention de l'espace*/
	D=D1;
}

void SolNim::set_Pos(MatrixXd Ps1)
{
	/*set_Pos: définis la matrice des positions des minimas*/
	Ps = Ps1;
}

void SolNim::set_Sig(MatrixXd Sigma1)
{
	/*set_Sig: définis la matrice des sigmas*/
	Sigma = Sigma1;
}

void SolNim::set_Hau(VectorXd H1)
{
	/*set_Hau: définis le vecteur des "Hauteurs" des minimas*/
	H = H1;
}

void SolNim::set_pos(VectorXd p1)
{
	/*set_pos: définis "notre" position*/
	p = p1;
}

/*=====================================Get====================================*/
void SolNim::get_Num(int& N2)
{
	/*get_Num: permet d'obtenir le nombre de minimas*/
	N2 = N;
}

void SolNim::get_Dim(int& D2)
{
	/*get_Dim: permet d'obtenir la dimention de l'espace*/
	D2 = D;
}

void SolNim::get_Pos(MatrixXd& Ps2)
{
	/*get_Pos: permet d'obtenir la matrice des positions des minimas*/
	Ps2 = Ps;
}

void SolNim::get_Sig(MatrixXd& Sigma2)
{
	/*get_Sig: permet d'obtenir la matrice des sigmas*/
	Sigma2 = Sigma;
}

void SolNim::get_A(MatrixXd& A2)
{
	/*get_A: permet d'obtenir la matrices des dérivées secondes*/
	A2 = A;
}

void SolNim::get_Hau(VectorXd& H2)
{
	/*get_Hau: permet d'obtenir le vecteur des "Hauteurs" des minimas*/
	H2 = H;
}

void SolNim::get_pos(VectorXd& p2)
{
	/*get_pos: permet d'obtenir "notre" position*/
	p2 = p;
}

void SolNim::get_b(VectorXd& b2)
{
	/*get_b: permet d'obtenir le vecteur des dérivées premières*/
	df();
	b2 = b;
}

/*==================================Functions=================================*/
double SolNim::landscape()
{
    /*landscape: retourne la valeur du paysage de gaussiennes en la position p*/

    double Sol = 0; //Variable de somme
    VectorXd P(D), sigma(D); //Vecteur position & Sigmas des minimas
    for(int i = 0; i < N; i++)//On parcours les minimas
    {
        P = Ps.col(i); //On récupère la position du minima
        sigma = Sigma.col(i); //On récupère les sigmas du minima
        Sol += DGauss(P, sigma, H(i));
        /*On calcule la valeur de la gaussiènne au point p et on la somme*/
    }
    return Sol;
}//landscape

double SolNim::dlandscape(int d)
{
    /*dlandscape: retourne la valeur de la dérivée du paysage de gaussiennes
	en la position p et selon l'axe d*/

    double Sol = 0; //Variable de somme
    VectorXd P(D), sigma(D);//Vecteur position & Sigmas des minimas
    for(int i = 0; i < N; i++)//On parcours les minimas
    {
        P = Ps.col(i);//On récupère la position du minima
        sigma = Sigma.col(i);//On récupère les sigmas du minima
        Sol += dDGauss(P, sigma, H(i), d);
        /*On calcule la valeur de la dérivée de la gaussiènne au point p et
        selon la dimension d puis on la somme*/
    }
    return Sol;
}//dlandscape

void SolNim::df()
{
	/*df: Génère le vecteur des dérivée premières*/
	for(int i = 0; i < D; i++)//On parcoure les dimensions
	{
    	b(i) = dlandscape(i);//On dérive au point p et celon la dimension i
    }
}//df

 /*==================================Méthodes=================================*/

void SolNim::SteepDescent(double dt)
{
    /*SteepDescent: On ce déplace d'un pas de longueur dt selon la dérivée*/

	df();			//On génère le vécteur de la dérivée première
	p -= b.normalized()*dt; //On normalise le vecteur et on ce déplace de dt
}//SteepDescent
