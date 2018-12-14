/* on importe nos Lib*/
#include "Class.hpp"



/*============================================================================*/
								/*Private func*/
/*============================================================================*/

double SolNim::alea()
{
	/*alea: retourne un double aléatoirement entre [0,1]*/
    return rand()/(RAND_MAX+1.0);
}//alea

/*============================================================================*/

double SolNim::DGauss(VectorXd& P, VectorXd& sigma, double A)
{
	/*DGauss: retourne la valeur d'une gaussienne à la position p a D dimensions
    dont le l'extremum se situe à là position P et à la hauteur A*/

    VectorXd r = (p-P).array()/sigma.array();
    r = r.array()*r.array();
    double arg = r.sum();

    return A*exp(-arg/2);

}//DGauss

/*============================================================================*/

double SolNim::dDGauss(VectorXd& P, VectorXd& sigma, double A, VectorXi& I)
{
	/*dDGauss: retourne la valeur de la dérivée d'une gaussienne à la position p
    a D dimensions dont l'extremum se situe à là position P et à la hauteur A*/
    double Prod = 1;

    VectorXd dr = (p-P);

    VectorXd r = dr.array()/sigma.array();
    r = r.array()*r.array();
    double arg = r.sum();

    dr = -dr.array()/(sigma.array()*sigma.array());

    for(int i = 0; i < D; i++)
    {
        Prod *= dr(I(i));
    }

    return Prod*A*exp(-arg/2);

}//DGauss

/*============================================================================*/
								/*Public func*/
/*============================================================================*/


/*====================================Init====================================*/
SolNim::SolNim()
{
}

void SolNim::Init()
{

	ifstream Initialiser("/home/yohan/Bureau/NiMFunder/Paysage.init");

	if(Initialiser.fail()) cout << "fichier non trouvé" << endl;

	char temp[100];


	Initialiser >> temp;

	N = atoi(temp);

	Initialiser >> temp;

	D = atoi(temp);

	Ps.resize(D,N);
	Sigma.resize(D,N);
	A.resize(D,D);
	H.resize(N);
	p.resize(D);
	b.resize(D);

	for(int i = 0; i < D; i++)
	{
	    for(int j = 0; j < N; j++)
	    {
	        Initialiser >> temp;
	        Ps(i,j) = atof(temp);
	    }
	}

	for(int i = 0; i < D; i++)
	{
	    for(int j = 0; j < N; j++)
	    {
	        Initialiser >> temp;
	        Sigma(i,j) = atof(temp);
	    }
	}

	for(int i = 0; i < N; i++)
	{
	    Initialiser >> temp;
	    H(i) = atof(temp);
	}

	return;
}

void SolNim::Rpos(double a, double b)
{
	for(int i = 0; i < D; i++)
	{
		p(i) = a*alea() + b;
	}
}

/*=====================================Set====================================*/
void SolNim::set_Num(int N1)
{
	N = N1;
}

void SolNim::set_Dim(int D1)
{
	D=D1;
}

void SolNim::set_Pos(MatrixXd Ps1)
{
	Ps = Ps1;
}

void SolNim::set_Sig(MatrixXd Sigma1)
{
	Sigma = Sigma1;
}

void SolNim::set_Hau(VectorXd H1)
{
	H = H1;
}

void SolNim::set_pos(VectorXd p1)
{
	p = p1;
}

/*=====================================Get====================================*/
void SolNim::get_Num(int& N2)
{
	N2 = N;
}

void SolNim::get_Dim(int& D2)
{
	D2 = D;
}

void SolNim::get_Pos(MatrixXd& Ps2)
{
	Ps2 = Ps;
}

void SolNim::get_Sig(MatrixXd& Sigma2)
{
	Sigma2 = Sigma;
}

void SolNim::get_A(MatrixXd& A2)
{
	A2 = A;
}

void SolNim::get_Hau(VectorXd& H2)
{
	H2 = H;
}

void SolNim::get_pos(VectorXd& p2)
{
	p2 = p;
}

void SolNim::get_b(VectorXd& b2)
{
	b2 = b;
}

/*==================================Functions=================================*/
double SolNim::landscape()
{
    /*landscape: retourne la valeur du paysage de gaussiennes en la position p*/

    double Sol = 0;
    VectorXd P(D), sigma(D);
    for(int i = 0; i < N; i++)
    {
        P = Ps.col(i);
        sigma = Sigma.col(i);
        Sol += DGauss(P, sigma, H(i));
    }
    return Sol;
}//landscape

double SolNim::dlandscape(VectorXi I)
{
    /*dlandscape: retourne la valeur de la dérivée du paysage de gaussiennes en la position p*/

    double Sol = 0;
    VectorXd P(D), sigma(D);
    for(int i = 0; i < N; i++)
    {
        P = Ps.col(i);
        sigma = Sigma.col(i);
        Sol += dDGauss(P, sigma, H(i), I);
    }
    return Sol;
}//dlandscape

void SolNim::df()
{
/*df: Génère le vecteur des dérivée premières*/
    MatrixXi I(D,D);
    I = MatrixXi::Identity(D,D);

	for(int i = 0; i < D; i++)
	{
    	b(i) = dlandscape(I.col(i));
	}
}

void SolNim::ddf()
{
/*ddf: Génère la matrice des dérivée secondes*/
	VectorXi I(D);
	I.fill(0);
	for(int i = 0; i < D ; i++)
	{
		for(int j = 0; j < D ; j++)
		{
			I(i) += 1;
			I(j) += 1;

			A(i,j) = dlandscape(I);
			I.fill(0);
		}
	}
}

