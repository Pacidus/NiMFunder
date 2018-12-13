#include "header.hpp"

/*============================================================================*/
                                /*Utilitaires*/
/*============================================================================*/
double alea()
{
    /*alea: retourne un double aléatoirement entre [0,1]*/
    return rand()/(RAND_MAX+1.0);
}//alea

/*============================================================================*/
                                /*Fonctions*/
/*============================================================================*/

double DGauss(VectorXd& p, VectorXd& P, VectorXd& sigma, double A)
{
	/*DGauss: retourne la valeur d'une gaussienne à la position p a D dimensions
    dont le l'extremum se situe à là position P et à la hauteur A*/

    VectorXd r = (p-P).array()/sigma.array();
    r = r.array()*r.array();
    double arg = r.sum();

    return A*exp(-arg/2);

}//DGauss

/*============================================================================*/

double dDGauss(VectorXd& p, VectorXd& P, VectorXd& sigma, double A, VectorXi& I)
{
	/*dDGauss: retourne la valeur de la dérivée d'une gaussienne à la position p
    a D dimensions dont le l'extremum se situe à là position P et à la hauteur A*/
    int n = I.size();

    double Prod = 1;

    VectorXd dr = (p-P);

    VectorXd r = dr.array()/sigma.array();
    r = r.array()*r.array();
    double arg = r.sum();

    dr = -dr.array()/(sigma.array()*sigma.array());

    for(int i = 0; i < n; i++)
    {
        Prod *= dr(I(i));
    }

    return Prod*A*exp(-arg/2);

}//DGauss

/*============================================================================*/


double landscape(VectorXd& p, MatrixXd& P, MatrixXd& sigma, VectorXd& A)
{
    /*landscape: retourne la valeur du paysage de gaussiennes en la position p*/

    int n = A.size();
    double Sol = 0;
    VectorXd Pi(n), sigmai(n);
    for(int i = 0; i < n; i++)
    {
        Pi = P.col(i);
        sigmai = sigma.col(i);
        Sol += DGauss(p, Pi, sigmai, A(i));
    }
    return Sol;
}//landscape

/*============================================================================*/

double dlandscape(VectorXd& p, MatrixXd& P, MatrixXd& sigma, VectorXd& A,VectorXi& I)
{
    /*dlandscape: retourne la valeur de la dérivée du paysage de gaussiennes en la position p*/

    int n = A.size();
    double Sol = 0;
    VectorXd Pi(n), sigmai(n);
    for(int i = 0; i < n; i++)
    {
        Pi = P.col(i);
        sigmai = sigma.col(i);
        Sol += dDGauss(p, Pi, sigmai, A(i), I);
    }
    return Sol;
}//dlandscape

/*============================================================================*/
                                /*Méthodes*/
/*============================================================================*/

int Descent(void ddf(MatrixXd&, VectorXd&),void df(VectorXd&, VectorXd&), VectorXd& p, int i_max, double epsilon)
{
    /*Descent: aplique une méthode classique pour trouver le minima
    (équlibre en quasi statique sum(F)=0)*/
    int D = p.size();
    VectorXd b(D);
    double df2;
    int i = 0;
    do
    {
	i++; 				//On implémente de 1
	df(b,p);			//On génère le vécteur de la dérivée première
	df2 = b.norm();			//On génère la norme de la dérivée première
	p -= b.normalized()/pow(i,1.5);
    }
    while(i < i_max and not(df2 < epsilon));
	return i;
}

void SteepDescent(void ddf(MatrixXd&, VectorXd&),void df(VectorXd&, VectorXd&), VectorXd& p, int i_max, double epsilon)
{
    /*SteepDescent: applique une méthode classique pour trouver le minima (en approximent la fonction à l'ordre 2)*/
    int D = p.size();
    VectorXd b(D);
    VectorXd r(D);
    MatrixXd A(D,D);

    df(b,p);
	ddf(A,p);

    double alpha;
    VectorXd q(D);

    int i = 0;
    r = b - A*p;

    double delta = r.transpose()*r;
    double delta0 = delta*epsilon*epsilon;

    while(i < i_max and  delta > delta0)
    {
    q = A*r;
    alpha = delta/(r.transpose()*q);
    p = p + alpha*r;
    if(i%1 == 0)
    {
		df(b,p);
    	ddf(A,p);
        r = b - A*p;
    }
	else r -= alpha*q;
	delta = r.transpose()*r;
	i++; 				//On implémente de 1
    }
}

