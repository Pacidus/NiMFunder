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

double DGauss(VectorXd& p, VectorXd& P, VectorXd& sigma, double a)
{
	/*DGauss: retourne la valeur d'une gaussienne à la position p a D dimensions
    dont le l'extremum se situe à là position P et à la hauteur a*/

    VectorXd r = (p-P).array()/sigma.array();

    double arg = r.norm();

    return a*exp(-arg*arg/2);

}//DGauss

/*============================================================================*/

double dDGauss(VectorXd& p, VectorXd& P, VectorXd& sigma, double A, VectorXi& I)
{
	/*dDGauss: retourne la valeur de la dérivée d'une gaussienne à la position p
    a D dimensions dont le l'extremum se situe à là position P et à la hauteur A*/
    double Prod = 1;

    VectorXd dr = (p-P);

    VectorXd r = dr.array()/sigma.array();

    double arg = r.norm();

    dr = -2*dr;

    for(int i = 0; i < N; i++)
    {
        Prod *= dr(I(i));
    }

    return Prod*A*exp(-arg*arg/2);

}//DGauss

/*============================================================================*/

double dlandscape(VectorXd& p, VectorXi& I)
{
    /*dlandscape: retourne la valeur de la dérivée du paysage de gaussiennes en la position p*/

    double Sol = 0;
    VectorXd Pi(N), sigmai(N);
    for(int i = 0; i < N; i++)
    {
        Pi = Pm.col(i);
        sigmai = Sigma.col(i);
        Sol += dDGauss(p, Pi, sigmai, H[i], I);
    }
    return Sol;
}//dlandscape

/*============================================================================*/
                                /*Méthodes*/
/*============================================================================*/ 

void Descent(VectorXd& p, int i_max, double epsilon)
{
    /*Descent: aplique une méthode classique pour trouver le minima 
    (équlibre en quasi statique sum(F)=0)*/
    VectorXd b(D);
    MatrixXd A(D,D);
    double df2,Sddfa,Sddf;
    int i = 0;
    do
    {
	i++; 				            //On implémente de 1
	dPaysage(p,b);			        //On génère le vécteur de la dérivée première
	df2 = b.norm();			        //On génère la norme de la dérivée première
	d2Paysage(p,A);			        //On génère la matrice des dérivées secondes
	Sddf = A.sum();		            //On génère la somme des derivées secondes
	Sddfa = -A.cwiseAbs().sum();	//On génère la somme de la valeur absolue des dérivées secondes
	p -= .1*b;
    }
    while(i < i_max and not(df2 < epsilon and Sddfa == Sddf));
}
