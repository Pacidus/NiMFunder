/* on importe nos Lib*/
#include <iostream>   //Permet de gérer les input output (io)
#include <fstream>    //Permet de gérer des fichier (f)
#include <cmath>      //Importe quelques fonctions mathématiques
#include <cstdlib>    //Permet d'utiliser la fonction rand()
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;
/*============================================================================*/
                                /*Utilitaires*/
/*============================================================================*/

/*============================================================================*/
                                /*Fonctions*/
/*============================================================================*/
double DGauss(VectorXd& p, VectorXd& P, VectorXd& sigma, double A)
{
	/*DGauss: retourne la valeur d'une gaussienne à la position p a D dimensions
    dont le l'extremum se situe à là position P et à la hauteur A*/

    VectorXd r = (p-P).array()/sigma.array();

    double arg = r.norm();

    return A*exp(-arg*arg/2);

}//DGauss

double dDGauss(VectorXd& p, VectorXd& P, VectorXd& sigma, double A, VectorXi& I)
{
	/*dDGauss: retourne la valeur de la dérivée d'une gaussienne à la position p
    a D dimensions dont le l'extremum se situe à là position P et à la hauteur A*/
    int n = I.size();

    double Prod = 1;

    VectorXd dr = (p-P);

    VectorXd r = dr.array()/sigma.array();

    double arg = r.norm();

    dr = -2*dr;

    for(int i = 0; i < n; i++)
    {
        Prod *= dr(I(i));
    }

    return Prod*A*exp(-arg*arg/2);

}//DGauss

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

void Descent(void ddf(MatrixXd&, VectorXd&),void df(VectorXd&, VectorXd&), VectorXd& p, int i_max, double epsilon)
{
    /*Descent: aplique une méthode classique pour trouver le minima 
    (équlibre en quasi statique sum(F)=0)*/
    int D = p.size();
    VectorXd b(D);
    MatrixXd A(D,D);
    double df2,Sddfa,Sddf;
    int i = 0;
    do
    {
	i++; 				//On implémente de 1
	df(b,p);			//On génère le vécteur de la dérivée première
	df2 = b.norm();			//On génère la norme de la dérivée première
	ddf(A,p);			//On génère la matrice des dérivées secondes
	Sddf = A.sum();		//On génère la somme des derivées secondes
	Sddfa = -A.cwiseAbs().sum();	//On génère la somme de la valeur absolue des dérivées secondes
	p -= b;
	cout << (Sddfa == Sddf) << endl;
    }
    while(i < i_max and not(df2 < epsilon and Sddfa == Sddf));
    cout << i << endl;
}

//void Gradient(void ddf(Ma
