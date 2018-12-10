/* on importe nos Lib*/
#include <iostream>   //Permet de gérer les input output (io)
#include <fstream>    //Permet de gérer des fichier (f)
#include <cmath>      //Importe quelques fonctions mathématiques
#include <cstdlib>    //Permet d'utiliser la fonction rand()
#include "Eigen/Dense"//Permet d'utiliser des matrices

using namespace Eigen;
using namespace std;

/*============================================================================*/
                                /*Utilitaires*/
/*============================================================================*/

double alea();
/*alea: retourne un double aléatoirement entre [0,1]*/


/*============================================================================*/
                                /*Fonctions*/
/*============================================================================*/

double DGauss(VectorXd& p, VectorXd& P, VectorXd& sigma, double A);
/*DGauss: retourne la valeur d'une gaussienne à la position p a D dimensions
dont le l'extremum se situe à là position P et à la hauteur A*/

double dDGauss(VectorXd& p, VectorXd& P, VectorXd& sigma, double A, VectorXi& I);
/*dDGauss: retourne la valeur de la dérivée d'une gaussienne à la position p
a D dimensions dont le l'extremum se situe à là position P et à la hauteur A*/

double dlandscape(VectorXd& p, MatrixXd& P, MatrixXd& sigma, VectorXd& A,VectorXi& I);
/*dlandscape: retourne la valeur de la dérivée du paysage de gaussiennes en la position p*/

/*============================================================================*/
                                /*Méthodes*/
/*============================================================================*/ 

void Descent(VectorXd& p, int i_max, double epsilon);
/*Descent: applique une méthode classique pour trouver le minima (équilibre en quasi statique sum(F)=0)*/


