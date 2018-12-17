/* on importe nos Lib*/
#include <iostream>   //Permet de gérer les input output (io)
#include <fstream>    //Permet de gérer des fichier (f)
#include <cmath>      //Importe quelques fonctions mathématiques
#include <cstdlib>    //Permet d'utiliser la fonction rand()
#include "Eigen/Dense"//Permet d'utiliser des matrices

using namespace Eigen;
using namespace std;

class SolNim
{
	private:
		/*Variable*/

		int N;
		int D;
		MatrixXd Ps;
		MatrixXd Sigma;
		MatrixXd A;
		VectorXd H;
		VectorXd p;
		VectorXd b;

		/*Functions*/

		double alea();
		/*alea: retourne un double aléatoirement entre [0,1]*/

		double DGauss(VectorXd& P, VectorXd& sigma, double A);
		/*DGauss: retourne la valeur d'une gaussienne à la position p a D dimensions
    	dont le l'extremum se situe à là position P et à la hauteur A*/

		double dDGauss(VectorXd& P, VectorXd& sigma, double A, VectorXi& I);
		/*dDGauss: retourne la valeur de la dérivée d'une gaussienne à la position p
	    a D dimensions dont le l'extremum se situe à là position P et à la hauteur A*/

	public:
		/*Init*/
		SolNim();
		void Init();
		void Rpos(double a, double b);

		/*Set*/
		void set_Num(int N1);
		void set_Dim(int D1);
		void set_Pos(MatrixXd Ps1);
		void set_Sig(MatrixXd Sigma1);
		void set_Hau(VectorXd H1);
		void set_pos(VectorXd p1);

		/*Get*/
		void get_Num(int& N2);
		void get_Dim(int& D2);
		void get_Pos(MatrixXd& Ps2);
		void get_Sig(MatrixXd& Sigma2);
		void get_A(MatrixXd& A2);
		void get_Hau(VectorXd& H2);
		void get_pos(VectorXd& p2);
		void get_b(VectorXd& b2);

		/*Functions*/
		double landscape();
		/*landscape: retourne la valeur du paysage de gaussiennes en la position p*/
		double dlandscape(VectorXi I);
    	/*dlandscape: retourne la valeur de la dérivée du paysage de gaussiennes en la position p*/
		void df();
		/*df: Génère le vecteur des dérivée premières*/
		void ddf();
		/*ddf: Génère la matrice des dérivée secondes*/

		/*Methodes*/
		void SteepDescent();
};
