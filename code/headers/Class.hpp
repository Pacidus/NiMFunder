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

		int N; 				//Nombre de minima
		int D;				//Nombre de Dimensions
		MatrixXd Ps;		//Matrice des Positions des minimas
		MatrixXd Sigma;		//Matrice des Minimas selon les axes
		MatrixXd A;			//Matrice des derivées secondes (inutilisées)
		VectorXd H;			//Vecteur des Hauteurs
		VectorXd p;			//Vecteur des positions
		VectorXd b;			//Vecteur des dérivées première

		/*Functions*/

		double alea();
		/*alea: retourne un double aléatoirement entre [0,1]*/

		double DGauss(VectorXd& P, VectorXd& sigma, double A);
		/*
		DGauss: retourne la valeur d'une gaussienne à la position p a D
		dimensions dont le l'extremum se situe à là position P et à la hauteur A
		*/

		double dDGauss(VectorXd& P, VectorXd& sigma, double A, int d);
		/*
		dDGauss: retourne la valeur de la dérivée d'une gaussienne à la position
		p a D dimensions dont le l'extremum se situe à là position P et à la
		hauteur A dérivée selon l'axe d
		*/

	public:
		/*Init*/
		SolNim();  						/* /!\ constructeur de la classe */
		void Init();					/* Initialise les valeurs de SolNim */
		void Rpos(double a, double b);	/* Genère une position p aléatoire */

		/*Set*/
		void set_Num(int N1);			/* Permet de choisir la valeur de N */
		void set_Dim(int D1);			/* Permet de choisir la valeur de D */
		void set_Pos(MatrixXd Ps1);		/* Permet de choisir la valeur de Pos */
		void set_Sig(MatrixXd Sigma1);	/* Permet de choisir la valeur de Sig */
		void set_Hau(VectorXd H1);		/* Permet de choisir la valeur de H*/
		void set_pos(VectorXd p1);		/* Permet de choisir la valeur de p*/

		/*Get*/
		void get_Num(int& N2);			/* Permet d'obtenir la valeur de N */
		void get_Dim(int& D2);			/* Permet d'obtenir la valeur de D */
		void get_Pos(MatrixXd& Ps2);	/* Permet d'obtenir la valeur de Pos */
		void get_Sig(MatrixXd& Sigma2);	/* Permet d'obtenir la valeur de Sig */
		void get_A(MatrixXd& A2);		/* Permet d'obtenir la valeur de A */
		void get_Hau(VectorXd& H2);		/* Permet d'obtenir la valeur de H */
		void get_pos(VectorXd& p2);		/* Permet d'obtenir la valeur de p */
		void get_b(VectorXd& b2);		/* Permet d'obtenir la valeur de b */

		/*Functions*/
		double landscape();
		/*
		landscape: retourne la valeur du paysage de gaussiennes en la position p
		*/
		double dlandscape(int d);
		/*
		dlandscape: retourne la valeur de la dérivée du paysage de gaussiennes
		en la position p et selon l'axe d
		*/
		void df();
		/*df: Génère le vecteur des dérivée premières*/

		/*Methodes*/
		void SteepDescent(double dt);
		/*SteepDescent: On ce déplace d'un pas de longueur dt selon la dérivée*/
};
