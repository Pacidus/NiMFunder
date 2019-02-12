/*On importe nos Lib*/
#include "Class.hpp" //On importe la classe SolNim

int Gabriele(VectorXd p[], double val[], int& Nbr, double L, double L0, double e, double l, double epsilon);
/*
Gabriele: Implémentation de la méthode de Gabriele
*/

void map(VectorXd p[], double val[], int Nbr, double L, double L0);
/*
map: 1ere étape de la méthode de Gabriele "mapper" l'espace aléatoirement
*/

void reject(VectorXd p[], double val[], int& Nbr, double e);
/*
reject: 2ème étape de la méthode de Gabriele rejeter tous les points qui sont
trop "haut" par rapport au point le plus bas
*/

void toclose(VectorXd p[], double val[], int& Nbr, double l);
/*
toclose: 3ème étape de la méthode de Gabriele rejeter tous les points qui sont
trop "proches" les un des autres
*/

int release(VectorXd p[], double val[], int Nbr, double epsilon);
/*
release: 4ème étape de la méthode de Gabriele faire une "release" des points
restants
*/

void GenPays(int D, int N, double L, double L0, double sigma, double sigma0, double H, double H0);
/*
GenPays: Permet de créer un paysage de manière procédurale
*/

int Steepdes(SolNim Release, double L, double L0, double epsilon);
/*
Steepdes: Permet de faire une "release" (pour la méthode AIRSS)
*/
