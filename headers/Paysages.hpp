#include "header.hpp"

/*============================================================================*/
                                /*Fonctions*/
/*============================================================================*/

double Paysage(VectorXd& p);
/*Paysage: prend en entrée un vecteur position et retourne la valeur de la fonction en ce point*/

void dPaysage(VectorXd& p, VectorXd& b);
/*dPaysage: prend en entrée un vecteur position et un vecteur des dérivées et le remplis*/

void d2Paysage(VectorXd& p);
/*d2Paysage: prend en entrée un vecteur position et une Matrice des dérivées secondes et le remplis*/
