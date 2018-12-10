#include "Paysages.hpp"

/*============================================================================*/
                                /*Fonctions*/
/*============================================================================*/

double Paysage(VectorXd& p)
{
    /*Paysage: prend en entrée un vecteur position et retourne la valeur de la fonction en ce point*/
    double Sol = 0;
    VectorXd Pi(N), sigmai(N);
    for(int i = 0; i < N; i++)
    {
        Pi = Pm.col(i);
        sigmai = Sigma.col(i);
        Sol += DGauss(p, Pi, sigmai, A(i));
    }
    return Sol;
}

/*============================================================================*
