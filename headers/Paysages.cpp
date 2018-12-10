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

/*============================================================================*/

void dPaysage(VectorXd& p, VectorXd& b)
{
    /*dPaysage: prend en entrée un vecteur position et un vecteur des dérivées et le remplis*/
    VectorXd I(D);
    for(int i = 0; i < D; i++)
    {
        I.fill(0);
        I(i) += 1;
        b(i) =  dlandscape(p, I)
    }
}

/*============================================================================*/

void d2Paysage(VectorXd& p, VectorXd& b)
{
    /*d2Paysage: prend en entrée un vecteur position et une Matrice des dérivées secondes et le remplis*/
    VectorXd I(D);
    for(int i = 0; i < D; i++)
    {
        for(int j = 0; j < D; j++)
        {
            I.fill(0);
            I(i) += 1;            
            I(j) += 1;
            b(i) =  dlandscape(p, I);
        }
    }
}
