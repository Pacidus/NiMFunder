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

