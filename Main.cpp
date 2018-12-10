#include "headers/header.hpp"

/*============================================================================*/
                                /*Variables*/
/*============================================================================*/

ifstream init;              /*On récupère D dans le fichier*/
int N,D;
double H;      	           	    /*H la "hauteur" des minimas H < 0*/
MatrixXd Sigma;	            /*Sigma la matrices des écarts types */
MatrixXd Pm;		            /*Pm la matrice des positions des minimas*/

int main()
{
    init.open("../Paysage.init");   /*Fichier d'initialisation*/
    init >> N;                          /*On récupère N dans le fichier*/
    init >> D;
    int H[N];

    for(int i = 0; i < N; i++)
    {
        init >> H[i];
        for(int j = 0; j < D; j++)
        {
            init >> Sigma(j,i);
            init >> Pm(j,i);
        }
    }
    cout << D << endl;
    VectorXd p(D);
    p << 0,0;
	Descent(p, 1000, 1e-5);
    cout << p << endl ;
	return 0;
}
