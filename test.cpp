/* on importe nos Lib*/
#include "header.hpp"

int main()
{

    VectorXd p(2);
    MatrixXd P(2,3), sigma(2,3);

    P << 1, .5, 1.5,
         1, .5,   1;

    sigma << .1, .3, .2,
             .1, .3, .5;

    VectorXd A(3);
    A << -2, -4, -4;

    ofstream fichier("map.val");
    
    for(int i = 0; i < 1000; i++)
    {

        p(0) = i/500.;
   
        for(int j = 0; j <1000; j++)
        {
            p(1) = j/500.;
            fichier << landscape(p, P, sigma, A) <<" ";
        }
    
        fichier << endl;    

    }
    
    fichier.close();
    return 0;
}
