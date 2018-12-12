/* on importe nos Lib*/
#include "headers/header.hpp"

VectorXd A(3);
VectorXd p(2);
MatrixXd P(2,3), sigma(2,3);

void df(VectorXd& b, VectorXd& p)
{
    VectorXi Ix(1);
    Ix << 0;

    VectorXi Iy(1);
    Iy << 1;

    b << dlandscape(p, P, sigma, A, Ix), dlandscape(p, P, sigma, A, Iy);
}

void ddf(MatrixXd& Ai, VectorXd& p)
{
    VectorXi Ixx(2);
    Ixx << 0,0;

    VectorXi Iyy(2);
    Iyy << 1,1;

    VectorXi Ixy(2);
    Ixy << 0,1;

    Ai << dlandscape(p, P, sigma, A, Ixx), dlandscape(p, P, sigma, A, Ixy),
          dlandscape(p, P, sigma, A, Ixy), dlandscape(p, P, sigma, A, Iyy);
}

int main()
{


    P << 1, .5, 1.5,
         1, .5,   1;

    sigma << .1, .3, .2,
             .1, .3, .5;

    A << -3, -4, -4;


    ofstream fichier("find.res");
    for(int i = 0; i < 10000 ; i++)
    {
        p << 2*alea(),2*alea();
        Descent(ddf, df, p, 200, 1e-5);
        fichier << p.transpose() << endl ;

    }
    fichier.close();
    return 0;
}
