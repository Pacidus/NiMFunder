/* on importe nos Lib*/
#include "header.hpp"



VectorXd p(2), A(3);
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

    A << -2, -4, -4;

    p << 0,2;
    Descent(ddf, df, p, 10000, 1e-10);
    cout << p << endl ;
    return 0;
}
