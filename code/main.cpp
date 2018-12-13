/* on importe nos Lib*/
#include "headers/Init.hpp"

int N;
int D;
MatrixXd Ps(1,1);
MatrixXd Sigma(1,1);
VectorXd H(1);

void df(VectorXd& b, VectorXd& p)
{
    VectorXi Ix(1);
    Ix << 0;

    VectorXi Iy(1);
    Iy << 1;

    b << dlandscape(p, Ps, Sigma, H, Ix), dlandscape(p, Ps, Sigma, H, Iy);
}

void ddf(MatrixXd& Ai, VectorXd& p)
{
    VectorXi Ixx(2);
    Ixx << 0,0;

    VectorXi Iyy(2);
    Iyy << 1,1;

    VectorXi Ixy(2);
    Ixy << 0,1;

    Ai << dlandscape(p, Ps, Sigma, H, Ixx), dlandscape(p, Ps, Sigma, H, Ixy),
          dlandscape(p, Ps, Sigma, H, Ixy), dlandscape(p, Ps, Sigma, H, Iyy);
}

int main()
{
	Init(N, D, Ps, Sigma, H);

	VectorXd p(D);

	ostringstream filename;

	filename << "/home/yohan/Bureau/NiMFunder/Results/Min" << N << "Dim" << D << ".res";

	ofstream file(filename.str().c_str());

	int n;
	for(int j = 0; j < 10000000; j++)
	{
		for(int i = 0; i < D; i++)
		{
			p(i) = alea()*2;
		}
		n = Descent(ddf,df,p, j, .1/j);
		file << j << " " << n << " " <<  p.transpose() << endl;
	}


    return 0;
}
