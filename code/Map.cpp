/* on importe nos Lib*/
#include "headers/Class.hpp"

int main()
{
	SolNim PayF;
	PayF.Init();

	ostringstream filename;

	filename << "/home/yohan/Bureau/NiMFunder/Results/Min" << N << "Dim" << D << "Map.res";

	ofstream file(filename.str().c_str());

	int n = 1000;
	double a, c;
	a = 2./n;
	c = 0;

	for(int j = 0; j <= n; j++)
	{
		p(0) = a*j+c;
		for(int i = 0; i <= n; i++)
		{
			p(1) = a*i+c;
			PayF.set_pos(p);
			file << p.transpose() << " " << PayF.landscape() << endl;
		}
		file << endl;
	}
    return 0;
}
