/* on importe nos Lib*/
#include "headers/Class.hpp"

void FindStep(int& NRe, int& NRa, int N_max)
{
	SolNim Release, Randomap;
	Release.Init();
	Randomap.Init();

	int D;
	Release.get_Dim(D);
	VectorXd b(D);

	ostringstream filename;
	filename << "/home/yohan/Bureau/NiMFunder/Results/FindStep.res";
	ofstream file(filename.str().c_str());

	if(NRe+NRa == 0)
	{
		ofstream file(filename.str().c_str());
	}
	else
	{
		fstream file;
		file.open(filename.str().c_str(), fstream::out | fstream::app);
	}

	bool Re, Ra;
	Re = 0;
	Ra = 0;

	int j = 0;
	int step = 0;
	double dt = 1;
	double epsilon = 1e-15;

	Release.Rpos(2,0);

	while(NRa*NRe < N_max*N_max)
	{
		step++;

		j++;
		dt = 0.2/j;
		Release.SteepDescent(dt);
		Release.get_b(b);
		Re = (b.norm() < epsilon);
		if(Re)
		{
		NRe++;
		j = 0;
		Release.Rpos(2,0);
		}

		Randomap.Rpos(2,0);
		Randomap.get_b(b);
		Ra = (b.norm() < epsilon);
		if(Ra) NRa++;

		file << step << " " << NRe << " " << NRa << endl;
		Re = 0;
		Ra = 0;
	}
}



int main()
{
	int NRe = 0;
	int NRa = 0;
	FindStep(NRe, NRa, 10);
    return 0;
}
