/* on importe nos Lib*/
#include "headers/Class.hpp"

void FindStep(int& NRe, int& NRa, int N_max)
{
	SolNim Release, Randomap;
	Release.Init();
	Randomap.Init();

	int D;
<<<<<<< HEAD
	Release.get_Dim(D);
	VectorXd b(D);

=======
	PayF.get_Dim(D);
	MatrixXd Ps(N,D);
	PayF.get_Pos(Ps);
	MatrixXd Sigma(N,D);
	PayF.get_Sig(Ps);
	VectorXd H(N);
	PayF.get_Hau(H);
	VectorXd p(D);
	VectorXd b(D);
>>>>>>> db2fd827ec57b64c87ee009ed955c86dcfd696cb
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

<<<<<<< HEAD
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
=======
	int n;
	cin >> n;
	double epsilon = 1e-5;
	int i = 0;
	double dt = 1;
	for(int j = 0; j < n; j++)
	{
		PayF.Rpos(1.5,.25);
		i = 1;

   		do
    	{
    	dt = 0.1/i;
		i++; 				//On implémente de 1
		PayF.SteepDescent(dt);
		PayF.get_b(b);
		PayF.get_pos(p);
		file << p.transpose() << " " << PayF.landscape() << endl;
		}
	    while(b.norm()*b.norm() > epsilon*epsilon and i < 1000);

>>>>>>> db2fd827ec57b64c87ee009ed955c86dcfd696cb
	}
}



int main()
{
	int NRe = 0;
	int NRa = 0;
	FindStep(NRe, NRa, 10);
    return 0;
}
