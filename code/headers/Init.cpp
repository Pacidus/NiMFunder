#include "Init.hpp"

void Init(int& N, int& D, MatrixXd& Ps, MatrixXd& Sigma, VectorXd& H)
{
	ifstream Initialiser("/home/yohan/Bureau/NiMFunder/Paysage.init");

	if(Initialiser.fail()) cout << "fichier non trouvé" << endl;

	char temp[100];


	Initialiser >> temp;

	N = atoi(temp);

	Initialiser >> temp;

	D = atoi(temp);

	Ps.resize(D,N);
	Sigma.resize(D,N);
	H.resize(N);


	for(int i = 0; i < D; i++)
	{
	    for(int j = 0; j < N; j++)
	    {
	        Initialiser >> temp;
	        Ps(i,j) = atof(temp);
	    }
	}

	for(int i = 0; i < D; i++)
	{
	    for(int j = 0; j < N; j++)
	    {
	        Initialiser >> temp;
	        Sigma(i,j) = atof(temp);
	    }
	}

	for(int i = 0; i < N; i++)
	{
	    Initialiser >> temp;
	    H(i) = atof(temp);
	}

	return;
}
