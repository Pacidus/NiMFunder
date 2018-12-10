#include "headers/header.hpp"

int main()
{
    cout << D << endl;
    VectorXd p(D);
    p << 0,0;
	Descent(p, 1000, 1e-5);
    cout << p << endl ;
	return 0;
}
