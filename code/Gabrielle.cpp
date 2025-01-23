/*On importe nos Lib*/
#include "headers/Methods.hpp"

int main()
{
	ostringstream filename;
	filename << "../Results/Gabrielle.csv";
	ofstream file(filename.str().c_str());

	return 0;
}
