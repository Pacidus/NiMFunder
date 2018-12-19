/* on importe nos Lib*/
#include "Class.hpp"

int Gabriele(VectorXd p[], double val[], int& Nbr, double L, double L0, double e, double l, double epsilon);

void map(VectorXd p[], double val[], int Nbr, double L, double L0);

void reject(VectorXd p[], double val[], int& Nbr, double e);

void toclose(VectorXd p[], double val[], int& Nbr, double l);

int release(VectorXd p[], double val[], int Nbr, double epsilon);
