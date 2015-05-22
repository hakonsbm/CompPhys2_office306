#include <iostream>
#include <stdlib.h>
//#include <cmath>

using namespace std;

int fac(int n)
{
	int fac = 1;
	for (int i = 2; i <= n; ++i) fac = fac * i;
	return fac;
}

int main(int argc, char const *argv[])
{
	int n = fac(atoi(argv[1]));
	cout << argv[1] << "! = " << n << endl;
	return 0;
}