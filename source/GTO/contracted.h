#ifndef Contracted_H
#define Contracted_H


class Contracted_H{
	
public:
    add_primitive(double alpha, double w, int i, int j, int k, const vec pos);
    get_contracted() {return contracted}

private:
	contracted;
	normalization(double alpha, int i, int j, int k);
	fac(int n);

};
#endif // Contracted_H