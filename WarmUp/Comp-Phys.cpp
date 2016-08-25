// Comp-Phys.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

double deriv_Forward(double x,double h, double f(double x));
double deriv_Three_point(double x, double h, double f(double x));
vector<double> deriv_Three_point(double x, vector<double> &h, double f(double x));
vector<double> error(vector<double> &computed, double exact);
void writeToFile(vector<double> &err, vector<double> &h, string name);
void writeVecToFile(ofstream & outFile, vector<double> & vec);



int main()
{

	vector<double> h_vec;
	for (int j = 0; j < 10; j++)
	{
		h_vec.push_back(pow(10,-j));
	}
	
	ofstream outFile("output/errors.bin", ios::out | ios::binary);

	double exact = 1. / 3;
	double x = sqrt(2.);
	double h = 0.0001;
	//cout << deriv_Forward(x,h,atan) << endl;
	//cout << deriv_Three_point(x, h, atan) << endl;
	vector<double> computed = deriv_Three_point(x, h_vec, atan);

	
	vector<double> err = error(computed, exact);

	for (int i = 0; i < computed.size(); i++) {
		cout << "For h= " << h_vec[i] << " We got error: " << computed[i] << "  And error: "<< err[i] << endl;
	}
	
	writeVecToFile(outFile, err);
	outFile.close();
    return 0;
}

double deriv_Forward(double x,double h, double f(double x)) {

	return (f(x+h) - f(x)) / h;
}

double deriv_Three_point(double x, double h, double f(double x)) {
	return (f(x+h) - f(x-h)) / (2*h);
}

vector<double> deriv_Three_point(double x, vector<double> &h, double f(double x)) {
	vector<double> ret;
	for (int i = 0; i < h.size();i++)
	{
		ret.push_back(deriv_Three_point(x, h[i], f));
	}
	return ret;
}


vector<double> error(vector<double> &computed, double exact) {
	vector<double> ret;
	for (int i = 0; i < computed.size(); i++) {
		ret.push_back(log10(abs((computed[i] - exact) / exact)));
	}
	return ret;
}

void writeVecToFile(ofstream & outFile, vector<double> & vec)
{
	
	outFile.write(reinterpret_cast<char*>(&vec[0]), vec.size()*sizeof(double));
}



