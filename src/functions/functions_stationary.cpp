#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>

#include "../classes/Ailette.h"
#include "../headers/headers_stationary.h"

using namespace std;


std::vector<double> calculate_a_stationary(int M,double  k,double h){

		std::vector<double> a;

		for(int i =0;i<M-2;i++)
			a.push_back(-k/(h*h));
		a.push_back(-k/h);

	return a;
}

std::vector<double> calculate_b_stationary(Ailette a1){

	std::vector<double> b;

	b.push_back(a1.k/a1.long_interval());

	for(int i=1;i<a1.m()-1;i++){
		b.push_back(((2*a1.k)/(a1.long_interval()*a1.long_interval()))+((a1.hc()*a1.perimetre())/a1.surface()));
	}
	b.push_back(a1.k/a1.long_interval());

	return b;
}

std::vector<double> calculate_c_stationary(int M,double  k,double h){
		std::vector<double> c;
		c.push_back(-k/h);
		for(int i =1;i<M-1;i++)
			c.push_back(-k/(h*h));

	return c;
}

std::vector<double> calculate_f_stationary(Ailette a1){

	std::vector<double> f;

	f.push_back(a1.phiP());
	for(int i=1;i<a1.m()-1;i++){
		f.push_back(a1.hc()*a1.perimetre()*a1.te()/a1.surface());
	}
	f.push_back(0);

	return f;
}

std::vector<double> solve_system_stationary(Ailette a1, bool verb){


	std::vector<double> a = calculate_a_stationary(a1.m(),a1.k,a1.long_interval());
	std::vector<double> b = calculate_b_stationary(a1);
	std::vector<double> c = calculate_c_stationary(a1.m(),a1.k,a1.long_interval());

	std::vector<double> c_star;
	std::vector<double> b_star;

	std::vector<double> f = calculate_f_stationary(a1);
	std::vector<double> y(a1.m());

	std::vector<double> x(a1.m());
	std::vector<double>::iterator en=x.end();
	std::vector<double>::iterator it=x.begin();

	int k = a1.m()-2;

	b_star.push_back(b[0]);
	c_star.push_back(c[0]/b_star[0]);

	for(int i=1;i<a1.m()-1;i++){
		b_star.push_back(b[i]-a[i-1]*c_star[i-1]);
		c_star.push_back(c[i]/b_star[i]);
	}
	b_star.push_back(b[a1.m()-1]-a[a1.m()-2]*c_star[a1.m()-2]);

	y = calculate_y(f,a,b_star);

	en--;
	*en = y[a1.m()-1];
	en--;

	for(;en!=it;en--){
		*en = y[k]-c_star[k]*x[k+1];
		k--;
	}
	*en = y[0]-c_star[0]*x[1];

	if(!verb)
		return x;
	cout << "Vecteur a pour la méthode stationnaire: \n";
	for(int i;i<a.size();i++)
		cout << a[i] << " ";
	cout << "\n\n";

	cout << "Vecteur b pour la méthode stationnaire: \n";
	for(int i;i<b.size();i++)
		cout << b[i] << " ";
	cout << "\n\n";

	cout << "Vecteur c pour la méthode stationnaire:\n";
	for(int i;i<c.size();i++)
		cout << c[i] << " ";
	cout << "\n\n";

	cout << "Vecteur b_star pour la méthode stationnaire:\n";
	for(int i;i<b_star.size();i++)
		cout << b_star[i] << " ";
	cout << "\n\n";

	cout << "Vecteur c_star pour la méthode stationnaire: \n";
	for(int i;i<c_star.size();i++)
		cout << c_star[i] << " ";
	cout << "\n\n";

	cout << "Vecteur f pour la méthode stationnaire:\n";
	for(int i;i<f.size();i++)
		cout << f[i] << " ";
	cout << "\n\n";

	cout << "Vecteur y pour la méthode stationnaire: \n";
	for(int i;i<y.size();i++)
		cout << y[i] << " ";
	cout << "\n\n";

	cout << "Vecteur X solution pour la méthode stationnaire:\n";
	for(int i;i<x.size();i++)
		cout << x[i] << " ";
	cout << "\n\n";

	return x;

}

std::vector<double> solve_system_exact(Ailette a1 ,bool verb){
	double 	a=(a1.hc()*a1.perimetre())/(a1.k*a1.surface());
	std::vector<double> t;
	double x=0,y,z;
	for(int i = 0; i<a1.m();i++){
		y = a1.phiP()*cosh(sqrt(a)*a1.lx())*cosh(sqrt(a)*(a1.lx()-x));
		z = a1.k * sqrt(a) * sinh(sqrt(a)*a1.lx()) * cosh(sqrt(a)*a1.lx());
		t.push_back(a1.te() + (y/z));
		x = x+a1.long_interval();
	}
	if (!verb)
		return t;
	cout << "Vecteur solution exact : \n";
	for(int i=0;i<a1.m();i++)
		cout << t[i] << " ";
	cout << "\n";

	return t;
}

void write_csv_file_stationary(double h, std::vector<double> sol1, std::vector<double> sol2){

	std::ofstream file;
	file.open("results/csv/result_stationary.csv");
	double x=0;
	file << "x,numeric,exact\n";
	for(int i=0;i<sol1.size();i++){
		file << x << ',' << sol1[i] << ',' << sol2[i] << '\n';
		x = x+h;
	}

	file.close();
}
