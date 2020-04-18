#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>

#include "../classes/Ailette.h"
#include "../headers/headers_instationary.h"

using namespace std;


std::vector<double> calculate_a_instationary(int M,double  k,double h){

		std::vector<double> a;

		for(int i =0;i<M-2;i++)
			a.push_back(-k/(h*h));
		a.push_back(k/h);

	return a;
}

std::vector<double> calculate_b_instationary(Ailette a1){

	std::vector<double> b;

	b.push_back(a1.k/a1.long_interval());

	for(int i=1;i<a1.m()-1;i++){
		b.push_back(((2*a1.k)/(a1.long_interval()*a1.long_interval()))+((a1.hc()*a1.perimetre())/a1.surface())+((a1.ro*a1.cp)/a1.delta_t()));
	}
	b.push_back(-a1.k/a1.long_interval());

	return b;
}

std::vector<double> calculate_c_instationary(int M,double  k,double h){
		std::vector<double> c;
		c.push_back(-k/h);
		for(int i =1;i<M-1;i++)
			c.push_back(-k/(h*h));

	return c;
}

std::vector<double>  calculate_y(std::vector<double> f,std::vector<double> a,std::vector<double> b_star){

	std::vector<double> y;

	y.push_back(f[0]/b_star[0]);

	for(int i=1;i<b_star.size();i++){
		y.push_back((f[i]-a[i-1]*y[i-1])/b_star[i]);
	}
	return y;

}

double ** solve_system_instationary_fluxConst(Ailette a1, bool verb){

	std::vector<double> a = calculate_a_instationary(a1.m(),a1.k,a1.long_interval());
	std::vector<double> b = calculate_b_instationary(a1);
	std::vector<double> c = calculate_c_instationary(a1.m(),a1.k,a1.long_interval());

	std::vector<double> c_star;
	std::vector<double> b_star;


	std::vector<double> f;
	std::vector<double> y;

	int tmp[6] = {15,30,60,90,150,210},z=0;
	float tmp_courrant = 1;

	double x[a1.m()];
	double** t;
	t = new double*[6];

	b_star.push_back(b[0]);
	c_star.push_back(c[0]/b_star[0]);

	for(int i=1;i<a1.m()-1;i++){
		b_star.push_back(b[i]-a[i-1]*c_star[i-1]);
		c_star.push_back(c[i]/b_star[i]);
	}
	b_star.push_back(b[a1.m()-1]-a[a1.m()-2]*c_star[a1.m()-2]);

	for(int i=0; i<a1.m();i++)
		x[i] = a1.te();


	for(int i=0;i<a1.n();i++){
		tmp_courrant= i*a1.delta_t();

		//Calcule de F de t
		f.push_back(a1.phiP());
		for(int j=1;j<a1.m()-1;j++){
			f.push_back((a1.hc()*a1.perimetre()*a1.te()/a1.surface())+((a1.ro*a1.cp*x[j])/a1.delta_t()));
		}
		f.push_back(0);

		y = calculate_y(f,a,b_star);

		x[a1.m()-1] = y[a1.m()-1];
		for(int j=a1.m()-2;j>=0;j--){
			x[j] = y[j]-c_star[j]*x[j+1];
		}
		f.clear();
		y.clear();

		if(tmp_courrant>=tmp[z]){
			t[z] = new double[a1.m()];
			for(int j=0;j<a1.m();j++)
			   t[z][j] = x[j];
				z++;
			if(z>=6)
				break;
		}
	}

  if(!verb)
    return t;

	cout << "Vecteur a pour la méthode instationnaire: \n";
	for(int i;i<a.size();i++)
		cout << a[i] << " ";
	cout << "\n\n";

	cout << "Vecteur b pour la méthode instationnaire: \n";
	for(int i;i<b.size();i++)
		cout << b[i] << " ";
	cout << "\n\n";

	cout << "Vecteur c pour la méthode instationnaire: \n";
	for(int i;i<c.size();i++)
		cout << c[i] << " ";
	cout << "\n\n";

	cout << "Vecteur b_star pour la méthode instationnaire: \n";
	for(int i;i<b_star.size();i++)
		cout << b_star[i] << " ";
	cout << "\n\n";

	cout << "Vecteur c_star pour la méthode instationnaire: \n";
	for(int i;i<c_star.size();i++)
		cout << c_star[i] << " ";
	cout << "\n\n";

	cout << "Vecteur t solution pour la méthode instationnaire: \n";
	for(int i=0;i<6;i++){
		for(int j=0;j<a1.m();j++)
			cout << t[i][j] << " ";
		cout << "\n\n";
	}
	cout << "\n\n";

	return t;
}

double ** solve_system_instationary_fluxInconst(Ailette a1, bool verb){

	std::vector<double> a = calculate_a_instationary(a1.m(),a1.k,a1.long_interval());
	std::vector<double> b = calculate_b_instationary(a1);
	std::vector<double> c = calculate_c_instationary(a1.m(),a1.k,a1.long_interval());

	std::vector<double> c_star;
	std::vector<double> b_star;


	std::vector<double> f;
	std::vector<double> y;

	int tmp[10] = {30,60,90,120,150,180,210,240,270,300},z=0;
	float tmp_courrant = 1;

	double x[a1.m()];
	double** t;
	t = new double*[3];
	t[0] = new double[a1.n()];
	t[1] = new double[a1.n()];
	t[2] = new double[a1.n()];

	b_star.push_back(b[0]);
	c_star.push_back(c[0]/b_star[0]);

	for(int i=1;i<a1.m()-1;i++){
		b_star.push_back(b[i]-a[i-1]*c_star[i-1]);
		c_star.push_back(c[i]/b_star[i]);
	}
	b_star.push_back(b[a1.m()-1]-a[a1.m()-2]*c_star[a1.m()-2]);

	for(int i=0; i<a1.m();i++)
		x[i] = a1.te();


	for(int i=0;i<a1.n();i++){
		tmp_courrant= i*a1.delta_t();

		//Calcule de F de t
		f.push_back(a1.phiP());
		for(int j=1;j<a1.m()-1;j++){
			f.push_back((a1.hc()*a1.perimetre()*a1.te()/a1.surface())+((a1.ro*a1.cp*x[j])/a1.delta_t()));
		}
		f.push_back(0);

		y = calculate_y(f,a,b_star);

		x[a1.m()-1] = y[a1.m()-1];
		for(int j=a1.m()-2;j>=0;j--){
			x[j] = y[j]-c_star[j]*x[j+1];
		}
		f.clear();
		y.clear();

		if(tmp_courrant>=tmp[z]){
			if (a1.phiP() == 0)
				a1.setPhiP(125000);
			else
				a1.setPhiP(0);
				z++;
			if(z>=10)
				break;
		}

		t[0][0] = 20;
		t[0][i+1] = x[0];
		t[1][i] = x[(a1.m()/2)-1];
		t[2][i] = x[a1.m()-1];
}

  if (!verb)
    return t;

	cout << "Vecteur a pour la méthode instationnaire: \n";
	for(int i;i<a.size();i++)
		cout << a[i] << " ";
	cout << "\n\n";

	cout << "Vecteur b pour la méthode instationnaire: \n";
	for(int i;i<b.size();i++)
		cout << b[i] << " ";
	cout << "\n\n";

	cout << "Vecteur c pour la méthode instationnaire: \n";
	for(int i;i<c.size();i++)
		cout << c[i] << " ";
	cout << "\n\n";

	cout << "Vecteur b_star pour la méthode instationnaire: \n";
	for(int i;i<b_star.size();i++)
		cout << b_star[i] << " ";
	cout << "\n\n";

	cout << "Vecteur c_star pour la méthode instationnaire: \n";
	for(int i;i<c_star.size();i++)
		cout << c_star[i] << " ";
	cout << "\n\n";


	cout << "Vecteur t solution pour la méthode instationnaire: \n";
	for(int i=0;i<3;i++){
		for(int j=0;j<a1.n();j++)
			cout << t[i][j] << " ";
		cout << "\n\n";
	}
	cout << "\n\n";
	return t;

}

void write_csv_file_instationary_constFlux(int m,double h,double **t){

	std::ofstream file;
	file.open("results/csv/result_instationary_constFlux.csv");
	double x = 0;
	file << "x,t=15s,t=30s,t=60s,t=90s,t=150s,t=210s\n";
	for(int i =0; i<m; i++){
		file << x << ',' << t[0][i] << ',' << t[1][i] << ',' << t[2][i] << ',' << t[3][i] << ',' << t[4][i] << ',' << t[5][i] << '\n';
		x = x+h;
	}
}

void write_csv_file_instationary_inconstFlux(int n,double delta_t,double **t){

	std::ofstream file;
	file.open("results/csv/result_instationary_inconstFlux.csv");
	double x = 0;
	file << "t,x(0),x(M/2),x(M)\n";
	for(int i =0; i<n; i++){
		file << x << ',' << t[0][i] << ',' << t[1][i] << ',' << t[2][i] << '\n';
		x = x+delta_t;
	}
}
