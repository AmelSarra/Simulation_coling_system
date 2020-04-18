#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>
#include <time.h>

#include "classes/Ailette.h"
#include "headers/headers_stationary.h"
#include "headers/headers_instationary.h"
#include "headers/headers_3d.h"

using namespace std;

//Local headers
std::vector<double> readCfg(char **args);

int main(int argc, char **argv){

  clock_t tStart = clock();

	cout << "\n\n";
	cout << "/***************************/\n";
	cout << "/*                         */\n";
	cout << "/* Nom : Amel Sarra Madani */\n";
	cout << "/*   Projet C++            */\n";
	cout << "/*   Class : M1 CSMI       */\n";
	cout << "/*                         */\n";
	cout << "/***************************/\n";


	cout << "\n\n";
	cout << "Ce programme a pour but d'étudier et de calculer des paramètres sur une Ailette du refroidisseur d'un processeur.\n";
	cout << "...\n";
	cout << "...\n";

	cout << "Lecture des paramètres du fichier simu.cfg...\n\n";
	std::vector<double> params;
	params = readCfg(argv);

	cout << "Initialisation des variables pour les calculs...\n\n";
	std::vector<double> x,exact,t,dest;
	char c;
	bool b=false;
	double ** t_constFlux;
	double ** t_inconstFlux;

	Ailette a1(params[0],params[1],params[2],params[3],params[4],params[6],params[5],params[8],params[9],params[10],params[11],params[12]);
	cout << "Voulez vous afficher les résultats dans la console ? (o/n)\n";
	cin >> c;

	if (c == 'O' || c == 'o')
		b=true;

	// stationary if params == 0 else instationary
	if(params[7] == 0){
		cout << "Calcul dans le modèle stationnaire...\n\n";

		x = solve_system_stationary(a1,b);
		exact = solve_system_exact(a1,b);

		cout << "Ecriture dans le fichier csv...\n\n";
		write_csv_file_stationary(a1.long_interval(),x,exact);

		cout << "Calcul des paramètres pour la modélisation 3D...\n";
		t = model_3d(a1,x,0);

		cout << "Création du fichier vtk...\n";
		write_jtk_file_3d(a1,t,0);
	}
	else{
		cout << "Calcul dans le modèle instationnaire...\n\n";
		t_constFlux = solve_system_instationary_fluxConst(a1,b);
		t_inconstFlux = solve_system_instationary_fluxInconst(a1,b);

		cout << "Ecriture dans le fichier csv...\n\n";
		write_csv_file_instationary_constFlux(a1.m(),a1.long_interval(),t_constFlux);
		cout << "Ecriture dans les fichiers vtk ...\n\n";

		for(int i=1;i<=6;i++){
			cout << "Ecriture dans le fichier vtk " << i << " ...\n\n";
			for(int j=0;j<a1.m();j++)
				dest.push_back(t_constFlux[i-1][j]);
			t = model_3d(a1,dest,i);
			write_jtk_file_3d(a1,t,i);
		}

		write_csv_file_instationary_inconstFlux(a1.n(),a1.delta_t(),t_inconstFlux);
	}

	cout << "\n\n";
	cout << "Fin du programme ...\n";
	printf("Temps d'éxécution : %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	return 0;
}

std::vector<double> readCfg(char **args){

	if(!args[1]){
		cout << "Erreur, veuillez entrez un fichier de configuration.\n";
		exit(-1);
	}

	ifstream file (args[1]);
	std::string line;
	std::vector<double> tab;
	double a;
	std::string word;

	if(file)
	{
		for(int i=0;i<9;i++){
			if(i == 0 || i == 8){
				for(int j=0;j<3;j++){
					file >> word;
					file >> a;
					tab.push_back(a);
				}
			}
			else{
				file >> word;
				file >> a;
				tab.push_back(a);
			}
		getline(file,line);
		}
	}
	file.close();

	return tab;
}
