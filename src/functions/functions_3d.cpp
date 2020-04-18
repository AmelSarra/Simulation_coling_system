#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>

#include "../classes/Ailette.h"
#include "../headers/headers_stationary.h"

using namespace std;

int find_xi00(double xioo,double h, int m){
  int i;
  double x;
  for(i=0; i<=m; i++){
    x = i * h;
    if(xioo<x)
      return i;
    }
}

std::vector<double> model_3d(Ailette a1, std::vector<double> x, int stationary){

  double a,b,h;
  std::vector<double> tf, xioo;
  int k;

  if(stationary == 0)
    h = a1.long_interval();
  else
    h = a1.delta_t();

  for(int i=0; i<=a1.mx(); i++){
    xioo.push_back((i*a1.lx())/a1.mx());
    k = find_xi00(xioo[i],h,a1.m());
    a = (x[k] - x[k-1])/((h * k) - (h * (k - 1)));
    b = x[k-1] - (h*(k-1)*a);
    tf.push_back((a*xioo[i])+b);
  }
  return tf;
}

void write_jtk_file_3d(Ailette a1, std::vector<double> t, int stationary){
  std::ofstream file;
  string str = to_string(stationary-1);
  if(stationary == 0)
	 file.open("results/vtk/stationary/stationary_parametres_3d.vtk");
   else
   file.open("results/vtk/instationary/instationary_parametres_3d."+str+".vtk");

  file << "# vtk DataFile Version 2.0\n";
  file << "vtk output\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_GRID\n";
  file << "DIMENSIONS " << a1.mx() << " " << a1.my() << " " << a1.mz() << "\n";
  file << "POINTS " << a1.nbr_point() << " double\n";
  for(int z = 0; z < a1.mz(); z++){
    for(int y = 0; y < a1.my(); y++){
      for(int x = 0; x<a1.mx(); x++)
        file << x << ' ' << y << ' ' << z << "\n";
    }
  }

  file << "POINT_DATA "<< a1.nbr_point() << "\n";
  file << "FIELD FieldData 1\n";
  file << "Temperature " << a1.nbr_point() << " double\n";

  for(int i=0; i<a1.mz();i++){
    for(int j= 0;j<a1.my();j++){
      for(int z=0;z<a1.mx();z++){
        file << t[z] << '\n';
      }}}

  file.close();
}
