//Headers for 3D

int find_xi00(double xioo,double h, int m);
std::vector<double> model_3d(Ailette a1, std::vector<double> x,int stationary);
void write_jtk_file_3d(Ailette a1, std::vector<double> t,int stationary);
