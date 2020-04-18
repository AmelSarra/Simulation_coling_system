//Headers of stationary functions

std::vector<double> calculate_a_stationary(int M,double  k,double h);
std::vector<double> calculate_b_stationary(Ailette a1);
std::vector<double> calculate_c_stationary(int M,double  k,double h);
std::vector<double> calculate_f_stationary(Ailette a1);
std::vector<double>  calculate_y(std::vector<double> f,std::vector<double> a,std::vector<double> b_star);
std::vector<double> solve_system_stationary(Ailette a1,bool verb);
std::vector<double> solve_system_exact(Ailette a1, bool verb);
void write_csv_file_stationary(double h, std::vector<double> sol1, std::vector<double> sol2);
