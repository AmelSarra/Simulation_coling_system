//Headers of instationary functions

std::vector<double> calculate_a_instationary(int M,double  k,double h);
std::vector<double> calculate_b_instationary(Ailette a1);
std::vector<double> calculate_c_instationary(int M,double  k,double h);
std::vector<double>  calculate_y(std::vector<double> f,std::vector<double> a,std::vector<double> b_star);
double ** solve_system_instationary_fluxConst(Ailette a1,bool verb);
double ** solve_system_instationary_fluxInconst(Ailette a1,bool verb);
void write_csv_file_instationary_constFlux(int m,double h,double **t);
void write_csv_file_instationary_inconstFlux(int n,double delta_t,double **t);
