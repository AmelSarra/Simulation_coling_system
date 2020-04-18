class Ailette{

	public:
	//Variables static
	static const int ro = 2700,cp = 940,k = 164;

  //Constructeur
  Ailette(double lx,double ly, double lz,int m,double phiP,int te,double hc);
  Ailette(double lx,double ly, double lz,int m,double phiP,int te,double hc,double tf, double n);
	Ailette(double lx,double ly, double lz,int m,double phiP,int te,double hc,double tf, double n,double M_mx,double M_my,double M_mz);

	//Accesseurs
	double lx();
	double ly();
	double lz();
	double phiP();
	double hc();
	double mx();
	double my();
	double mz();

	int tf();
	int n();
	int m();
	int te();

	//Mutateurs
	void setLx(double lx);
	void setLy(double ly);
	void setLz(double lz);
	void setPhiP(double phiP);
	void setHc(double hc);
	void setM(int m);
	void setTe(int te);
	void setTF(double tf);
	void setN(double n);
	void setMx(double mx);
	void setMy(double my);
	void setMz(double mz);

	//Methodes
	double surface();
	double perimetre();
	double long_interval();
	double delta_t();
	double nbr_point();

	private:
	   double M_lx,M_ly,M_lz,M_phiP,M_hc,M_tf,M_n,M_mx,M_my,M_mz;
	   int M_m,M_te;

};
