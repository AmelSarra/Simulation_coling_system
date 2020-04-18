#include "Ailette.h"

	//ructeur
  Ailette::Ailette(double lx,double ly, double lz,int m,double phiP,int te,double hc):M_lx(lx),M_ly(ly),M_lz(lz),M_m(m),M_phiP(phiP),M_te(te),M_hc(hc){};
	Ailette::Ailette(double lx,double ly, double lz,int m,double phiP,int te,double hc,double tf, double n):M_lx(lx),M_ly(ly),M_lz(lz),M_m(m),M_phiP(phiP),M_te(te),M_hc(hc),M_tf(tf),M_n(n){};
  Ailette::Ailette(double lx,double ly, double lz,int m,double phiP,int te,double hc,double tf, double n,double mx,double my,double mz):M_lx(lx),M_ly(ly),M_lz(lz),M_m(m),M_phiP(phiP),M_te(te),M_hc(hc),M_tf(tf),M_n(n),M_mx(mx),M_my(my),M_mz(mz){};


	//Accesseurs
	double Ailette::lx(){return M_lx;}
	double Ailette::ly(){return M_ly;}
	double Ailette::lz(){return M_lz;}
	double Ailette::phiP(){return M_phiP;}
	double Ailette::hc(){return M_hc;}
  double Ailette::mx(){return M_mx;}
  double Ailette::my(){return M_my;}
  double Ailette::mz(){return M_mz;}

	int Ailette::tf(){return M_tf;}
	int Ailette::n(){return M_n;}
	int Ailette::m(){return M_m;}
	int Ailette::te(){return M_te;}

	//Mutateurs
	void Ailette::setLx(double lx){M_lx=lx;}
	void Ailette::setLy(double ly){M_ly=ly;}
	void Ailette::setLz(double lz){M_lz=lz;}
	void Ailette::setPhiP(double phiP){M_phiP=phiP;}
	void Ailette::setHc(double hc){M_hc=hc;}
	void Ailette::setM(int m){M_m=m;}
	void Ailette::setTe(int te){M_te=te;}
	void Ailette::setTF(double tf){M_tf=tf;}
	void Ailette::setN(double n){M_n=n;}
  void Ailette::setMx(double mx){M_mx=mx;}
  void Ailette::setMy(double my){M_mx=my;}
  void Ailette::setMz(double mz){M_mx=mz;}

	//Methodes
	double Ailette::surface(){
		return M_ly * M_lz;
	}
	double Ailette::perimetre(){
		return	2*(M_ly + M_lz);
	}
	double Ailette::long_interval(){
		return M_lx / M_m;
	}
	double Ailette::delta_t(){
		return M_tf / M_n;
	}
  double Ailette::nbr_point(){
    return M_mx*M_my*M_mz;
  }
