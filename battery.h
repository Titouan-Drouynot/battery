#include "setup.h"
#include <stdbool.h>

class Battery{
  private:
    std::vector<double> intersect(double P, double dt);
    //paramètre intrinsèques
    double Vmin;
    double Vmax;
    double ac;
    double ad;
    double Rid;
    double Ric;
    double a1(double I);
    double a2(double I);
    // pour calcul de a1 et a2 précalcul en certains I
    std::vector<double> I1;
    std::vector<double> I2;
    std::vector<double> A1;
    std::vector<double> A2;
  public:
    Battery(double vmin,double vmax, double ac, double ad, double rid, double nominal_capacity,double charge, int series, int para, const std::string discharging_filename);
    void visualize(double P,double dt,std::vector<double>& Iv,std::vector<double>& Pv, std::vector<double>& bv, std::vector<double>& V, std::vector<double>& M, bool domaine);
    void update(double P, double dt);
    bool canupdate(double P, double dt);
    //variable d'états
    double b;
    double U;
    double I;
    surface M;
};
