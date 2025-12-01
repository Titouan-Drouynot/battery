#include "battery.h"
#include <limits>

Battery::Battery(double Vmin,double Vmax, double ac, double ad, double Rid, double nominal_capacity,double b, int series, int para, const std::string discharging_filename)
  : Vmin(Vmin), Vmax(Vmax), ac(ac), ad(ad), Rid(Rid), Ric(Rid), b(b) {
  this->M = M_init(discharging_filename, nominal_capacity, Rid, this->A1, this->I1, this->A2, this->I2, series, para);
  this->I = 0;
  this->U = this->Vmax;
}  

double interpolateLinear(double I_query, const std::vector<double>& I_values, const std::vector<double>& a_values) {
  // Cas particuliers
  if (I_query <= I_values.front()) return a_values.front();
  if (I_query >= I_values.back()) return a_values.back();

  // Trouver les deux points encadrants
  for (size_t i = 0; i < I_values.size() - 1; i++) {
    if (I_query >= I_values[i] && I_query <= I_values[i+1]) {
      // Interpolation linéaire
      double t = (I_query - I_values[i]) / (I_values[i+1] - I_values[i]);
      return a_values[i] + t * (a_values[i+1] - a_values[i]);
    }
  }
  return 0.0; // Ne devrait pas arriver ici
}

double Battery::a1(double I){
  return interpolateLinear(I,this->I1,this->A1);
}

double Battery::a2(double I){
  return interpolateLinear(I,this->I2,this->A2);
}

bool Battery::canupdate(double P, double dt){
  std::cout << "recherche des solutions\n";
  std::vector<double> inter = this->intersect(P,dt);
  std::cout << inter.size() << " solutions trouvées\n";
  double bestDiff = std::numeric_limits<double>::infinity();
  if( inter.size()==0){
    std::cout << "puissance non autorisée\n";
    return false;
  }
  return true;
}

void Battery::update(double P, double dt){
  std::cout << "recherche des solutions\n";
  std::vector<double> inter = this->intersect(P,dt);
  std::cout << inter.size() << " solutions trouvées\n";
  double bestDiff = std::numeric_limits<double>::infinity();
  if( inter.size()==0){
    std::cout << "puissance non autorisée\n";
    return;
  }
  int bestIndex=-1;
  for (int i = 0; i < (int)inter.size(); i++) {
    double diff = std::abs(inter[i] - this->U);
    if (diff < bestDiff) {
      bestDiff = diff;
      bestIndex = i;
    }
  }
  this->U = inter[bestIndex];
  this->I = P/this->U;
  if (P>=0){
    this->b+= (1 - this->Ric*this->I/this->U)*P*dt;
  }
  else{
    this->b+= (1 - this->Rid*this->I/this->U)*P*dt;
  }
}

std::vector<double> Battery::intersect(double P,double dt){
   int range =500;
   double U = this->Vmin;
   double C;
   double b;
   double Pprec;
   bool indef = false;
   double dU = (this->Vmax - this->Vmin)/range;
   std::vector<double> intersections;
   for (int i = 0; i <= range; i++){
     U = U+dU;
     C = P/U;
     if(C<=this->ac && C>=this->ad){
       b=this->b;
       if (P>=0){
         b+= (1 - this->Ric*C/U)*P*dt;
       }
       else{
         b+= (1 - this->Rid*C/U)*P*dt;
       }
       if(this->a1(C)<=b && this->a2(C)>=b){
         if(indef){ 
           double Pcand = C*interpolate(C,b,this->M);
           if(((P-Pcand) >=0)!= ((P-Pprec)>=0)){
             double t = (P - Pcand)/(Pprec - Pcand);
             intersections.push_back(U-dU*t);
             Pprec = Pcand;
           }
         }
         else{
           indef = true;
           Pprec = C*interpolate(C,b,this->M);
         }
       }
     }
   }
   return intersections;
 }

void Battery::visualize(double P,double dt,
                std::vector<double>& Iv,std::vector<double>& Pv, std::vector<double>& bv,
                std::vector<double>& Uv, std::vector<double>& Mv, bool domaine){
  int range =100;
  double U = this->Vmin;
  double dU = (this->Vmax - this->Vmin)/range;
  double C = P/U;
  double b = this->b;
  U = U - dU;
  for (int i = 0; i <= range; i++){
    U = U+dU;
    C = P/U;
    if(domaine||(C<=this->ac && C>=this->ad)){
      b=this->b;
      if (P>=0){
        b+= (1 - this->Ric*C/U)*P*dt;
      }
      else{
        b+= (1 - this->Rid*C/U)*P*dt;
      }
      if(domaine||(this->a1(C)<=b && this->a2(C)>=b)){
        Iv.push_back(C);
        Pv.push_back(C*interpolate(C,b,this->M));
        Uv.push_back(U);
        bv.push_back(b);
        Mv.push_back(interpolate(C,b,this->M));
      }
    }
  }
}
