#include <Rcpp.h>
#include <iostream>     // std::cout
#include <cmath>        // std::abs
#include <numeric>
using namespace Rcpp;
// [[Rcpp::export]]
List effort_sim(NumericVector Bini, NumericVector K, NumericVector r, NumericMatrix qflt, int Time, NumericVector Effort,
		  NumericVector lambda, NumericMatrix Pprey, NumericVector Xp, NumericVector mu, NumericVector Eff, int np,
		  NumericVector Ispredation){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  //
  // ~               Creating vectors and matrices            ~  //
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  //
  int Nsp  = K.size();
  int Nflt = Effort.size();
  NumericMatrix predation(Time, Nsp);
  NumericMatrix Biomass(Time, Nsp);
  NumericMatrix Report(2, Nsp);
  NumericMatrix prey(Time, Nsp);
  NumericMatrix Production(Time, Nsp);
  NumericMatrix Tot_Catch(Time, Nsp);
  // double Catch[Time][Nsp][Nflt] = {};
  int t = 0;
  double change = 100.0;
  for (int i = 0; i < Nsp; i++) {
    Biomass(0, i) = Bini(i);
  }

  while(change > 1){
    // Gettin predation
    for(int sp = 0; sp < Nsp; sp++){
      for(int pr = 0; pr < Nsp; pr++){
	// Calculating biomass removed by predation and competition
	predation(t, sp) += (Pprey(sp, pr) * Xp(pr) * Biomass(t, pr)) *
	  (pow(Biomass(t, sp) / K(sp), np) / (pow(mu(sp), np) + pow((Biomass(t, sp) / K(sp)), np)));
	// Calculating biomass added by predation
	prey(t, pr) += Eff(pr) * (Pprey(sp, pr) * Xp(pr) * Biomass(t, pr)) *
	  (pow(Biomass(t, sp) / K(sp), np) / (pow(mu(sp), np) + pow((Biomass(t, sp) / K(sp)), np)));
      }
    }


    for(int sp = 0; sp < Nsp; sp++){
      for(int flt = 0; flt < Nflt; flt++){
	// Catch(t, sp, flt) =  qflt(flt, sp) * Effort(flt) * Biomass(t, sp);
	Tot_Catch(t, sp) +=  qflt(flt, sp) * Effort(flt) * Biomass(t, sp);
      }
      // Checking that values of predation and prey are positive
      if(predation(t, sp) < 0){
	predation(t, sp) = 0;
      }
      if(prey(t, sp) < 0){
	prey(t, sp) = 0;
      }
      // Biomass calculation
      Biomass(t + 1, sp) = Biomass(t, sp) + (r(sp) / lambda(sp)) * Biomass(t, sp) * (1 - pow((Biomass(t, sp) / K(sp)), lambda(sp))) - predation(t, sp) + prey(t, sp) - Tot_Catch(t, sp);
      if(Biomass(t + 1, sp) <= 0){
	Biomass(t + 1, sp) = 1;  // this is to avoid extinsion
      }
      Production(t + 1, sp)  =   (r(sp) / lambda(sp)) * Biomass(t + 1, sp) * (1 - pow((Biomass(t + 1, sp) / K(sp)), lambda(sp))) - predation(t, sp) + prey(t, sp);
      if(Production(t + 1, sp) < 0){
	Production(t + 1, sp) = 0;
      }
    }
    //  checking the change statement to stop the run when necesary
    if(t > 100){
      double m1 = 0.0;
      double m2 = 0.0;
      for(int sp = 0; sp < Nsp; sp++){
	m1 += 0.5 * (Biomass(t - 1, sp) + Biomass(t - 2, sp));
	m2 += 0.5 * (Biomass(t, sp) + Biomass(t + 1, sp));
	Report(0, sp) =  0.5 * (Biomass(t, sp) + Biomass(t + 1, sp));
	Report(1, sp) =  0.5 * (Production(t, sp) + Production(t + 1, sp));
      }
      change  =  abs(m1/Nsp  - m2/Nsp);
    }
    if(t == (Time-1)) {
      change = 1;
    }
    t++;
    }

  return List::create(Report,
		      Production, Biomass, predation, prey, t);

}
