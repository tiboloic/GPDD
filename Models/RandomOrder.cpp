// LT 26/08/2019

// TMB Hierarchical model: adds random effect on taxonomic order


#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  // DATA
  
  // observations
  DATA_MATRIX(obs);
  // length of each time-series 
  DATA_IVECTOR(lengths);
  // indicator of each observation: is NA or not
  DATA_IMATRIX(nas);
  // order each time-series belongs to
  DATA_FACTOR(orders);
  
  //META parameters
  // average b
  PARAMETER(b);

  // Time-series level parameters
  PARAMETER_VECTOR(mu);
  PARAMETER_VECTOR(lsig);
  PARAMETER_VECTOR(ltau);

  //log-stdev of random effect on b between orders
  PARAMETER(lnu0);
  
  // log-stdev of random effect on b within order. smae variance for all orders, as supported by AIC
  PARAMETER(lnu1);

  // random effects on b between order
  PARAMETER_VECTOR (bdevs0);

  // random effects on b within order
  PARAMETER_VECTOR (bdevs1);
  
  // REPORT
  //REPORT(exp(2*lnu0));
  //REPORT(exp(2*lnu1));
  //REPORT(bdevs0);
  
  
  Type nLL=0.0;

  Type m, v, cste;

  Type thisb;
  Type thisa;

  // add likelihood contribution for order level random effects
  nLL -= dnorm(bdevs0, Type(0), Type(1), true).sum();

  // add likelihood contribution for time-series level random effects
  nLL -= dnorm(bdevs1, Type(0), Type(1), true).sum();
    
  // likelihood contribution per time-series
  for(int its=0;its<lengths.size();its++){

    // order time-series belongs to
    int iorder = orders(its);
    
    // gaussian random effect with mean b + order specific ranef + time-series specific ranef
    thisb = b + exp(lnu0) * bdevs0(iorder) + exp(lnu1) * bdevs1(its);
    
    thisa = mu(its)*(1-thisb);
      
    // initial conditions
    m = obs(its,1);
    v = 10 + exp(ltau(its));

    // contribution to likelihood
    nLL-= dnorm(obs(its,0), m, sqrt(v), true);

    for (int i=1; i<lengths(its); i++) {

          // is this observation missing ?
          if (nas(its,i-1)==1) {
              // observation is missing, just project ahead
              m = thisa + thisb * m;
              v = thisb*thisb*(v-exp(ltau(its))) + exp(lsig(its)) + exp(ltau(its));
          }
          else
          {
              cste = (v-exp(ltau(its)))/v;
              m =  thisa + thisb * (m + cste*(obs(its,i-1)-m));
	      v = thisb*thisb * cste * exp(ltau(its)) + exp(lsig(its)) + exp(ltau(its));
          }
          if (nas(its,i)==0)
              // contribution to the likelihood
 	          nLL-= dnorm(obs(its,i),m,sqrt(v),true);
      }
  }
  return nLL;
}

