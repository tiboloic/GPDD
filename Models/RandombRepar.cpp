// LT 14/10/2014


#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(obs);
  DATA_IVECTOR(lengths);
  DATA_IMATRIX(nas);

  PARAMETER(b);
  PARAMETER_VECTOR(mu);
  PARAMETER_VECTOR(lsig);
  PARAMETER_VECTOR(ltau);

  //log-stdev of random effect on b
  PARAMETER(lnu)

  // random effects on b
  PARAMETER_VECTOR (bdevs)


  Type nLL=0;

  Type m, v, cste;

  Type thisb;
  Type thisa;

  for(int its=0;its<lengths.size();its++){

     // contribution of random effect to likelihood
    nLL -= dnorm(bdevs(its), Type(0), Type(1),true);
    
    // gaussian random effect with mean b and stdev exp(lnu)
    thisb = b + exp(lnu) * bdevs(its);

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

