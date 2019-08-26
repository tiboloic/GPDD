// LT 14/10/2014


#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(obs);
  DATA_IVECTOR(lengths);
  DATA_IMATRIX(nas);

  PARAMETER(b);
  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(lsig);
  PARAMETER_VECTOR(ltau);

  Type nLL=0;

  Type m, v, cste;

  for(int its=0;its<lengths.size();its++){

    
      // initial conditions
      m = obs(its,1);
      v = 10 + exp(ltau(its));

      // contribution to likelihood
      nLL-= dnorm(obs(its,0), m, sqrt(v), true);

      for (int i=1; i<lengths(its); i++) {

          // is this observation missing ?
          if (nas(its,i-1)==1) {
              // observation is missing, just project ahead
              m = a(its) + b * m;
              v = b*b*(v-exp(ltau(its))) + exp(lsig(its)) + exp(ltau(its));
          }
          else
          {
              cste = (v-exp(ltau(its)))/v;
              m =  a(its) + b * (m + cste*(obs(its,i-1)-m));
	      v = b*b * cste * exp(ltau(its)) + exp(lsig(its)) + exp(ltau(its));
          }
          if (nas(its,i)==0)
              // contribution to the likelihood
 	    nLL-= dnorm(obs(its,i),m,sqrt(v),true);
      }
  }
  return nLL;
}

