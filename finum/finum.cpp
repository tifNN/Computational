#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <complex>
#include "finum.hpp"

using namespace std;

/* One-dimensional Normal Law. Density function */
double nd(double x)
{
  return(exp(-SQR(x)/2.)/(sqrt(2.*PI)));
}

/*One-Dimensional Normal Law. Cumulative distribution function. */
/*Abramowitz, Milton et Stegun, Handbook of MathematicalFunctions, 1968, Dover Publication, New York, page 932 (26.2.18).Precision 10-7*/
double N(double x)
{
  const double p= 0.2316419;
  const double b1= 0.319381530;
  const double b2= -0.356563782;
  const double b3= 1.781477937;
  const double b4= -1.821255978;
  const double b5= 1.330274429;
  const double one_over_twopi= 0.39894228;
  
  double t;
  
  if(x >= 0.0)
    {
      t = 1.0 / ( 1.0 + p * x );
      return (1.0 - one_over_twopi * exp( -x * x / 2.0 ) * t * ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    } 
  else 
    {/* x < 0 */
      t = 1.0 / ( 1.0 - p * x );
      return ( one_over_twopi * exp( -x * x / 2.0 ) * t * ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}





// Closed formulae for the BlackScholes Call price and delta.
// return OK
int finum::Call_BlackScholes_73(double s, double Maturity, double Strike, double r,double divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma)
{
  double sigmasqrt,d1,d2,delta;
  double t,k;

  t=Maturity;
  k=Strike;

  if(k>0)
    {
  sigmasqrt=sigma*sqrt(t);
  d1=(log(s/k)+(r-divid)*t)/sigmasqrt+sigmasqrt*0.5;
  d2=d1-sigmasqrt;
  delta=exp(-divid*t)*N(d1);

  /*Price*/

  *ptprice=s*delta-exp(-r*t)*k*N(d2);

	
  /*Delta*/
  *ptdelta=delta;

  /*Gamma*/
  *ptgamma = exp(-0.5*d1*d1 - divid*t)/sqrt(2.*PI)/sigmasqrt/s;
     }
  else
    {
      *ptprice = s*exp(-divid*t);
      *ptdelta = 1;
      *ptgamma = 0.;
    }
  return 0;
}

// Closed formulae for the BlackScholes Put price and delta.
// return OK
int finum::Put_BlackScholes_73(double s, double Maturity, double Strike, double r,double divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma)
{
  double sigmasqrt,d1,d2,delta;
  double t,k;

  t=Maturity;
  k=Strike;

  if(k>0)
    {
  sigmasqrt=sigma*sqrt(t);
  d1=(log(s/k)+(r-divid)*t)/sigmasqrt+sigmasqrt*0.5;
  d2=d1-sigmasqrt;
  delta=-exp(-divid*t)*N(-d1);

  /*Price*/

  *ptprice=exp(-r*t)*k*N(-d2)+delta*s;

	
  /*Delta*/
  *ptdelta=delta;

  /*Gamma*/
  *ptgamma = exp(-0.5*d1*d1 - divid*t)/sqrt(2.*PI)/sigmasqrt/s;
     }
  else
    {
      *ptprice = 0.;
      *ptdelta = -1;
      *ptgamma = 0.;
    }
  return 0;
}

/* BS implied vol*/
int finum::ImpliedVol(double price,double s0,double Maturity, double Strike,double eqr, double eqdiv,double sigma_ref,double *volimpl, double *impliedDelta, double *impliedGamma)
{

  double result1,minor,major,price_major,price_minor,eps=PRICE_ACCURACY, dummy, dummy1, K,sigma;

  minor=major=sigma_ref;

  Put_BlackScholes_73(s0, Maturity, Strike, eqr, eqdiv, sigma_ref, &result1, &dummy, &dummy1);

  K=Strike;



  // Dichotomy to find the implied volatility. If the length of the interval (minor, major) is < eps*FAC_SLOPE, the Dichotomy is breaked.
  if (result1 > price)
    {
      do {
        minor=minor*0.5;
        Put_BlackScholes_73(s0, Maturity, Strike, eqr, eqdiv, minor, &result1, &dummy, &dummy1);

        //cout << "sigma:"<<minor<<"   prix:"<<result1<<"   prixm:"<<Price<<"   intr:"<< Payoff(s0,Strike) <<"  bound:"<<s0-exp(-eqr*Maturity)*Strike <<endl;
        //getchar();

        //if (minor<0.0001) return -1;

      } while( result1 >price && fabs(result1-price)> eps);
      price_minor=result1;
    }
  else
    {
      do {
        major=major*2;
        Put_BlackScholes_73(s0, Maturity, Strike, eqr, eqdiv, major, &result1, &dummy, &dummy1);

      } while( result1 <price && fabs(result1-price)> eps);
      price_major=result1;
    }

  sigma=(minor+major)*0.5;

  Put_BlackScholes_73(s0, Maturity, Strike, eqr, eqdiv, sigma, &result1, &dummy, &dummy1);

  while((fabs(result1- price)>eps)&&(major-minor>eps*FAC_SLOPE))
    {
      if(result1>price)
        {
          major=sigma;
        }
      else
        {
          minor=sigma;
        }

      sigma=(minor+major)*0.5;
      Put_BlackScholes_73(s0, Maturity, Strike, eqr, eqdiv, sigma, &result1, &dummy, &dummy1);


    }
  *volimpl= sigma;
  *impliedDelta = dummy;
  *impliedGamma = dummy1;

  return 0;
}

static int formula(double s,double k,double r,double divid,double sigma,double t,double l, double rebate,int phi,int eta,double *A,double *B,double *C,double *D,double *E,double *F,double *dA,double *dB,double *dC,double *dD,double *dE,double *dF)
    {
    double b,x1,x2,y1,y2,z,mu,lambda,sigmasqrt;

    sigmasqrt=sigma*sqrt(t);
    b=r-divid;
    mu=(b-SQR(sigma)/2.)/SQR(sigma);
    lambda=sqrt(SQR(mu)+2.*r/SQR(sigma));
    x1=log(s/k)/sigmasqrt + (1+mu)*sigmasqrt;
    x2=log(s/l)/sigmasqrt + (1+mu)*sigmasqrt;
    y1=log(SQR(l)/(s*k))/sigmasqrt+(1+mu)*sigmasqrt;
    y2=log(l/s)/sigmasqrt + (1+mu)*sigmasqrt;
    z=log(l/s)/sigmasqrt + lambda*sigmasqrt;
    *A=phi*s*exp((b-r)*t)*N(phi*x1)-phi*k*exp(-r*t)*N(phi*x1-phi*sigmasqrt);
    *B=phi*s*exp((b-r)*t)*N(phi*x2)-phi*k*exp(-r*t)*N(phi*x2-phi*sigmasqrt);
    *C=phi*s*exp((b-r)*t)*pow(l/s,2.*(1.+mu))*N(eta*y1)-
        phi*k*exp(-r*t)*pow(l/s,2.*mu)*N(eta*y1-eta*sigmasqrt);
    *D=phi*s*exp((b-r)*t)*pow(l/s,2.*(1.+mu))*N(eta*y2)-
        phi*k*exp(-r*t)*pow(l/s,2.*mu)*N(eta*y2-eta*sigmasqrt);
    *E=rebate*exp(-r*t)*(N(eta*x2-eta*sigmasqrt)-pow(l/s,2.*mu)*N(eta*y2-eta*sigmasqrt));
    *F=rebate*(pow(l/s,mu+lambda)*N(eta*z)+pow(l/s,mu-lambda)*N(eta*z-2.*eta*lambda*sigmasqrt));

    *dA=phi*exp(-divid*t)*N(phi*x1);

    *dB=phi*exp(-divid*t)*N(phi*x2)+exp(-divid*t)*nd(x2)/(sigma*sqrt(t))*(1.-k/l);

    *dC=-phi*2.*mu*pow(l/s,2.*mu)*(1./s)*(s*exp(-divid*t)*SQR(l/s)*N(eta*y1)
        -k*exp(-r*t)*N(eta*y1-eta*sigma*sqrt(t)))-
        phi*pow(l/s,2.*mu+2.)*exp(-divid*t)*N(eta*y1);

    *dD=-2.*mu*(phi/s)*pow(l/s,2.*mu)*(s*exp(-divid*t)*SQR(l/s)*N(eta*y2)-
        k*exp(-r*t)*N(eta*(y2-sigma*sqrt(t))))-
        phi*pow(l/s,2.*mu+2.)*exp(-divid*t)*N(eta*y2)-phi*eta*exp(-divid*t)*
        pow(l/s,2.*mu+2.)*nd(y2)/(sigma*sqrt(t))*(1.-k/l);

    *dE=2.*(rebate/s)*exp(-r*t)*pow(l/s,2.*mu)*(N(eta*(y2-sigma*sqrt(t)))*mu+
        eta*nd(y2-sigma*sqrt(t))/(sigma*sqrt(t)));

    *dF=-pow(l/s,mu+lambda)*(rebate/s)*((mu+lambda)*N(eta*z)+(mu-lambda)*pow(s/l,2.*lambda)*N(eta*(z-2.*lambda*sigma*sqrt(t))))-
        2.*eta*rebate*pow(l/s,mu+lambda)*nd(z)/(s*sigma*sqrt(t));

    return 0;
    }

int finum::CallUpOut_ReinerRubinstein(double s,double k,double l,double rebate,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
    {
    int phi,eta;
    double A,B,C,D,E,F;
    double dA,dB,dC,dD,dE,dF;

    phi=1;
    eta=-1;
    formula(s,k,r,divid,sigma,t,l,rebate,phi,eta,&A,&B,&C,&D,&E,&F,
        &dA,&dB,&dC,&dD,&dE,&dF);
    if (k>=l)
        {
        *ptprice=F;
        *ptdelta=dF;
        }
    else
        {
        *ptprice=A-B+C-D+F;
        *ptdelta=dA-dB+dC-dD+dF;
        }
    return 0;
    }


int finum::Fixed_CallLookback_ConzeWiswanathan(double s, double s_max, double k, double t, double r,
              double divid, double sigma, double *ptprice, double *ptdelta)
{
  double  b,sigmasqrt,a1,a2,esp,disc;

  if (s_max < s)
    {
      *ptprice=0.;
      *ptdelta=0.;
    }
  else
    {
      b=r-divid;
      sigmasqrt=sigma*sqrt(t);
      esp=2.*b/SQR(sigma);
      disc=exp(-r*t);

      if (k>s_max)
    {
      a1=(log(s/k)+ (b+SQR(sigma)/2.)*t)/sigmasqrt;
      a2=a1-sigmasqrt;
      if (b == 0)
        {
          *ptprice = s*disc*(1.+SQR(sigma)*t/2.+log(s/k))*N(a1) +
                 s*disc*sigmasqrt*nd(a1) - k*disc*N(a2);

          *ptdelta = disc*N(a1)*(2.+SQR(sigma)*t/2.+log(s/k)) +
                disc*nd(a1)*(1.+SQR(sigma)*t)/sigmasqrt   -
                disc*(k/s)*nd(a2)/sigmasqrt;
        }
      else
        {
          *ptprice=s*exp(-divid*t)*N(a1)-k*exp(-r*t)*N(a2)+
        s*exp(-r*t)*(SQR(sigma)/(2.*b))*
        (-pow(s/k,-esp)*N(a1-(2.*b/sigma)*sqrt(t))+exp(b*t)*N(a1));

          *ptdelta=exp(-divid*t)*N(a1)*(1.+SQR(sigma)/(2.*b))+
        exp(-divid*t)*nd(a1)/(sigma*sqrt(t))-
        exp(-r*t)*(k/s)*nd(a2)/sigmasqrt+
        exp(-r*t)*pow(s/k,-esp)*N(a1-2.*(b/sigma)*sqrt(t))*(1.-SQR(sigma)/(2*b));
        }
    }
      else
    {
      a1=(log(s/s_max)+ (b+SQR(sigma)/2.)*t)/sigmasqrt;
      a2=a1-sigmasqrt;
      if (b == 0)
        {
          *ptprice = disc*(s_max-k) + s*disc*(1.+SQR(sigma)*t/2.+log(s/s_max))*N(a1) +
                 s*disc*sigmasqrt*nd(a1) - s_max*disc*N(a2) ;

          *ptdelta = disc*N(a1)*(2.+SQR(sigma)*t/2.+log(s/s_max)) +
                 disc*nd(a1)*(1.+SQR(sigma)*t)/sigmasqrt   -
                 disc*(s_max/s)*nd(a2)/sigmasqrt;
        }
      else
        {

          *ptprice=exp(-r*t)*(s_max-k)+s*exp(-divid*t)*N(a1)-
        s_max*exp(-r*t)*N(a2)+
        s*exp(-r*t)*(SQR(sigma)/(2.*b))*
        (-pow(s/s_max,-esp)*N(a1-(2.*b/sigma)*sqrt(t))+exp(b*t)*N(a1));

          *ptdelta=exp(-divid*t)*N(a1)*(1.+SQR(sigma)/(2.*b))+
        exp(-divid*t)*nd(a1)/(sigma*sqrt(t))-
        exp(-r*t)*(s_max/s)*nd(a2)/sigmasqrt+
        exp(-r*t)*pow(s/s_max,-esp)*N(a1-2.*(b/sigma)*sqrt(t))*(1.-SQR(sigma)/(2*b));
        }
    }
    }


  return 0;
}

int finum::CallDownOut_ReinerRubinstein(double s,double k,double l,double rebate,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
    {
    int phi,eta;
    double A,B,C,D,E,F;
    double dA,dB,dC,dD,dE,dF;

    phi=1;
    eta=1;
    formula(s,k,r,divid,sigma,t,l,rebate,phi,eta,&A,&B,&C,&D,&E,&F,
        &dA,&dB,&dC,&dD,&dE,&dF);
    if (k>=l)
        {
        *ptprice=A-C+F;
        *ptdelta=dA-dC+dF;
        }
    else
        {
        *ptprice=B-D+F;
        *ptdelta=dB-dD+dF;
        }

    return 0;
    }


static int CallOut_KunitomoIkeda_91(double s,double l,double  u,double  Rebate,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
    {
    double d1,d2,d3,d4,mu1,mu2,mu3,F,sum1,sum2,time,delta1,delta2;
    int n;

    sum1=0.0;
    sum2=0.0;
    time=0.;
    delta1=0.0;
    delta2=0.0;
    F=u*exp(delta1*t);
    for(n=-5;n<=5;n++)
	  {
        mu1=2.*(r-divid-delta2-(double)n*(delta1-delta2))/SQR(sigma)+1.0;
        mu2=2.*(double)n*(delta1-delta2)/SQR(sigma);
        mu3=mu1;
        d1=(log(s*pow(u,2.0*(double)n)/(k*pow(l,2.0*(double)n)))+
            (r-divid+SQR(sigma)/2.0)*t)/(sigma*sqrt(t));
        d2=(log(s*pow(u,2.0*n)/(F*pow(l,2.0*(double)n)))+
            (r-divid+SQR(sigma)/2.0)*t)/(sigma*sqrt(t));
        d3=(log(pow(l,2.0*(double)n+2.)/(k*s*pow(u,2.0*(double)n)))+
            (r-divid+SQR(sigma)/2.0)*t)/(sigma*sqrt(t));
        d4=(log(pow(l,2.0*(double)n+2.)/(F*s*pow(u,2.0*(double)n)))+
            (r-divid+SQR(sigma)/2.0)*t)/(sigma*sqrt(t));
        sum2+=pow(pow(u/l,(double)n),mu1-2.)*pow(l/s,mu2)*(N(d1-sigma*sqrt(t))-
            N(d2-sigma*sqrt(t)))-
            pow(pow(l,(double)(n+1))/(s*pow(u,(double)n)),mu3-2.)*(N(d3-sigma*sqrt(t))-N(d4-sigma*sqrt(t)));
        sum1+=pow(pow(u/l,(double)n),mu1)*pow(l/s,mu2)*(N(d1)-N(d2))-
          pow(pow(l,(double)(n+1))/(s*pow(u,(double)n)),mu3)*(N(d3)-N(d4));
     }

    /*Price*/
    *ptprice=s*exp(-divid*t)*sum1-k*exp(-r*t)*sum2;

    /*Delta*/
    *ptdelta=0.;

    return 0;
    }

int finum::Call_KunitomoIkeda_91(double s,double L,double U,double Rebate,double Strike,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta){
    int dummy;
    double out_price,out_delta,price_plus,price_minus;

    dummy=CallOut_KunitomoIkeda_91(s,L,U,Rebate,Strike,t,r,divid,sigma,&out_price,&out_delta);

    /*Price*/
    *ptprice=out_price;
	
    dummy=CallOut_KunitomoIkeda_91(s*(1.+INC),L,U,Rebate,Strike,t,r,divid,sigma,&out_price,&out_delta);
    price_plus=out_price;
    dummy=CallOut_KunitomoIkeda_91(s*(1.-INC),L,U,Rebate,Strike,t,r,divid,sigma,&out_price,&out_delta);
    price_minus=out_price;

	 /*Delta*/
    *ptdelta=(price_plus-price_minus)/(2.*s*INC);

    return 0;
    }


struct heston_parms {double K; double x; double r;  double v; double tau;  
			double lambda; double theta; 
                     double rho; double eta;  int j;};


extern "C"{

double heston_integrand_j(double u, void *p){ 
    struct heston_parms* parms = (struct heston_parms*)p;
    double K = (parms->K); double x = (parms->x); 
    double v = (parms->v); double r = (parms->r);
    double lambda = (parms->lambda);
    double theta = (parms->theta);
    double rho = (parms->rho);
    double eta = (parms->eta);
    double tau = (parms->tau);
    double j = (parms->j);
    double eta_sqr = pow(eta,2);
    double uj;    double bj;
    if (j==1){	uj=0.5;  bj=lambda-rho*eta; }
    else {      uj=-0.5; bj=lambda;   };
    std::complex <double> i(0,1);
    double a = lambda*theta;
    std::complex<double> d  = sqrt( pow(rho*eta*u*i-bj,2) - 
		eta_sqr*(2*uj*u*i-pow(u,2)) );
    std::complex<double> g  = (bj - rho*eta*u*i+d)/(bj-rho*eta*u*i-d);
    std::complex<double> C  = r*u*i*tau+(a/eta_sqr)*((bj-rho*eta*u*i+d)*tau
		-2.0*log((1.0-g*exp(d*tau))/(1.0-g)));
    std::complex<double> D  = (bj-rho*eta*u*i+d)/eta_sqr 
		* ( (1.0-exp(d*tau))/(1.0-g*exp(d*tau)) );
    std::complex<double> f1 = exp(C+D*v+i*u*x);
    std::complex<double> F  = exp(-u*i*log(K))*f1/(i*u);
	return std::real(F); 
};};


////Gauss-Laguerre
 static double x[150] = {
0.542626,0.370395,0.230993,0.124405,0.050618,0.009607,
0.747706,0.985657,1.256506,1.560281,
1.897016,2.266749,2.669519,3.105372,3.574356,4.076522,4.611926,5.180628,5.782691,6.418183,
7.087175,7.789743,8.525966,9.295927,10.099714,10.937420,11.809140,12.714974,13.655029,14.629413,
15.638241,16.681630,17.759705,18.872593,20.020428,21.203347,22.421494,23.675016,24.964068,26.288808,
27.649400,29.046015,30.478829,31.948023,33.453785,34.996309,36.575795,38.192451,39.846489,41.538130,
43.267601,45.035137,46.840980,48.685380,50.568593,52.490886,54.452532,56.453814,58.495024,60.576461,
62.698436,64.861268,67.065287,69.310833,71.598257,73.927922,76.300199,78.715476,81.174150,83.676631,
86.223344,88.814727,91.451232,94.133325,96.861491,99.636228,102.458051,105.327495,108.245112,
111.211471,114.227165,117.292806,120.409026,123.576484,126.795859,130.067858,133.393212,136.772683,
140.207058,143.697158,147.243835,150.847975,154.510500,158.232369,162.014582,165.858181,169.764254,
173.733934,177.768406,181.868909,186.036738,190.273249,194.579863,198.958070,203.409435,207.935600,
212.538294,217.219337,221.980648,226.824250,231.752283,236.767009,241.870826,247.066279,252.356069,
257.743073,263.230360,268.821204,274.519111,280.327841,286.251434,292.294241,298.460962,304.756684,
311.186934,317.757729,324.475644,331.347894,338.382420,345.588004,352.974402,360.552509,368.334558,
376.334374,384.567687,393.052539,401.809802,410.863871,420.243586,429.983502,440.125674,450.722256,
461.839398,473.563384,486.010824,499.346798,513.820397,529.844189,548.212947,570.989411};

static double w[150] = {
0.188652,0.155814,0.122992,0.090185,0.057392,0.024654,0.221512,0.254395,0.287307,0.320249,0.353228,
0.386245,0.419304,0.452410,0.485566,0.518775,0.552043,0.585372,0.618766,0.652230,0.685767,0.719382,
0.753078,0.786859,0.820731,0.854696,0.888760,0.922927,0.957201,0.991586,1.026088,1.060711,1.095460,
1.130339,1.165353,1.200509,1.235809,1.271261,1.306869,1.342638,1.378575,1.414685,1.450973,1.487446,
1.524110,1.560971,1.598036,1.635311,1.672802,1.710518,1.748464,1.786649,1.825079,1.863763,1.902708,
1.941923,1.981416,2.021196,2.061272,2.101654,2.142349,2.183370,2.224725,2.266426,2.308483,2.350908,
2.393712,2.436908,2.480508,2.524525,2.568974,2.613867,2.659220,2.705048,2.751367,2.798193,2.845543,
2.893437,2.941891,2.990927,3.040563,3.090823,3.141728,3.193301,3.245567,3.298552,3.352283,3.406789,
3.462099,3.518244,3.575258,3.633176,3.692034,3.751871,3.812729,3.874651,3.937683,4.001874,4.067277,
4.133945,4.201939,4.271320,4.342156,4.414519,4.488484,4.564134,4.641557,4.720849,4.802110,4.885451,
4.970991,5.058860,5.149198,5.242157,5.337903,5.436618,5.538501,5.643771,5.752670,5.865463,5.982447,
6.103950,6.230339,6.362026,6.499475,6.643209,6.793823,6.951998,7.118515,7.294276,7.480333,7.677916,
7.888481,8.113766,8.355859,8.617311,8.901266,9.211659,9.553501,9.933297,10.359677,10.844413,11.404090,
12.063035,12.858816,13.853485,15.159514,17.010660,20.015419,26.688256};

inline double heston_Pj(double S, double K, double r, double v, double tau, double eta, 
                 double lambda,  double rho, double theta, int j){

	//FILE   *out_abs= fopen("abs.res","w"), *out_wexp= fopen("wexp.res","w");
    double y=log(S);
    struct heston_parms parms = { K, y, r, v, tau, lambda, theta, rho,eta, j};
    size_t n=10000;

	int i;

	const int       glsteps = 150; /* <= 150 for double precision */
   
	/*integral Gauss-Laguerre*/
	double sum = 0, u;
	for (i=glsteps-1;i>=0;--i)
	  {
		u = x[i];
		sum += w[i] * heston_integrand_j(u, &parms);
	  }

    return 0.5 + sum/PI;
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
				       const double& K,       // strike (exercise) price,
				       const double& r,       // interest rate
				       const double& sigma,   // volatility 
				       const double& time) {  // time to maturity 
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt; 
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

double finum::heston_call_option_price(const double& S, const double& K, const double& r, const double& v, 
				const double& tau, const double& rho,   
				const double& lambda,  const double& theta, const double& eta){
    double P1 = heston_Pj(S,K,r,v,tau,eta, lambda,rho,theta,1);
    double P2 = heston_Pj(S,K,r,v,tau,eta, lambda,rho,theta,2);
    double C=S*P1-K*exp(-r*tau)*P2;
    return C;
};



double finum::option_price_call_merton_jump_diffusion( const double& S, 
						const double& X,
						const double& r,
						const double& sigma, 
						const double& time_to_maturity,
						const double& gamma,
						const double& Jbar,
						const double& beta) {
    const int MAXN=50;
    double tau=time_to_maturity;
    double sigma_sqr = sigma*sigma; 
    double gammaprime = gamma * (1+Jbar);
    double alpha_plus_halfbeta = log(1+Jbar);
    double c = exp(-gammaprime*tau)*option_price_call_black_scholes(S,X,r-gamma*Jbar,sigma,tau);
    double log_n = 0;
    for (int n=1;n<=MAXN; ++n) {
	log_n += log(double(n));
	double sigma_n = sqrt( sigma_sqr+n*beta/tau );
	double r_n = r-gamma*Jbar+n*alpha_plus_halfbeta/tau;
	c += exp(-gammaprime*tau+n*log(gammaprime*tau)-log_n)*
	    option_price_call_black_scholes(S,X,r_n,sigma_n,tau);
    };
    return c;
};


//================== P R I C E === N O N E ==============================================================================
void Progonka_in_place(const  std::vector<double>& low, 
					   const  std::vector<double>& diag, const  std::vector<double>& up,  std::vector<double>& rightside_result)
{
      int N = diag.size();
      if((low.size()!=N-1)||(low.size()!=N-1)||(rightside_result.size()!=N)) {
		  std::cout << "les tailles des vecteurs sont incorrectes!" << std::endl;
		  std::exit(1);
      }

      std::vector<double> A(N-1), B(N-1);

      // descente
      const std::vector<double> & rightside = rightside_result;
      A[0] = -up[0]/diag[0];
      B[0] = rightside[0]/diag[0];
      double denom;
      for(int i=1;i<N-1;i++)
      {
        denom = diag[i] + low[i-1]*A[i-1];
        A[i] = -up[i]/denom;
        B[i] = (rightside[i]-low[i-1]*B[i-1])/denom;
      };

      // remontee
      std::vector<double> & result = rightside_result;
      result[N-1] = (rightside[N-1] - low[N-2]*B[N-2])/(diag[N-1] + low[N-2]*A[N-2]);
      for(int i=N-2;i>=0;i--)
      {
        result[i] = A[i]*result[i+1] + B[i];
      };
}



void Progonka( std::vector<double>::const_iterator low_begin,  
			  std::vector<double>::const_iterator diag_begin,  std::vector<double>::const_iterator up_begin, 
			std::vector<double>::const_iterator rs_begin,  std::vector<double>::iterator result_begin,  std::vector<double>::iterator result_end)
{
      int N = result_end - result_begin;

      std::vector<double> A(N-1), B(N-1);

      // descente
       std::vector<double>::const_iterator low_iter = low_begin;
       std::vector<double>::const_iterator diag_iter = diag_begin;
       std::vector<double>::const_iterator up_iter = up_begin;
       std::vector<double>::const_iterator rs_iter = rs_begin;


      A[0] = -(*up_iter)/(*diag_iter);
      B[0] = (*rs_iter)/(*diag_iter);
      double denom;
      for(int i=1;i<N-1;i++)
      {
        diag_iter++;
        up_iter++;
        rs_iter++;

        denom = (*diag_iter) + (*low_iter)*A[i-1];
        A[i] = -(*up_iter)/denom;
        B[i] = ((*rs_iter) - (*low_iter)*B[i-1])/denom;

        low_iter++;
      };


      // remontee
       std::vector<double>::iterator result_iter = --result_end;
      diag_iter++;
      rs_iter++;
      *result_iter = ((*rs_iter) - (*low_iter)*B[N-2])/((*diag_iter) + (*low_iter)*A[N-2]);
      result_iter--;

      for(int i = N - 2; i >= 0; i--)
      {
         *(result_iter) = A[i]*(*(result_iter+1)) + B[i];
         result_iter--;
      };
  }
 

int finum::ImplPriceNone(double s0, Vanilla Type, double Maturity,  double Strike, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace)
{ 
  return 0;
}


int finum::TreePriceNone(double s0, Vanilla Type, double Maturity,  double Strike, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime, double lambda)
{ 
  return 0;
}


int finum::CNPriceNone(double s0, Vanilla Type, double Maturity,  double Strike, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace)
{ 
  return 0;
}

//===================== U P === A N D === O U T ============================================================
int finum::TreePriceUpOut(double s0,  Vanilla Type, double Maturity,  double Strike, double barrier,  double rebate, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime, double lambda)
{ 
  return 0;
}
int finum::ImplPriceUpOut(double s0,  Vanilla Type, double Maturity,  double Strike, double barrier,  double rebate, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace)
{ 
	return 0;
}



/*-------------------------------------------------------------------*/
#define TailleTable 32 
 
/* pour le generateur direct */  
#define m1 2147483563 
#define m2 2147483399 
#define a1 40014 
#define a2 40692 
#define qa1 53668 
#define qa2 52774 
#define ra1 12211 
#define ra2 3791 
#define hom (1+(m1-1)/TailleTable) 
 


#define MAXBIT 30
#define MAXDIM 160

/**
    @brief numerical inverse of the Gaussian c.d.f.
*/
double finum::Inverse_erf_Moro(double x)    
{
  static double a[4] = {
    2.50662823884, 
    -18.61500062529, 
    41.39119773534, 
    -25.44106049636};
  static double b[4] = {
    -8.47351093090,
    23.08336743743 , 
    -21.06224101826, 
    3.13082909833};
  static double c[9] ={
    0.337475482276147,
    0.9761690190917186,
    0.1607979714918209,
    0.0276438810333863,
    0.0038405729373609,
    0.0003951896511919,
    0.0000321767881768,
    0.000002888167364,
    0.0000003960315187};
  
  double u,r;
  
  u=x-0.5;
  if (fabs(u)<0.42) 
    {
      r=u*u;
      r=u*(((a[3]*r+a[2])*r+a[1])*r+a[0])/((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.0);
      return(r);
    }

  r=x;
  if(u>0.0) 
    r=1.0-x;
  r=log(-log(r));
  r=c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*c[6]+r*(c[7]+r*c[8]))))));
  if (u<0.0) 
    r=-r;
  return(r);
}

	

void finum::tsb2_seq(int *n, float x[])
{
    int                     j, k, l;
    unsigned long           i, im, ipp;

    static unsigned long    in;
    static float            fac;
    static unsigned long    ix[MAXDIM + 1];
    static unsigned long    *iu[MAXBIT + 1];
    /*
        static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
        static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4};
        static unsigned long iv[MAXDIM*MAXBIT+1]={
                0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
	*/

        static unsigned long  mdeg[MAXDIM+1] = {0 ,1 ,2 ,3 ,3 ,4 ,4 ,5 ,5 ,5 ,5 ,5 ,5 ,6 ,6 ,6 ,6 ,6 ,6 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10};
 
        static unsigned long  ip[MAXDIM+1] = {0 ,0 ,1 ,1 ,2 ,1 ,4 , 2 ,4 ,7 ,11 ,13 ,14 ,1 ,13 ,16 ,19 ,22 ,25 ,1 ,4 ,7 ,8 ,14 ,19 ,21 ,28 ,31 ,32 ,37 ,41 ,42 ,50 ,55 ,56 ,59 ,62 ,14 ,21 ,22 ,38 ,47 ,49 ,50 ,52 ,56 ,67 ,70 ,84 ,97 ,103 ,115 ,122 ,8 ,13 ,16 ,22 ,25 ,44 ,47 ,52 ,55 ,59 ,62 ,67 ,74 ,81 ,82 ,87 ,91 ,94 ,103 ,104 ,109 ,122 ,124 ,137 ,138 ,143 ,145 ,152 ,157 ,167 ,173 ,176 ,181 ,182 ,185 ,191 ,194 ,199 ,218 ,220 ,227 ,229 ,230 ,234 ,236 ,241 ,244 ,253 ,4 ,13 ,19 ,22 ,50 ,55 ,64 ,69 ,98 ,107 ,115 ,121 ,127 ,134 ,140 ,145 ,152 ,158 ,161 ,171 ,181 ,194 ,199 ,203 ,208 ,227 ,242 ,251 ,253 ,265 ,266 ,274 ,283 ,289 ,295 ,301 ,316 ,319 ,324 ,346 ,352 ,361 ,367 ,382 ,395 ,398 ,400 ,412 ,419 ,422 ,426 ,428 ,433 ,446 ,454 ,457 ,472 ,493 ,505 ,508};
 
 
        static unsigned long  iv[MAXDIM*MAXBIT+1]={0 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,3 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,5 ,7 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,15 ,11 ,13 ,15 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,17 ,13 ,5 ,23 ,25 ,31 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,51 ,61 ,47 ,53 ,9 ,47 ,57 ,61 ,53 ,55 ,51 ,59 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,85 ,67 ,113 ,29 ,11 ,109 ,25 ,29 ,21 ,71 ,75 ,83 ,97, 101 ,127 ,111 ,119 ,99 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,63 ,255 ,79 ,145 ,149 ,191 ,173 ,27 ,151 ,191 ,229 ,73 ,81 ,33 ,37 ,191 ,143 ,167 ,155 ,193 ,241 ,209 ,249 ,233 ,221 ,205 ,245 ,213 ,255 ,207 ,199 ,231 ,227 ,211 ,251 ,219 ,235 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,127 ,257 ,465 ,147 ,183 ,449 ,181 ,303 ,147 ,491 ,101 ,425 ,209 ,35 ,335 ,445 ,461 ,389 ,153 ,65 ,113 ,81 ,121 ,105 ,93 ,77 ,117 ,85 ,383 ,303 ,311 ,343 ,347 ,267 ,371 ,275 ,323 ,465 ,409 ,473 ,477 ,429 ,389 ,453 ,485 ,501 ,447 ,479 ,487 ,387 ,419 ,443 ,459 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,255 ,771 ,721 ,445 ,411 ,193 ,949 ,295 ,947 ,235 ,109 ,169 ,913 ,39 ,923 ,701 ,589 ,133 ,889 ,67 ,115 ,83 ,123 ,107 ,727 ,647 ,607 ,767 ,893 ,1005 ,981 ,821 ,345 ,265 ,369 ,273 ,321 ,209 ,153 ,217 ,221 ,173 ,133 ,197 ,229 ,245 ,575 ,671 ,727 ,635 ,539 ,563 ,675 ,993 ,801 ,1009 ,945 ,785 ,985 ,857 ,969 ,841 ,873 ,937 ,893 ,925 ,781 ,909 ,845 ,877 ,941 ,837 ,997 ,805 ,917 ,981 ,799 ,927 ,863 ,783 ,1007 ,815 ,839 ,807 ,1015 ,823 ,951 ,791 ,855 ,899 ,835 ,915 ,979 ,891 ,827 ,955 ,923 ,987 ,779 ,971 ,811 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511 ,511};
 

    if (*n < 0) 
    {
       for (k = 1; k <= MAXDIM; k++)  ix[k] = 0;

       in=0;

       if (iv[1] != 1) return;

       fac=1.0/(1L << MAXBIT);

       for (j = 1, k = 0; j <= MAXBIT; j++, k += MAXDIM) 
          iu[j] = &iv[k];

       for (k = 1; k <= MAXDIM; k++) 
       {
          for (j = 1; j <= mdeg[k]; j++) 
             iu[j][k] <<= (MAXBIT - j);

          for (j = mdeg[k] + 1; j <= MAXBIT; j++) 
          {
             ipp = ip[k];

             i = iu[j-mdeg[k]][k];

             i ^= (i >> mdeg[k]);

             for (l = mdeg[k] - 1; l >= 1; l--) 
             {
                if (ipp & 1)  i ^= iu[j-l][k];

                ipp >>= 1;
             }

             iu[j][k]=i;
          }

       }

    } 
    else 
    {
       im=in;

       for (j = 1; j <= MAXBIT; j++) 
       {
          if (!(im & 1)) break;

          im >>= 1;
       }

       if (j > MAXBIT) return;

       im=(j-1)*MAXDIM;

       for (k = 1; k <= MIN(*n, MAXDIM); k++) 
       {
          ix[k] ^= iv[im+k];

          x[k]=ix[k]*fac;
       }

       in++;
    }


}

//===================== TESTS============================================================
void test_heston(){
    double S=100;
    double K=100;
    double r=0.01;
    double v=0.01;
    double tau=0.5;
    double rho=0;
    double lambda=2;
    double theta=0.01;
    double eta=0.01;
    cout << " Heston call price " << 
		heston_call_option_price( S, K, r, v, tau, rho, lambda, theta, eta) 
		<< endl;
};


void test_merton_jump_diff_call(){
    double S=100; double K=100; double r=0.05;
    double sigma=0.3; 
    double time_to_maturity=1;
    double gamma=0.5;
    double Jbar=0.5;
    double beta=0.5*0.5;
    cout << " Merton Jump diffusion call = " 
	<< option_price_call_merton_jump_diffusion(S,K,r,sigma,time_to_maturity,
			gamma,Jbar,beta)
	 << endl;
};

//}///int main(int argc, char *argv[])extern "C" {
int main() 
 {

  double price,s0, Maturity,  Strike, eqr,  eqdiv, sigma_ref, volimpl,impliedDelta,  impliedGamma, barrierup, barrierdown, rebate;
  price=.2;
  s0=1.;
  Maturity=1.;
  Strike=1.;
  rebate=0.;
  eqr=0.;
  eqdiv=0.;
  sigma_ref=.2;
  volimpl=0.3;
  barrierup=1.1;
  barrierdown=0.9;

  ImpliedVol(price, s0, Maturity,  Strike, eqr,  eqdiv, sigma_ref, &volimpl,  &impliedDelta,  &impliedGamma);

  Put_BlackScholes_73(s0,  Maturity,  Strike,  eqr,  eqdiv,  volimpl, &price, &impliedDelta,  &impliedGamma) ;

  cout << volimpl << " " << price << " " <<  impliedDelta << " " << impliedGamma << " " << endl ;
	
    Call_BlackScholes_73(s0,  Maturity,  Strike,  eqr,  eqdiv,  sigma_ref, &price, &impliedDelta,  &impliedGamma) ;
 cout << price << " " <<  impliedDelta <<endl ;

 
 CallUpOut_ReinerRubinstein(s0,Strike,barrierup,rebate,Maturity,eqr,eqdiv, sigma_ref, &price, &impliedDelta);
	  cout << "CallUpOut_ReinerRubinstein " <<  price << " " <<  impliedDelta <<endl ;

CallDownOut_ReinerRubinstein(s0,Strike,barrierdown,rebate,Maturity,eqr,eqdiv, sigma_ref, &price, &impliedDelta);
	  cout << "CallDownOut_ReinerRubinstein " <<  price << " " <<  impliedDelta <<endl ;


	 Call_KunitomoIkeda_91(s0,barrierdown,barrierup,rebate,Strike,Maturity,eqr,eqdiv, sigma_ref, &price, &impliedDelta);
	  cout << "Call_KunitomoIkeda_91 " <<  price << " " <<  impliedDelta <<endl ;

	   price=.2;
  s0=100;
  Maturity=1.;
  Strike=100.;
  eqr=0.1;
  eqdiv=0.;
  sigma_ref=.2;
  volimpl=0.3;
  
Fixed_CallLookback_ConzeWiswanathan(s0,s0, Strike, Maturity,eqr,eqdiv, sigma_ref, &price, &impliedDelta);
	  cout << "Fixed_CallLookback_ConzeWiswanathan " <<  price << " " <<  impliedDelta <<endl ;


	test_heston();
	test_merton_jump_diff_call();    

 int nSpace = 3;

  std::vector<double> Price(nSpace + 1);

  // ces vecteurs representent la matrice tridiagonale du systeme:
  std::vector<double> low(nSpace);      // diagonale inferieure
  std::vector<double> diag(nSpace + 1); // diagonale principale
  std::vector<double> up(nSpace);       // diagonale superieure


  // conditions au bord: propagation des valeurs initiales:
  up[0] = 0;
  diag[0] = 1;

  low[nSpace - 1] = 0;
  diag[nSpace] = 1;

  for(int i = 1; i < nSpace; i++)
    {
      Price[i] = 0;
    }

  for(int i = 1; i < nSpace; i++)
    {
      up[i] = 0;
      diag[i] = 1;
      low[i-1] = 0;
    }

  Progonka_in_place(low, diag, up, Price);

  for(int i = 1; i < nSpace; i++)
    {
		cout << i << " "<< up[i]<< " " << diag[i]<< " " << low[i-1]  << " " <<  
		Price[i]  << endl;
    }

  cout << Inverse_erf_Moro(0.5)<< endl;

  std::ofstream outfile("sortie.dat");
  float tsb2[MAXDIM+1];
  int N=5;
  //initialization of tsb2
  int m=-1;
  tsb2_seq(&m, tsb2);
  m=MAXDIM;

  for(int i=0;i<N;i++)//draw N quasi-uniform point in [0;1]^m
    {
      tsb2_seq(&m,tsb2);//tsb2[1..m]=quasi-uniform point in [0;1]^m
	  outfile <<  tsb2[1] << " " << tsb2[2] << std::endl;
    }
  
  return 0;
}

