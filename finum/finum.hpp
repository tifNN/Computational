#include <cmath>
#include <ctime>
#include <vector>
#include <iterator>

/*! \mainpage FINUM (S Crépey, Evry University)
 *
 A dll of C/C++ basic financial functions for pedagogical purposes
 */

#ifndef _MAFIN_H
#define _MAFIN_H

typedef unsigned char boolean;

#define false 0
#define true 1

#define PI 3.14159265358979
#define DEUXPI 6.283185306
#define PI2 1.570796326794
#define INC 0.01

#define PRECISION 1.0e-7 /*Precision for the bivariate Gaussian cumulative distribution function*/
#define PRICE_ACCURACY 1e-5 /*Precision for the implied volatility search*/
#define FAC_SLOPE 0.01


#define MAXLOOPS 5000

#define MAX(A,B) ( (A) > (B) ? (A):(B) )
#define MIN(A,B) ( (A) < (B) ? (A):(B) )
#define SQR(X) ((X)*(X))

/**
    @brief Enumerated type for Vanilla options
*/
enum Vanilla {CallEuro=1,PutEuro=2,CallAmer=3,PutAmer=4};


/**
	@brief  Gaussian density
*/
double nd(double x);


/**
    @brief Gaussian c.d.f.
 */
double N(double x);
/**
    @brief Library Functions
*/

/**
    @brief Thomas algorithm for tridiagonal systems Ax=b
	
	@param low	lower diagonal of A
    @param diag	 diagonal of A
    @param up	upper diagonal of A
 */
void Progonka_in_place(const std::vector<double>& low, const std::vector<double>& diag, const
 std::vector<double>& up, std::vector<double>& rightside_result) ;

/**
    @brief Sub-routine called by Progonka_in_place
*/
void Progonka( std::vector<double>::const_iterator low_begin,  
			  std::vector<double>::const_iterator diag_begin,  std::vector<double>::const_iterator up_begin,
               std::vector<double>::const_iterator rs_begin,  
			   std::vector<double>::iterator result_begin,  std::vector<double>::iterator result_end);

namespace finum 
{
/**
    @brief BS call formula
*/
	int Call_BlackScholes_73(double s, double Maturity, double Strike, double r,double divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma);
/**
    @brief BS put formula
*/
	int Put_BlackScholes_73(double s, double Maturity, double Strike, double r,double divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma);
/**
    @brief BS implied volatility solver (dichotomy)
*/
	int ImpliedVol(double price,double s0,double Maturity, double Strike,double eqr, double eqdiv,double sigma_ref,double *volimpl, double *impliedDelta, double *impliedGamma);

/**
    @brief BS Fixed_CallLookback_ConzeWiswanathan. Source: Premia5
*/
	int Fixed_CallLookback_ConzeWiswanathan(double s, double s_max, double k, double t, double r,
              double divid, double sigma, double *ptprice, double *ptdelta);
/**
    @brief BS Call with Barrier Up Out. Source: Premia5
*/
	int CallUpOut_ReinerRubinstein(double s,double k,double l,double rebate,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);

/**
    @brief BS Call with Barrier Down Out. Source: Premia5
*/
	int CallDownOut_ReinerRubinstein(double s,double k,double l,double rebate,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);

/**
    @brief BS Call with Double Out Barrier. Source: Premia5
*/
	int Call_KunitomoIkeda_91(double s,double L,double U,double Rebate,double Strike,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);

/**
    @brief Heston call formula (Source: Financial Numerical Recipes, 
	Bernt Odegaard)
*/
	double heston_call_option_price(const double& S, const double& K, const double& r, const double& v, 
				const double& tau, const double& rho,   
				const double& lambda,  const double& theta, const double& eta);
	
	/**
    @brief Merton call formula (Source: Financial Numerical Recipes, 
	Bernt Odegaard)
*/
	double option_price_call_merton_jump_diffusion( const double& S, 
						const double& X,
						const double& r,
						const double& sigma, 
						const double& time_to_maturity,
						const double& gamma,
						const double& Jbar,
						const double& beta);
	/**
    @brief BS vanilla implicit pricer
	*/
	int  ImplPriceNone(double s0, Vanilla Type, double Maturity,  double Strike, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace);

	/**
    @brief BS vanilla tree pricer
*/
	int TreePriceNone(double s0, Vanilla Type, double Maturity,  double Strike, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime, double lambda);
/**
    @brief BS vanilla Crank-Nicholson pricer
*/
	int  CNPriceNone(double s0, Vanilla Type, double Maturity,  double Strike, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace);

/**
    @brief BS up-and-out barrier tree pricer
*/
	int TreePriceUpOut(double s0,  Vanilla Type, double Maturity,  double Strike, double barrier,  double rebate, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime, double lambda);

/**
    @brief BS up-and-out barrier implicit pricer
*/
	int  ImplPriceUpOut(double s0,  Vanilla Type, double Maturity,  double Strike, double barrier,  double rebate, double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace);
/**
    @brief numerical inverse of the Gaussian c.d.f.
*/
	double Inverse_erf_Moro(double x);
/**
    @brief Sobol sequence in dimension 1 to 160
*/
	void tsb2_seq(int *n, float x[]);
}

/**
    @brief Tests (for console compilation)
*/
//int main();

#define int_dll  __declspec(dllexport) int __stdcall
#define void_dll  __declspec(dllexport) void  __stdcall
#define double_dll __declspec(dllexport) double __stdcall
extern "C" {


int_dll Call_BlackScholes_73(double s, double Maturity, double Strike, double r,double divid, 
							 double sigma, double *ptprice, double *ptdelta, 
							 double *ptgamma){return finum::Call_BlackScholes_73( s,  
							 Maturity,  Strike,  r,divid, sigma, ptprice, ptdelta, ptgamma) ;} 

int_dll Put_BlackScholes_73(double s, double Maturity, double Strike, 
							double r,double divid, double sigma, double *ptprice,
							double *ptdelta, double *ptgamma)
							{return finum::Put_BlackScholes_73( s,  Maturity,  Strike,  r,divid, sigma, ptprice, ptdelta, ptgamma) ;} 

int_dll ImpliedVol(double price,double s0,double Maturity, double Strike,double eqr, 
				   double eqdiv,double sigma_ref,double *volimpl, 
				   double *impliedDelta, double *impliedGamma)
	{return finum::ImpliedVol(  price, s0, Maturity,   
	Strike, eqr,   eqdiv,  sigma_ref, volimpl,  impliedDelta,  impliedGamma) ;} 

int_dll Fixed_CallLookback_ConzeWiswanathan(double s, double s_max, double k, double t, double r,
              double divid, double sigma, double *ptprice, double *ptdelta)
     {return finum::Fixed_CallLookback_ConzeWiswanathan(s, s_max, k, t, r, divid, sigma, ptprice, ptdelta) ;} 

	int_dll CallUpOut_ReinerRubinstein(double s,double k,double l,double rebate,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
	{return finum::CallUpOut_ReinerRubinstein( s, k, l, rebate, t, r, divid, sigma, ptprice, ptdelta) ;} 

int_dll CallDownOut_ReinerRubinstein(double s,double k,double l,double rebate,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
	{return finum::CallDownOut_ReinerRubinstein( s, k, l, rebate, t, r, divid, sigma, ptprice, ptdelta) ;} 

int_dll Call_KunitomoIkeda_91(double s,double L,double U,double Rebate,double Strike,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta){return finum::Call_KunitomoIkeda_91( s, L, U, Rebate, Strike, t, r, divid, sigma,ptprice,ptdelta) ;} 

double_dll heston_call_option_price(const double& S, const double& K, const double& r, const double& v, 
				const double& tau, const double& rho,  
				const double& lambda,  const double& theta, const double& eta)
					{return finum::heston_call_option_price(S, K, r, v, 
				tau, rho,  
				lambda,  theta, 
				eta) ;} 

	
double_dll option_price_call_merton_jump_diffusion( const double& S, 
						const double& X,
						const double& r,
						const double& sigma, 
						const double& time_to_maturity,
						const double& gamma,
						const double& Jbar,
						const double& beta)
	{return finum::option_price_call_merton_jump_diffusion(  S, 
						 X,
						 r,
						 sigma, 
						 time_to_maturity,
						 gamma,
						 Jbar,
						 beta) ;} 

int_dll ImplPriceNone(double s0, Vanilla Type, double Maturity,  
					  double Strike, double r,double  divid, double sigma,
					  double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace){return finum::ImplPriceNone(s0, Type, Maturity, Strike,  r,
					  divid, sigma, ptprice, ptdelta, 
					 ptgamma, nTime, nSpace) ;} 

int_dll TreePriceNone(double s0, Vanilla Type, double Maturity,  double Strike, double r,
					  double  divid, double sigma, double *ptprice, double *ptdelta, 
					  double *ptgamma, int nTime, double lambda){return finum::TreePriceNone(s0, Type, Maturity, Strike,  r,
					  divid, sigma, ptprice, ptdelta, 
					 ptgamma, nTime,  lambda) ;} 

int_dll  CNPriceNone(double s0, Vanilla Type, double Maturity,  double Strike, double r,double  divid, double sigma, 
					 double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace)
					{return finum::CNPriceNone(s0, Type, Maturity, Strike,  r,
					  divid, sigma, ptprice, ptdelta, 
					 ptgamma, nTime, nSpace) ;} 

int_dll TreePriceUpOut(double s0,  Vanilla Type, double Maturity,  double Strike, double barrier,  double rebate, 
					   double r,double  divid, double sigma, double *ptprice, double *ptdelta, double *ptgamma, 
					   int nTime, double lambda){return finum::TreePriceUpOut
					   ( s0,   Type, Maturity,   Strike,  barrier,   rebate,  r, divid,  
					   sigma, ptprice, ptdelta, ptgamma,  nTime,  lambda) ;} 

int_dll  ImplPriceUpOut(double s0,  Vanilla Type, double Maturity,  double Strike, double barrier,  
						double rebate, double r,double  divid, double sigma, 
						double *ptprice, double *ptdelta, double *ptgamma, int nTime,int nSpace)
						{return finum::ImplPriceUpOut( s0,    Type, Maturity,   Strike,  barrier,   rebate,  r, divid,  
					   sigma, ptprice, ptdelta, ptgamma,  nTime, nSpace) ;} 

double_dll Inverse_erf_Moro(double x){return finum::Inverse_erf_Moro(x) ;} 

void_dll tsb2_seq(int *n, float x[]){finum::tsb2_seq(n, x) ;} 

//int_dll main(){return finum::main() ;} 
}

#endif