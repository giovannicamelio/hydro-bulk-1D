/// @file
/// @brief hydro-bulk-1D: 1-dimensional general relativistic hydrodynamics
///        with bulk viscosity
/// @details
///    Unless otherwise specified, the code units are c = G = M☉ = kB = 1.
///
///    The output profile has the following format:
///    r|x,   ρ,   Wv,   u,   [Ye|Π],   [Yμ],   [mg],   cs
///
///    The logfile has the following format:
///    t,   ρ,   Wv,   u,   [Ye|Π],   [Yμ],   cs,   [Mg],   [Mb],   ρend
/// @author Giovanni Camelio
/// @date 2022
/// @copyright 2022 Giovanni Camelio.
///            This project is released under the MIT License.

// INCLUDES ********************************************************************

// for fopen, fclose, printf, fprintf
#include <stdio.h>

// for sqrt, pow, cbrt, fmax, fmin, fabs, exp, log, copysign, HUGE_VAL
#include <math.h>

// for bool, true, false
#include <stdbool.h>

// contains the code settings
#include "parameters.h"

// MATHEMATICAL CONSTANT *******************************************************

#ifndef M_PI
/// @brief M_PI is not defined in standard C
#define M_PI (3.14159265358979323846)
#endif

// PHYSICAL CONSTANTS IN CGS UNITS *********************************************

/// @brief nuclear saturation density in cgs units [g/cm³]
#define RHON_CGS 2.68e14

/// @brief solar mass in cgs units [g]
#define MSOL_CGS 1.98847e33

/// @brief neutron mass in cgs units [g]
#define MN_CGS 1.67492749804e-24

/// @brief Boltzmann constant in cgs units [erg/K]
#define KB_CGS 1.380649e-16

/// @brief 1 MeV in gcs units [erg]
#define MEV_CGS 1.602176634e-6

/// @brief speed of light in cgs units [cm/s]
#define C_CGS 2.99792458e10

/// @brief gravitational constant in cgs units [cm³·g⁻¹·s⁻²]
#define G_CGS 6.67430e-8

// PHYSICAL CONSTANTS IN CODE UNITS [c=G=M☉=kB=1] ******************************

/// @brief neutron mass in code units [M☉]
#define MN (MN_CGS/MSOL_CGS)

/// @brief 1 cm in code units [c²·G⁻¹·M☉⁻¹]
#define CM (C_CGS*C_CGS/G_CGS/MSOL_CGS)

/// @brief 1 sec in code units [c³·G⁻¹·M☉⁻¹]
#define SEC (C_CGS*CM)

/// @brief nuclear saturation density in code units [M☉⁴·G³/c⁶]
#define RHON (RHON_CGS/(CM*CM*CM*MSOL_CGS))

/// @brief 1 erg in code units [M☉·c²]
#define ERG (MN/(MN_CGS*C_CGS*C_CGS))

/// @brief 1 K in code units [M☉·c²/kB]
#define KELVIN (KB_CGS*ERG)

// DERIVED PARAMETERS **********************************************************
// BEWARE: the constants in this section are set automatically
// don't touch them unless you know what you are doing

/// @brief control-flow flag, true if the equations have gravity
/// @details true (gravitating) / false (non gravitating)
#define GRAVITATING (PROBLEM >= CLOUD)

/// @brief number of ghost points on each side
/// @details integer, NG > 0
#define NG (ORDER_SPACE)

/// @brief number of independent equations/variables
/// @details integer, N0 > 0
#define N0 (FLUID == BULK ? 4 : 2 + NSPECIES)

/// @brief total number of variables in each grid point
/// @details integer, D0 > 0
#define D0 (N0 + (GRAVITATING ? 1 : 0))

/// @brief total dimension of the grid vector (including ghosts)
/// @details integer, D1 = @ref N1 + 2 * @ref NG
#define D1 (N1 + 2*NG)

/// @brief minimal value for the internal energy [M☉⁻²]
#define MIN_U (FLUID == DUST ? 0. : K0*pow(MIN_RHO, GAM - 1.))

/// @brief minimal value for the energy density [M☉⁻²]
#define MIN_EPS (MIN_RHO*(1. + MIN_U))

// @brief minimal value for the pressure
#define MIN_P ((GAM - 1.)*K0*pow(MIN_RHO, GAM))

// FUNCTIONS *******************************************************************

/// @brief error handling: log and return the error
/// @param[in] str error string to log
/// @return 1 (which means failure since it is ≠ 0)
int error (const char* str) {
   printf("ERROR: %s\n",str);
   return 1;
}

/// @brief root finder with Brent's method
/// @param[in] minx minimal value for the search
/// @param[in] maxx maximal value for the search
/// @param[in] f function to find the root of
/// @param[in] par conservative and auxiliary variables
/// @param[out] out primitive variables
/// @param[in] x first guess
/// @return error code (0 on success)
int find_root (
   const double minx,
   const double maxx,
   double (*f)(const double, const double*, double*),
   const double *par,
   double *out,
   const double x)
{
   // initial guess
   double a= fmax(minx, fmin(maxx, x) - TOL);
   double fa= f(a, par, out);
   if (fa == 0.) return 0;
   double b= fmin(maxx, a + TOL);
   double fb= f(b, par, out);
   if (fb == 0.) return 0;

   // bracket the root
   const int max_iter= 100;
   for (int i= 0; i < max_iter; i++) {
      if (fa == 0. || fb == 0.) return 0;
      if (fa*fb <= 0.) break;
      else if (fabs(fa) < fabs(fb)) {
         if (a <= minx) return error("find_root: min value reached");
         a= fmax(minx, 2.*a - b);
         fa= f(a, par, out);
      } else {
         if (b >= maxx) return error("find_root: max value reached");
         b= fmin(maxx, 2.*b - a);
         fb= f(b, par, out);
      }
   } // end for (bracket)
   if (fa*fb > 0.) return error("find_root: bracketing failed");

   // Brent's method
   const double abs_tol= 0.5*TOL*(fabs(a) + fabs(b));
   double tmp, s, fs, d= 0.; // I initialize d only to avoid compilation warnings
   // b is the closest to the root
   if (fabs(fa) < fabs(fb)) {
      tmp= b;
      b= a;
      a= tmp;
      tmp= fb;
      fb= fa;
      fa= tmp;
   }
   double c= a, fc= fa;
   enum FLAG {SET, CLEAR} flag= SET;
   for (int i= 0; i < max_iter; i++) {
      if (fabs(b - a) < abs_tol || fb == 0.)
         return 0;
      if (fa != fc && fb != fc) {
         // inverse quadratic interpolation
         s= a*fb*fc/(fa - fb)/(fa - fc)
          + b*fa*fc/(fb - fa)/(fb - fc)
          + c*fa*fb/(fc - fa)/(fc - fb);
      } else {
         // secant
         s= b - fb*(b - a)/(fb - fa);
      }
      // note that in the first iteration flag == SET, hence d is not used
      if (   (4.*s - 3.*a - b)*(s - b) > 0.
          || (flag == SET && fabs(s - b) >= 0.5*fabs(b - c))
          || (flag == CLEAR && fabs(s - b) >= 0.5*fabs(c - d))
          || (flag == SET && fabs(b - c) < abs_tol)
          || (flag == CLEAR && fabs(c - d) < abs_tol)  ) {
         // bisection
         s= 0.5*(a + b);
         flag= SET;
      } else {
         flag= CLEAR;
      }
      fs= f(s, par, out);
      // old old value of the root
      d= c;
      // old value of the root
      c= b;
      fc= fb;
      // new bracket
      if (fa*fs < 0.) {
         b= s;
         fb= fs;
      } else {
         a= s;
         fa= fs;
      }
      // b is the closest to the root
      if (fabs(fa) < fabs(fb)) {
         tmp= b;
         b= a;
         a= tmp;
         tmp= fb;
         fb= fa;
         fa= tmp;
      }
   } // end for
   return error("find_root: brent didn't converge");
}

/// @brief rest mass density from the equilibrium equation of state
/// @param[in] p pressure
/// @param[in] s entropy per baryon
/// @return equilibrium rest mass density
double p2rho(const double p, const double s) {
   // if there is no thermal contribution, directly invert the EOS
   double rho= pow(p/K0/(GAM - 1.), 1./GAM);
   if (s <= 0.) return rho;
   // otherwise, solve the system with Newton-Raphson
   else {
      const int max_iter= 1000;
      double tmp1, tmp2, drho;
      for (int i= 0; i < max_iter; i++) {
         tmp1= K0*(GAM - 1.)*pow(rho, GAM - 1.);
         tmp2= KTH*(GAMTH - 1.)*s*s*pow(rho, GAMTH - 1.);
         drho= (p - (tmp1 + tmp2)*rho)/(GAM*tmp1 + GAMTH*tmp2);
         if (fabs(drho) < TOL*rho) return rho;
         else if (rho == MIN_RHO && drho <= 0.) return rho;
         else rho= fmax(MIN_RHO, rho + drho);
      }
      // being the function to invert well behaved (positive with
      // positive derivative), I do not expect failure.
      // Therefore, and also because I use this function only at
      // initialization where is easier to spot errors, I do not
      // implement error handling, but just resort to the atmosphere
      // in the (theoretically possible) case of failure.
      return MIN_RHO;
   }
}

/// @brief total energy density from the equilibrium equation of state
/// @param[in] rho rest mass density
/// @param[in] s entropy per baryon
/// @return equilibrium total energy density
double rho2eps(const double rho, const double s) {
   if (FLUID == DUST) return rho;
   else return rho + K0*pow(rho, GAM) + KTH*s*s*pow(rho, GAMTH);
}

/// @brief pressure from the equilibrium equation of state
/// @param[in] rho rest mass density
/// @param[in] s entropy per baryon
/// @return equilibrium pressure
double rho2p(const double rho, const double s) {
   if (FLUID == DUST) return 0.;
   else return K0*(GAM - 1.)*pow(rho, GAM) + KTH*(GAMTH - 1.)*s*s*pow(rho, GAMTH);
}

/// @brief compute the inverse bulk viscous coefficient
/// @param[in] rho rest mass density
/// @return 1/χ inverse bulk viscous coefficient
double get_1_over_chi (const double rho) {
   if (NSPECIES == 1) return 1./CHI_BULK;
   const double yeqe= rho/RHON*YE0;
   const double yeqmu= (NSPECIES == 3 ? rho/RHON*YMU0 : 0.);
   return 2.*rho*(KE*yeqe*yeqe + KMU*yeqmu*yeqmu);
}

/// @brief `ultraviolet' sound speed from equation of state
/// @param[in] rho rest mass density
/// @param[in] u internal specific (per unit mass) energy
/// @param[in] ye electron fraction
/// @param[in] ymu muon fraction
/// @return `ultraviolet' speed of sound
double eos2c(const double rho, const double u, const double ye, const double ymu) {
   if (FLUID == DUST) return 0.;
   // partial specific energies
   const double u0= K0*pow(rho, GAM - 1.);
   double ue= 0., umu= 0., tmp1= 0., tmp2= 0.;
   if (FLUID == PERFECT && NSPECIES >= 2) {
      const double yeqe= rho*YE0/RHON;
      const double dye= ye - yeqe;
      ue= KE*dye*dye;
      tmp1= -2.*KE*yeqe*(2.*dye - yeqe);
      tmp2= KE*dye*(dye - 2.*yeqe);
      if (NSPECIES == 3) {
         const double yeqmu= rho*YMU0/RHON;
         const double dymu= ymu - yeqmu;
         umu= KMU*dymu*dymu;
         tmp1= tmp1 - 2.*KMU*yeqmu*(2.*dymu - yeqmu);
         tmp2= tmp2 + KMU*dymu*(dymu - 2.*yeqmu);
      }
   }
   const double uth= fmax(0., u - u0 - ue - umu);
   const double heq= 1. + GAM*u0 + GAMTH*uth;
   const double dpdrho= GAM*(GAM - 1.)*u0 + GAMTH*(GAMTH - 1.)*uth + tmp1;
   const double depsdrho= heq + tmp2;
   double c2= dpdrho / depsdrho;
   if (FLUID == BULK) c2= c2 + get_1_over_chi(rho)/(rho*heq);
   return sqrt(fmax(0., fmin(1., c2)));
}

/// @brief pressure from equation of state
/// @param[in] rho rest mass density
/// @param[in] u internal specific (per unit mass) energy
/// @param[in] ye electron fraction
/// @param[in] ymu muon fraction
/// @return pressure
double eos2p (const double rho, const double u, const double ye, const double ymu) {
   if (FLUID == DUST) return 0.;
   // particle fraction contribution
   double epse= 0., epsmu= 0., pe= 0., pmu= 0.;
   if (FLUID == PERFECT && NSPECIES >= 2) {
      const double yeqe= YE0/RHON*rho;
      const double dye= ye - yeqe;
      epse= KE*dye*dye*rho;
      pe= -2.*rho*yeqe*KE*dye;
      if (NSPECIES == 3) {
         const double yeqmu= YMU0/RHON*rho;
         const double dymu= ymu - yeqmu;
         epsmu= KMU*dymu*dymu*rho;
         pmu= -2.*rho*yeqmu*KMU*dymu;
      }
   }
   // energy density
   const double eps0= K0*pow(rho, GAM);
   const double epsth= fmax(0., rho*u - eps0 - epse - epsmu);
   // return total pressure
   return fmax(MIN_P, eps0*(GAM - 1.) + epsth*(GAMTH - 1.) + pe + pmu);
}

/// @brief electron affinity from hot equation of state
/// @details this function is called only if @ref FLUID == @ref PERFECT && @ref NSPECIES >= 2
/// @param[in] rho rest mass density
/// @param[in] ye electron fraction
/// @return affinity of the electron A_e = μn - μp - μe [c=G=M☉=kB=1]
double eos2affinity_e (const double rho, const double ye) {
   return -2.*MN*KE*(ye - YE0/RHON*rho);
}

/// @brief muon affinity from hot equation of state
/// @details this function is called only if @ref FLUID == @ref PERFECT && @ref NSPECIES == 3
/// @param[in] rho rest mass density
/// @param[in] ymu muon fraction
/// @return affinity of the muon A_μ = μn - μp - μμ [c=G=M☉=kB=1]
double eos2affinity_mu (const double rho, const double ymu) {
   return -2.*MN*KMU*(ymu - YMU0/RHON*rho);
}

/// @brief temperature at equilibrium
/// @param[in] rho rest mass density
/// @param[in] u internal specific (per unit mass) energy
/// @return equilibrium temperature T [c=G=M☉=kB=1]
double eos2teq (const double rho, const double u) {
   const double epsth= fmax(0., rho*u - K0*pow(rho, GAM));
   return 2.*MN*sqrt(fmax(0., KTH*epsth*pow(rho, GAMTH - 2.)));
}

/// @brief entropy at equilibrium
/// @param[in] rho rest mass density
/// @param[in] u internal specific (per unit mass) energy
/// @return equilibrium entropy per baryon s [c=G=M☉=kB=1]
double eos2seq (const double rho, const double u) {
   const double epsth= fmax(0., rho*u - K0*pow(rho, GAM));
   return sqrt(fmax(0., epsth/pow(rho, GAMTH)/KTH));
}

/// @brief get the net number reaction rates over the affinity
/// @param[in] rho rest mass density
/// @param[in] yeq particle fraction at equilibrium
/// @param[in] t temperature
/// @return net number reaction rate over affinity
double get_nrate_over_aff(const double rho, const double yeq, const double t) {
   // BEWARE: I am dividing μ and π by 10⁹K to avoid overflow of the doubles (KELVIN ~ 7e71)
   const double t0= M_PI*t/(1e9*KELVIN);
   const double tmp= t0*t0;
   const double nfactor= 8.86e31/(CM*CM*CM*SEC); // from cgs to code units
   const double cbrty= cbrt(yeq*rho/RHON);
   return nfactor*cbrty*17./30.*tmp*tmp/(1e9*KELVIN);
}

/// @brief get the total energy reaction rates
/// @param[in] rho rest mass density
/// @param[in] yeq particle fraction at equilibrium
/// @param[in] y particle fraction
/// @param[in] t temperature
/// @return total energy reaction rate
double get_erate (
   const double rho,
   const double yeq,
   const double y,
   const double t)
{
   // BEWARE: I am dividing μ and π by 10⁹K to avoid overflow of the doubles (KELVIN ~ 7e71)
   const double t0= M_PI*t/(1e9*KELVIN);
   const double tmp= t0*t0;
   const double efactor= 1.22e25/(CM*CM*CM*SEC)*ERG; // from cgs to code units
   const double cbrty= cbrt(yeq*rho/RHON);
   return efactor*cbrty*(1. + (y/yeq - 1.)/3.)*457./1260.*tmp*tmp*tmp;
}

/// @brief compute the inverse bulk viscous timescale
/// @param[in] rho rest mass density
/// @param[in] u internal specific (per unit mass) energy
/// @return 1/τ inverse bulk viscous timescale
double get_1_over_tau (const double rho, const double u) {
   if (NSPECIES == 1) return 1./TAU_BULK;
   const double rh= rho/RHON;
   const double t= eos2teq(rho, u);
   const double tmp= t*M_PI/(1e9*KELVIN);
   const double cons= 8.86e31/(CM*CM*CM*SEC)*cbrt(rh*rh)*17./30.*tmp*tmp*tmp*tmp/(1e9*KELVIN);
   const double num= KE*YE0*YE0 + (NSPECIES == 3 ? KMU*YMU0*YMU0 : 0.);
   const double den= pow(YE0, 5./3.) + (NSPECIES == 3 ? pow(YMU0, 5./3.) : 0.);
   return 2.*cons*MN*MN/rho*num/den;
}

/// @brief gravitational mass integrand
/// @param[in] r radius
/// @param[in] var conservative variables
/// @return gravitational mass integrand
double get_dmg (const double r, const double var[N0]) {
   const double tmp= var[2] + var[0];
   if (tmp > MIN_EPS)
      return 4.*M_PI*tmp*r*r;
   else return 0.;
}

/// @brief compute the gravitational mass
/// @param[in] x1 radial grid
/// @param[inout] var conservative variables (inout) and gravitational mass (out): [D, Sr, τɛ, ..., mg]
void get_mg (const double x1[D1], double var[D1][D0]) {
   double m= 0., fold= 0., xold= 0., f;
   for (int i= NG; i < NG+N1; i++) {
      f= get_dmg(x1[i], var[i]);
      m= m + 0.5*(x1[i] - xold)*(f + fold);
      fold= f;
      xold= x1[i];
      var[i][N0]= m;
   }
   // the boundary conditions will be set later in cons2prim_bc
   // note that conservative variables are not updated in the ghost cells
}

/// @brief baryonic mass integrand
/// @param[in] r radius
/// @param[in] var primitive variables
/// @return baryonic mass integrand
double get_dmb (const double r, const double var[D0]) {
   const double x= sqrt(r/(r - 2.*var[N0]));
   const double wv= var[1];
   const double rho= var[0];
   if (rho >= THR_RHO)
      return 4.*M_PI*rho*x*sqrt(1. + wv*wv)*r*r;
   else return 0.;
}

/// @brief compute the baryonic mass
/// @param[in] x1 radial grid
/// @param[inout] var primitive variables and gravitational mass: [ρ, Wv, u, ..., mg]
/// @return baryonic mass
double get_mb (const double x1[D1], const double var[D1][D0]) {
   double m= 0., fold= 0., xold= 0., f;
   for (int i= NG; i < NG+N1; i++) {
      f= get_dmb(x1[i], var[i]);
      m= m + 0.5*(x1[i] - xold)*(f + fold);
      fold= f;
      xold= x1[i];
   }
   return m;
}

/// @brief gravitational field potential integrand
/// @param[in] r radius
/// @param[in] var primitive variables
/// @return integrand dϕ
double get_dphi (const double r, const double var[D0]) {
   const double m= var[N0];
   const double rho= var[0];
   if (rho < THR_RHO) return m/r/(r - 2.*m);
   const double wv= var[1];
   const double u= var[2];
   const double ye= (FLUID == PERFECT && NSPECIES >= 2 ? var[3] : 0.);
   const double ymu= (FLUID == PERFECT && NSPECIES == 3 ? var[4] : 0.);
   const double pi= (FLUID == BULK ? var[3] : 0.);
   const double p= eos2p(rho, u, ye, ymu);
   const double srv= (rho*(1. + u) + p + pi)*wv*wv;
   return (m/r + 4.*M_PI*r*r*(srv + p + pi)) / (r - 2.*m);
}

/// @brief compute the metric field α
/// @details called only if @ref GRAVITATING == true
/// @param[in] x1 radial grid
/// @param[in] var primitive variables
/// @param[out] alp metric field α
void get_alpha (const double x1[D1], const double var[D1][D0], double alp[D1]) {
   double f, phi0= 0., xold= 0., fold= 0., phi[D1];
   for (int i= NG; i < D1; i++) {
      f= get_dphi(x1[i], var[i]);
      phi0= phi0 + 0.5*(x1[i] - xold)*(f + fold);
      xold= x1[i];
      fold= f;
      phi[i]= phi0;
   }
   // Dirichlet boundary conditions
   const int j= NG + N1 - 1;
   phi0= 0.5*log(1. - 2.*var[j][N0]/x1[j]) - phi0;
   /*
   // Newmann boundary conditions
   const int j= NG + N1 - 2;
   const double fp= exp(2.*phi[j+1]) - 1.;
   const double fm= exp(2.*phi[j-1]) - 1.;
   const double f0= -x1[j]*(fp - fm)/(x1[j+1] - x1[j-1]);
   phi0= 0.5*log(f0 + 1.) - phi[j];
   */
   // shift the potential and take case of the external boundary
   for (int i= NG; i < D1; i++) alp[i]= exp(phi[i] + phi0);
   // internal boundary condition
   for (int i= 0; i < NG; i++) alp[i]= alp[2*NG-i-1];
}

/// @brief set the boundary conditions
/// @param[inout] vec primitive variables
void set_boundaries (double vec[D1][D0]) {
   double dy;
   for (int i= 0; i < NG; i++) {
      for (int l= 0; l < D0; l++)
         // outflow
         if (!GRAVITATING) {
            vec[i][l]= vec[NG][l];
            vec[NG+N1+i][l]= vec[NG+N1-1][l];
         // spherical
         } else { // GRAVITATING
            // beware: the speed [l=1] is mirrored at the center!
            vec[i][l]= vec[2*NG-i-1][l] * (l == 1 ? -1. : 1.);
            // beware: this assume evenly spaced grid
            if (vec[NG+N1-1][0] <= MIN_RHO)
               dy= 0.0;
            else
               dy= vec[NG+N1-1][l] - vec[NG+N1-2][l];
            vec[NG+N1+i][l]= vec[NG+N1-1][l] + (i + 1)*dy;
         }
      // sanitize extrapolation
      if (GRAVITATING) {
         if (vec[NG+N1+i][0] < THR_RHO || vec[NG+N1+i][2] < MIN_U) {
            vec[NG+N1+i][0]= MIN_RHO;
            vec[NG+N1+i][1]= 0.;
            vec[NG+N1+i][2]= MIN_U;
            for (int l= 3; l < N0; l++) vec[NG+N1+i][l]= 0.;
         } else if (FLUID == PERFECT && NSPECIES >=2)
            for (int l= 3; l < N0; l++)
               vec[NG+N1+i][l]= fmax(0.0, fmin(1.0, vec[NG+N1+i][l]));
      }
   }
}

/// @brief get the primitives and the function to root for the final step
/// @param[in] wv velocity times Lorentz factor
/// @param[in] par modified conservative variables
/// @param[out] out primitive
/// @result function to root
double cons2prim1_final (const double wv, const double par[N0], double out[N0]) {
   const double wrho= par[0];
   const double sr= par[1];
   const double tau_plus_d= par[2];
   const double ye= (FLUID == PERFECT && NSPECIES >= 2 ? par[3] : 0.);
   const double ymu= (FLUID == PERFECT && NSPECIES == 3 ? par[4] : 0.);
   const double w= sqrt(1. + wv*wv);
   const double rho= fmax(MIN_RHO, wrho/w);
   const double eps= fmax(MIN_EPS, tau_plus_d - sr*wv/w);
   const double u= fmax(MIN_U, eps/rho - 1.);
   const double p= eos2p(rho, u, ye, ymu);
   out[0]= rho;
   out[1]= wv;
   out[2]= u;
   const double pi= (FLUID == BULK ?  fmax(-p, fmin(p, par[3]/w)) : 0.);
   if (FLUID == BULK) {
      out[3]= pi; // Π
   } else if (FLUID == PERFECT) {
      if (NSPECIES >= 2) out[3]= ye; // Ye
      if (NSPECIES == 3) out[4]= ymu; // Yμ
   }
   return sr - (eps + p + pi)*wv*w;
}

/// @brief invert conservative to primitive and set boundary conditions for the final step
/// @param[in] first_guess first guess of the primitive variables
/// @param[in] x1 spatial grid
/// @param[inout] var conservative (input) or primitive (output) variables
/// @param[in] time time
/// @return error code (0 on success)
int cons2prim_bc_final (
   const double first_guess,
   const double x1[D1],
   double var[D1][D0],
   const double time)
{
   double wv= first_guess, one_over_x= 1.;
   double par[N0];
   // compute the mass profile
   if (GRAVITATING) get_mg(x1, var);
   for (int i= NG; i < NG+N1; i++) {
      // prepare the conservartives
      if (PROBLEM == BJORKEN) one_over_x= 1./time;
      else if (GRAVITATING) one_over_x= sqrt((x1[i] - 2.*var[i][N0])/x1[i]);
      par[0]= fmax(MIN_RHO, var[i][0]*one_over_x); // Wρ ≥ ρmin
      par[1]= var[i][1]; // Sr
      par[2]= fmax(MIN_EPS, var[i][2] + var[i][0]); // τ + D ≥ εmin
      if (FLUID == BULK) {
         par[3]= var[i][3]*one_over_x; // ΠW
      } else if (FLUID == PERFECT) {
         if (NSPECIES >= 2) par[3]= fmax(0., fmin(1., var[i][3]/var[i][0])); // 0 ≤ Ye ≤ 1
         if (NSPECIES == 3) par[4]= fmax(0., fmin(1., var[i][4]/var[i][0])); // 0 ≤ Yμ ≤ 1
      }

      // invert from conservatives to primitives
      if (par[0] < THR_RHO || par[2] < MIN_EPS) {
         var[i][0]= MIN_RHO;
      } else if (FLUID == DUST) { // ε = ρ, p = 0
         wv= par[1]/par[0];
         var[i][0]= fmax(MIN_RHO, par[0]/sqrt(1. + wv*wv)); // ρ
         var[i][1]= wv; // Wv
         var[i][2]= 0.; // u
      } else if (PROBLEM == BJORKEN) {
         cons2prim1_final(0., par, var[i]);
      } else {
         if (find_root(-HUGE_VAL, HUGE_VAL, &cons2prim1_final, par, var[i], wv))
            return error("cons2prim_bc");
      }

      // set atmosphere
      if (var[i][0] < THR_RHO) {
         var[i][0]= MIN_RHO;
         var[i][1]= 0.;
         var[i][2]= MIN_U;
         for (int l= 3; l < N0; l++) var[i][l]= 0.;
      }
   }
   set_boundaries(var);
   return 0;
}

/// @brief get the primitives and the function to root for the implicit step
/// @param[in] wv velocity times Lorentz factor
/// @param[in] par modified conservative variables
/// @param[out] out primitive
/// @result function to root
double cons2prim1_implicit (const double wv, const double par[N0+5], double out[N0]) {
   const double wrho= par[0];
   const double sr= par[1];
   const double tau_plus_d= par[2];
   const double w= sqrt(1. + wv*wv);
   const double rho= fmax(MIN_RHO, wrho/w);
   const double eps= fmax(MIN_EPS, tau_plus_d - sr*wv/w);
   const double u= fmax(MIN_U, eps/rho - 1.);
   const double alp_dt= par[N0+1]; // α=1 for BJORKEN -> alp_dt= dt
   out[0]= rho;
   out[1]= wv;
   out[2]= u;
   double ye= 0., ymu= 0.;
   if (FLUID == PERFECT && NSPECIES >= 2) {
      const double t= eos2teq(rho, u);
      const double yeqe= YE0/RHON*rho;
      const double facte= 2.*get_nrate_over_aff(rho, yeqe, t)*MN*MN*alp_dt*KE;
      // implicit integration
      ye= fmax(0., fmin(1., (par[3] + facte*yeqe)/(wrho + facte))); // 0 ≤ Ye ≤ 1
      out[3]= ye; // Ye
      if (NSPECIES == 3) {
         const double yeqmu= YMU0/RHON*rho;
         const double factmu= 2.*get_nrate_over_aff(rho, yeqmu, t)*MN*MN*alp_dt*KMU;
         // implicit integration
         ymu= fmax(0., fmin(1., (par[4] + factmu*yeqmu)/(wrho + factmu))); // 0 ≤ Yμ ≤ 1
         out[4]= ymu; // Yμ
      }
   }
   const double p= eos2p(rho, u, ye, ymu);
   double pi= 0.;
   if (FLUID == BULK) {
      const double one_over_chi= get_1_over_chi(rho);
      const double one_over_tau= get_1_over_tau(rho, u);
      const double log_chi_over_t= -log(one_over_chi*fmax(eos2teq(rho, u), TINY_T));
      // implicit integration
      if (PROBLEM == BJORKEN) {
         if (HISCOCK_LINDBLOM)
            pi= 2.*par[3]/(2. + alp_dt*(2.*one_over_tau - 1./par[N0])
                          + par[N0+2]*log_chi_over_t - par[N0+3]);
         else
            pi= par[3]/(1. + alp_dt*(one_over_tau - 1./par[N0]));
      } else {
         if (HISCOCK_LINDBLOM) {
            pi= 2.0*(par[3] - one_over_chi*(w*par[N0+2] - par[N0]))/
                (w*(2. - par[N0+2] - par[N0+3] + par[N0+2]*log_chi_over_t)
                + 2.*alp_dt*one_over_tau + par[N0] + wv*par[N0+4]);
         } else {
            pi= (par[3] - one_over_chi*(w*par[N0+2] - par[N0]))/
                (w*(1. - par[N0+2]) + alp_dt*one_over_tau + par[N0]);
         }
      }
      pi= fmax(-p, fmin(p, pi));
      out[3]= pi; // Π
   }
   return sr - (eps + p + pi)*wv*w;
}

/// @brief invert conservative to primitive and set boundary conditions for the implicit step
/// @param[in] first_guess first guess of W·v
/// @param[in] x1 spatial grid
/// @param[inout] var conservative (input) or primitive (output) variables
/// @param[in] xw0 X·W at time t0
/// @param[in] alp metric function α
/// @param[in] der 1/r²·d(r²αWv)/dr
/// @param[in] ln0 log(χ/T) at time t0
/// @param[in] adl α·d(log(χ/T))/dr
/// @param[in] time time
/// @param[in] coef_dt timestep times Butcher coefficient
/// @param[in] ratio ratio of the Butcher coefficients
/// @return error code (0 on success)
int cons2prim_bc_implicit (
   const double first_guess,
   const double x1[D1],
   double var[D1][D0],
   const double xw0[D1],
   const double alp[D1],
   const double der[D1],
   const double ln0[D1],
   const double adl[D1],
   const double time,
   const double coef_dt,
   const double ratio
) {
   double wv= first_guess, one_over_x= 1.;
   double par[N0+5];
   // compute the mass profile
   if (GRAVITATING) get_mg(x1, var);
   for (int i= NG; i < NG+N1; i++) {
      // prepare the conservartives
      if (PROBLEM == BJORKEN) one_over_x= 1./time;
      else if (GRAVITATING) one_over_x= sqrt((x1[i] - 2.*var[i][N0])/x1[i]);
      par[0]= fmax(MIN_RHO, var[i][0]*one_over_x); // Wρ ≥ ρmin
      par[1]= var[i][1]; // Sr
      par[2]= fmax(MIN_EPS, var[i][2] + var[i][0]); // τ + D ≥ εmin
      if (FLUID == BULK) {
         par[3]= var[i][3]*one_over_x; // ΠW
         if (PROBLEM == BJORKEN) {
            par[N0]= time;
            par[N0+1]= coef_dt;
            par[N0+2]= ratio;
            if (HISCOCK_LINDBLOM) par[N0+3]= ratio*ln0[i];
         } else {
            par[N0]= (xw0[i]*ratio - der[i]*coef_dt)*one_over_x;
            par[N0+1]= alp[i]*coef_dt;
            par[N0+2]= ratio;
            if (HISCOCK_LINDBLOM) par[N0+3]= ratio*ln0[i];
            if (HISCOCK_LINDBLOM) par[N0+4]= coef_dt*adl[i]*one_over_x;
         }
      } else if (FLUID == PERFECT && NSPECIES >= 2) {
         par[3]= var[i][3]*one_over_x; // WρYe
         if (NSPECIES == 3) par[4]= var[i][4]*one_over_x; // WρYμ
         par[N0+1]= alp[i]*coef_dt; // α(t)dt
      }

      // invert from conservatives to primitives
      if (par[0] < THR_RHO || par[2] < MIN_EPS) {
         var[i][0]= MIN_RHO;
      } else if (FLUID == DUST) { // ε = ρ, p = 0
         wv= par[1]/par[0];
         var[i][0]= fmax(MIN_RHO, par[0]/sqrt(1. + wv*wv)); // ρ
         var[i][1]= wv; // Wv
         var[i][2]= 0.; // u
      } else if (PROBLEM == BJORKEN) {
         cons2prim1_implicit(0., par, var[i]);
      } else {
         if (find_root(-HUGE_VAL, HUGE_VAL, &cons2prim1_implicit, par, var[i], wv))
            return error("cons2prim_bc_implicit");
      }

      // set atmosphere
      if (var[i][0] < THR_RHO) {
         var[i][0]= MIN_RHO;
         var[i][1]= 0.;
         var[i][2]= MIN_U;
         for (int l= 3; l < N0; l++) var[i][l]= 0.;
      }
   }
   set_boundaries(var);
   return 0;
}

/// @brief primitive to conservative inversion for a single point
/// @param[inout] var primitive (input) or conservative (output) variables
/// @param[in] xt metric function X (for @ref GRAVITATING) or time t (for @ref BJORKEN)
void prim2cons1(double var[N0], const double xt) {
   const double rho= var[0];
   const double u= var[2];
   const double ye= (FLUID == PERFECT && NSPECIES >= 2 ? var[3] : 0.);
   const double ymu= (FLUID == PERFECT && NSPECIES == 3 ? var[4] : 0.);
   const double pi= (FLUID == BULK ? var[3] : 0.);
   const double wv= (PROBLEM == BJORKEN ? 0. : var[1]);
   const double w= (PROBLEM == BJORKEN ? 1. : sqrt(1. + wv*wv));
   const double wx= w*(GRAVITATING || (PROBLEM == BJORKEN) ? xt : 1.);
   const double d= rho*wx;
   const double p= eos2p(rho, u, ye, ymu);
   const double hrho= rho*(1. + u) + p;
   var[0]= d; // D
   var[1]= (hrho + pi)*wv*w; // Sr
   var[2]= (hrho + pi)*w*w - (p + pi) - d; // τɛ
   if (FLUID == BULK) {
      var[3]= wx*pi; // XWΠ
   } else if (FLUID == PERFECT) {
      if (NSPECIES >= 2) var[3]= d*ye; // DYe
      if (NSPECIES == 3) var[4]= d*ymu; // DYμ
   }
}

/// @brief second order limiter of the first order derivative
/// @param[in] xm position x0 - dx
/// @param[in] x0 position x0
/// @param[in] xp position x0 + dx
/// @param[in] cm quantity at x0 - dx
/// @param[in] c0 quantity at x0
/// @param[in] cp quantity at x0 + dx
/// @return slope limited derivative
double slope_limit_der (
   const double xm,
   const double x0,
   const double xp,
   const double cm,
   const double c0,
   const double cp
) {
   const double d= (cp - cm) / (xp - xm);
   if (SLOPE_LIMITER == NONE) return d;
   else {
      const double dm= (c0 - cm) / (x0 - xm);
      const double dp= (cp - c0) / (xp - x0);
      // flux limiter
      if (SLOPE_LIMITER == MINMOD)
         return dm*dp > 0.
                ? (dm*dm < dp*dp ? dm : dp)
                : 0.;
      else // SLOPE_LIMITER == MCLIM
         return dm*dp > 0.
                ? copysign(fmin(fabs(d), 2.*fmin(fabs(dm), fabs(dp))), d)
                : 0.;
   }
}

/// @brief reconstruct the left and right states
/// @param[in] i grid index
/// @param[in] x1 spatial grid
/// @param[in] var variables to reconstruct
/// @param[out] left reconstructed left interface
/// @param[out] right reconstructed right interface
void reconstruct (
   const int i,
   const double x1[D1],
   const double var[D1][D0],
   double left[N0],
   double right[N0])
{
   // reconstruct
   double d;
   for (int l= 0; l < N0; l++) {
      if (ORDER_SPACE == 1) {
         right[l]= var[i][l];
         left[l]=  var[i][l];
      } else { // ORDER_SPACE == 2
         d= slope_limit_der(x1[i-1], x1[i], x1[i+1], var[i-1][l], var[i][l], var[i+1][l]);
         right[l]= var[i][l] - 0.5*d*(x1[i] - x1[i-1]);
         left[l]=  var[i][l] + 0.5*d*(x1[i+1] - x1[i]);
      }
   }

   // sanitaze the output
   // rest mass density
   if (right[0] < THR_RHO || right[2] < MIN_U) {
      right[0]= MIN_RHO;
      right[1]= 0.;
      right[2]= MIN_U;
      for (int l= 3; l < N0; l++) right[l]= 0.;
   }
   // internal specific energy
   if (left[0] < THR_RHO || left[2] < MIN_U) {
      left[0]= MIN_RHO;
      left[1]= 0.;
      left[2]= MIN_U;
      for (int l= 3; l < N0; l++) left[l]= 0.;
   }
   // number fractions
   if (FLUID == PERFECT && NSPECIES >= 2) {
      right[3]= fmax(0., fmin(1., right[3]));
      left[3]= fmax(0., fmin(1., left[3]));
      if (NSPECIES == 3) {
         right[4]= fmax(0., fmin(1., right[4]));
         left[4]= fmax(0., fmin(1., left[4]));
      }
   // bulk viscous stress
   } else if (FLUID == BULK) {
      double p= eos2p(right[0], right[2], 0., 0.);
      right[3]= fmax(-p, fmin(p, right[3]));
      p= eos2p(left[0], left[2], 0., 0.);
      left[3]= fmax(-p, fmin(p, left[3]));
   } // end if (FLUID)
} // end reconstruct

/// @brief get the implicit source
/// @param[in] x1 radius
/// @param[in] var primitive variables
/// @param[in] alp metric field α
/// @param[in] xw0 X·W at time t0
/// @param[in] ln0 log(χ/T) at time t0
/// @param[in] t0 initial time of the timestep
/// @param[in] dt time step
/// @param[out] der implicit RHS for the time evolution
void get_der_implicit (
   const double x1[D1],
   const double var[D1][D0],
   const double alp[D1],
   const double xw0[D1],
   const double ln0[D1],
   const double t0,
   const double dt,
   double der[D1][N0])
{
   // init the derivative
   for (int i= NG; i < NG+N1; i++)
      for (int l= 0; l < N0; l++)
         der[i][l]= 0.;

   if (FLUID == DUST) return;

   double pi= 0., x= 1., alp0= 1.;
   double rho, u, t, w, cm, c0, cp, der0, adl, lnt[D1];
   const double time= t0 + dt;
   const double fact= (HISCOCK_LINDBLOM ? 0.5 : 1.);

   // init log(chi/T)
   if (FLUID == BULK && HISCOCK_LINDBLOM)
      for (int i= NG-1; i < NG+N1+1; i++)
         lnt[i]= -log(get_1_over_chi(var[i][0])*fmax(eos2teq(var[i][0], var[i][2]), TINY_T));

   for (int i= NG; i < NG+N1; i++) {
      rho= var[i][0];
      if (rho < THR_RHO) return;
      if (GRAVITATING) {
         x= sqrt(x1[i]/(x1[i] - 2.*var[i][N0]));
         alp0= alp[i];
      }
      // else: alp0=x=1
      // primitives have been already sanitized: if u < threeshold, also ρ < threeshold
      u= var[i][2];
      if (FLUID == BULK) {
         pi= var[i][3];
         if (PROBLEM == BJORKEN) {
            der[i][3]= -pi*(time*get_1_over_tau(rho, u) - fact);
            if (HISCOCK_LINDBLOM)
               der[i][3]= der[i][3] - 0.5*pi*time*(lnt[i] - ln0[i])/dt;
         } else {
            cm= x1[i-1]*x1[i-1]*alp[i-1]*var[i-1][1]; // r²αWv-
            c0= x1[i]*x1[i]*alp[i]*var[i][1];         // r²αWv0
            cp= x1[i+1]*x1[i+1]*alp[i+1]*var[i+1][1]; // r²αWv+
            der0= slope_limit_der(x1[i-1], x1[i], x1[i+1], cm, c0, cp)/(x1[i]*x1[i]);
            w= sqrt(1. + var[i][1]*var[i][1]);
            // BEWARE: this assumes evenly spaced spatial grid
            der[i][3]= -alp0*x*pi*get_1_over_tau(rho, u) \
                       -(get_1_over_chi(rho) - fact*pi)*((x*w - xw0[i])/dt + der0);
            if (HISCOCK_LINDBLOM) {
               adl= alp0*slope_limit_der(x1[i-1], x1[i], x1[i+1], lnt[i-1], lnt[i], lnt[i+1]);
               der[i][3]= der[i][3] - 0.5*pi*(w*x*(lnt[i] - ln0[i])/dt + var[i][1]*adl);
            }
         }
      } else if (FLUID == PERFECT && NSPECIES >= 2) {
         t= eos2teq(rho, u);
         der[i][3]= alp0*x*get_nrate_over_aff(rho, rho*YE0/RHON, t)*eos2affinity_e(rho, var[i][3])*MN;
         if (NSPECIES == 3)
            der[i][4]= alp0*x*get_nrate_over_aff(rho, rho*YMU0/RHON, t)*eos2affinity_mu(rho, var[i][4])*MN;
      } // end if (FLUID)
   } // end for
}

/// @brief compute the flux
/// @note flux[i]= flux_{i-1/2}
/// @param[in] var primitive variables
/// @param[out] flux flux
/// @param[in] x metric function
void get_flux(const double var[N0], double flux[N0], const double x) {
   const double rho= var[0];
   const double wv= var[1];
   const double u= var[2];
   const double pi= (FLUID == BULK ? var[3] : 0.);
   const double ye= (FLUID == PERFECT && NSPECIES >= 2 ? var[3] : 0.);
   const double ymu= (FLUID == PERFECT && NSPECIES == 3 ? var[4] : 0.);
   const double p= eos2p(rho, u, ye, ymu);
   const double hrho= rho*(1. + u) + p;
   const double rhox= rho*(GRAVITATING ? x : 1.);
   flux[0]= rhox*wv;
   flux[1]= (hrho + pi)*wv*wv + p + pi;
   flux[2]= (hrho + pi)*wv*sqrt(1. + wv*wv) - rhox*wv;
   if (FLUID == BULK) {
      flux[3]= wv*x*pi;
   } else if (FLUID == PERFECT) {
      if (NSPECIES >= 2) flux[3]= rhox*wv*ye;
      if (NSPECIES == 3) flux[4]= rhox*wv*ymu;
   }
}

/// @brief solve the Riemann problem on the interface
/// @details called only if @ref PROBLEM != @ref BJORKEN
/// @param[in] left primitive variables on the left interface
/// @param[in] right primitive variables on the right interface
/// @param[in] x metric field X at the interface
/// @param[out] flux the flux at the interface
void solve_riemann (
   const double left[D0],
   const double right[D0],
   const double x,
   double flux[N0])
{
   double fl[N0], fr[N0], consl[D0], consr[D0];
   double vl, vr, cl, cr, vm, vp;
   const double yel= (FLUID == PERFECT && NSPECIES >= 2 ? left[3] : 0.);
   const double ymul= (FLUID == PERFECT && NSPECIES == 3 ? left[4] : 0.);
   const double yer= (FLUID == PERFECT && NSPECIES >= 2 ? right[3] : 0.);
   const double ymur= (FLUID == PERFECT && NSPECIES == 3 ? right[4] : 0.);

   // compute maximal and minimal signal speeds
   cl= eos2c(left[0], left[2], yel, ymul);
   cr= eos2c(right[0], right[2], yer, ymur);
   vl= left[1]/sqrt(1. + left[1]*left[1]); // v= Wv/W
   vr= right[1]/sqrt(1. + right[1]*right[1]); // v= Wv/W
   vm= fmin(0., fmin((vl - cl)/(1. - vl*cl), (vr - cr)/(1. - vr*cr)));
   vp= fmax(0., fmax((vl + cl)/(1. + vl*cl), (vr + cr)/(1. + vr*cr)));

   // compute the fluxes
   get_flux(left,  fl, x);
   get_flux(right, fr, x);

   // convert the primitive variables to conservatives
   for (int l= 0; l < D0; l++) {
      consl[l]= left[l];
      consr[l]= right[l];
   }
   prim2cons1(consl, x);
   prim2cons1(consr, x);

   // solve the Riemann problem
   for (int l= 0; l < N0; l++) {
      if (USE_LLF_RIEMANN_SOLVER)
         flux[l]= 0.5*(fl[l] + fr[l] - fmax(fabs(vm), fabs(vp))*(consr[l] - consl[l]));
      else if (vm >= 0.)
         flux[l]= fl[l];
      else if (vp <= 0.)
         flux[l]= fr[l];
      else
         flux[l]= (vp*fl[l] - vm*fr[l] + vp*vm*(consr[l] - consl[l]))/(vp - vm);
   }
}

/// @brief compute the explicit source and the fluxes
/// @param[in] x1 spatial grid
/// @param[in] var primitive variables
/// @param[in] alp metric field α
/// @param[in] time time
/// @param[out] der explicit RHS for the time evolution
void get_der_explicit (
   const double x1[D1],
   const double var[D1][D0],
   const double alp[D1],
   const double time,
   double der[D1][N0])
{
   double x= 1., r, m;

   // init the derivative
   for (int i= NG; i < NG+N1; i++)
      for (int l= 0; l < N0; l++)
         der[i][l]= 0.;

   // compute the source of the momentum equation
   if (GRAVITATING) {
      double rho, u, r2, p, wv, teq, yeq, y;
      double pi= 0., ye= 0., ymu= 0., q= 0.;
      for (int i= NG; i < NG+N1; i++) {
         rho= var[i][0];
         if (rho < THR_RHO) break;
         wv= var[i][1];
         // primitives have been already sanitized: if u < threeshold, also ρ < threeshold
         u= var[i][2];
         if (FLUID == BULK) {
            pi= var[i][3];
            if (INCLUDE_LUMINOSITY && NSPECIES == 2) {
               yeq= rho*YE0/RHON;
               y= fmax(0., fmin(1., yeq - 0.5*pi/(rho*KE*yeq)));
               q= get_erate(rho, yeq, y, eos2teq(rho, u));
            }
         } else if (FLUID == PERFECT && NSPECIES >= 2) {
            ye= var[i][3];
            if (INCLUDE_LUMINOSITY) {
               teq= eos2teq(rho, u);
               q= get_erate(rho, YE0/RHON*rho, ye, teq);
            }
            if (NSPECIES == 3) {
               ymu= var[i][4];
               if (INCLUDE_LUMINOSITY) q= q + get_erate(rho, YMU0/RHON*rho, ymu, teq);
            }
         }
         m= var[i][N0];
         r= x1[i];
         r2= r*r;
         p= eos2p(rho, u, ye, ymu);
         x= sqrt(r/(r - 2.*m));
         der[i][1]= alp[i]*(x*(-rho*(1. + u)*(8.*M_PI*r*(p + pi) + m/r2) \
            + (p + pi)*m/r2) + 2.*(p + pi)/x/r);
         if (INCLUDE_LUMINOSITY && NSPECIES >= 2) {
            der[i][1]= der[i][1] - alp[i]*wv*q;
            der[i][2]= -alp[i]*sqrt(1. + wv*wv)*q;
         }
      }
   } else if (PROBLEM == BJORKEN) {
      double rho, u, p, pi;
      for (int i= NG; i < NG+N1; i++) {
         rho= var[i][0];
         u= var[i][2];
         pi= var[i][3];
         p= eos2p(rho, u, 0., 0.);
         der[i][2]= -(rho*(1. + u) + p + pi)/time;
         der[i][3]= -get_1_over_chi(rho);
      }
      // the Bjorken flow does not have fluxes
      return;
   }

   double tmp[D0], left[D0], right[D0], flux[D0], fold[D0]; // flux_i= flux_{i-1/2}
   double den, oldfactor, factor= 1.;

   // start from one point left and stop one point right
   for (int i= NG-1; i < N1+NG+1; i++) {
      // load the left state, save the old flux
      for (int l= 0; l < D0; l++) {
         left[l]= tmp[l];
         fold[l]= flux[l];
         oldfactor= factor;
      }

      // compute the current right state and the next left state (tmp)
      reconstruct(i, x1, var, tmp, right);

      // we don't need to compute the derivative for ghost cells
      if (i < NG) continue;

      // update the mass and potential at half grid
      if (GRAVITATING) {
         if (i == NG) {
            x= 1.;
            factor= 0.;
         } else {
            r= 0.5*(x1[i] + x1[i-1]);
            m= 0.5*(var[i][N0] + var[i-1][N0]);
            x= sqrt(r/(r - 2.*m));
            // BEWARE: the form in which you write it greatly change the result!
            // alp[i-1/2] = exp(phi[i-1/2])= sqrt(alp[i]*alp[i-1])
            factor= r*sqrt(alp[i]*alp[i-1]*r*(r - 2.*m));
         }
      }
      // solve the Riemann problem to determine the flux
      solve_riemann(left, right, x, flux);

      // we can compute the derivative once we know the right flux
      if (i < NG+1) continue;

      // compute the derivative
      for (int l= 0; l < N0; l++) {
         // radial coordinate
         if (GRAVITATING) {
            den= 0.5*(x1[i] - x1[i-2])*x1[i-1]*x1[i-1];
            der[i-1][l]= der[i-1][l] - (factor*flux[l] - oldfactor*fold[l])/den;
         // cartesian coordinate
         } else {
            den= 0.5*(x1[i] - x1[i-2]);
            der[i-1][l]= der[i-1][l] - (flux[l] - fold[l])/den;
         }
      }
      // null velocity in the atmosphere
      if (var[i-1][0] < THR_RHO && GRAVITATING) der[i-1][1]= 0.;
   }
}

/// @brief evolve the equations for one time step with IMplicit-Explicit integration
/// @param[in] t time
/// @param[in] dt timestep
/// @param[in] x1 spatial grid
/// @param[inout] var primitive variables at t (input) and t + dt (output)
/// @return error code (0 on success)
int imex_integration (
   const double t,
   const double dt,
   const double x1[D1],
   double var[D1][D0])
{
   const double first_guess= var[NG][1]; // Wv
   const double g= 1. - sqrt(0.5); // γ
   double der[D1][N0], xw0[D1], alp[D1], der0[D1], ln0[D1], adl[D1], x= 1.;

   get_alpha(x1, var, alp);
   get_der_explicit(x1, var, alp, t, der);
   // preparation
   if (FLUID == BULK && HISCOCK_LINDBLOM)
      for (int i= NG-1; i < NG+N1+1; i++)
         ln0[i]= -log(get_1_over_chi(var[i][0])*fmax(eos2teq(var[i][0], var[i][2]), TINY_T));
   if (GRAVITATING && FLUID == BULK) {
      double cm, c0, cp;
      for (int i= NG; i < NG+N1; i++) {
         cm= x1[i-1]*x1[i-1]*alp[i-1]*var[i-1][1];
         c0= x1[i]*x1[i]*alp[i]*var[i][1];
         cp= x1[i+1]*x1[i+1]*alp[i+1]*var[i+1][1];
         der0[i]= slope_limit_der(x1[i-1], x1[i], x1[i+1], cm, c0, cp)/(x1[i]*x1[i]);
         if (HISCOCK_LINDBLOM)
            adl[i]= alp[i]*slope_limit_der(x1[i-1], x1[i], x1[i+1], ln0[i-1], ln0[i], ln0[i+1]);
      }
   }
   // primitive to conservative
   for (int i= NG; i < NG+N1; i++) {
      if (GRAVITATING) x= sqrt(x1[i]/(x1[i] - 2.*var[i][N0]));
      if (GRAVITATING && FLUID == BULK) xw0[i]= x*sqrt(1. + var[i][1]*var[i][1]);
      prim2cons1(var[i], (PROBLEM == BJORKEN ? t : x));
   }
   // now var is conservative at t

   // semi-implicit Euler
   // 1st:  0 | 0        0 | 0
   // 2nd:  0 | 0  0     1 | 0  1
   //       --|-----     --|-----
   // last:   | 1  0       | 0  1
   if (ORDER_TIME == 1) {
      for (int i= NG; i < N1+NG; i++)
         for (int l= 0; l < N0; l++)
            var[i][l]= var[i][l] + dt*der[i][l];
      if (cons2prim_bc_implicit(first_guess, x1, var, xw0, alp, der0, ln0, adl, t+dt, dt, 1.))
         return error("imex_integration: 1st order, last step");
   }
   // 2nd order IMEX RK
   // 1st:  0 | 0            γ   | γ
   // 2nd:  1 | 1    0       1-γ | 1-2γ  γ
   //       --|---------     ----|---------
   // last:   | 1/2  1/2         | 1/2   1/2
   else if (ORDER_TIME == 2) {
      double d1e[D1][N0], d2e[D1][N0], d1i[D1][N0], d2i[D1][N0], tmp[D1][D0];
      // 1st step
      for (int i= NG; i < N1+NG; i++)
         for (int l= 0; l < N0; l++)
            tmp[i][l]= var[i][l];
      if (cons2prim_bc_implicit(first_guess, x1, tmp, xw0, alp, der0, ln0, adl, t+g*dt, g*dt, 1.))
         return error("imex_integration: 2nd order, 1st step");
      get_der_explicit(x1, tmp, alp, t, d1e);
      get_der_implicit(x1, tmp, alp, xw0, ln0, t, g*dt, d1i);
      // 2nd step
      for (int i= NG; i < N1+NG; i++)
         for (int l= 0; l < N0; l++)
            tmp[i][l]= var[i][l] + dt*(d1e[i][l] + (1. - 2.*g)*d1i[i][l]);
      if (cons2prim_bc_implicit(first_guess, x1, tmp, xw0, alp, der0, ln0, adl, t+(1.-g)*dt, g*dt, g/(1.-g)))
         return error("imex_integration: 2nd order, 2nd step");
      get_der_explicit(x1, tmp, alp, t + dt, d2e);
      get_der_implicit(x1, tmp, alp, xw0, ln0, t, (1. - g)*dt, d2i);
      // last step
      for (int i= NG; i < N1+NG; i++)
         for (int l= 0; l < N0; l++)
            var[i][l]= var[i][l] + 0.5*dt*(d1e[i][l] + d2e[i][l] + d1i[i][l] + d2i[i][l]);
      if (cons2prim_bc_final(first_guess, x1, var, t + dt))
         return error("imex_integration: 2nd order, last step");
   }
   // 3nd order IMEX RK
   // 1st:  0   | 0                  γ   | γ
   // 2nd:  1   | 1    0             1-γ | 1-2γ   γ
   // 3rd:  1/2 | 1/4  1/4  0        1/2 | 1/2-γ  0    γ
   //       ----|--------------      ----|----------------
   // last:     | 1/6  1/6  2/3          | 1/6    1/6  2/3
   else {
      double d1e[D1][N0], d2e[D1][N0], d3e[N1][N0];
      double d1i[D1][N0], d2i[D1][N0], d3i[N1][N0];
      double tmp[D1][D0];
      // 1st step
      for (int i= NG; i < N1+NG; i++)
         for (int l= 0; l < N0; l++)
            tmp[i][l]= var[i][l];
      if (cons2prim_bc_implicit(first_guess, x1, tmp, xw0, alp, der0, ln0, adl, t+g*dt, g*dt, 1.))
         return error("imex_integration: 3rd order, 1st step");
      get_der_explicit(x1, tmp, alp, t, d1e);
      get_der_implicit(x1, tmp, alp, xw0, ln0, t, g*dt, d1i);
      // 2nd step
      for (int i= NG; i < N1+NG; i++)
         for (int l= 0; l < N0; l++)
            tmp[i][l]= var[i][l] + dt*(d1e[i][l] + (1. - 2.*g)*d1i[i][l]);
      if (cons2prim_bc_implicit(first_guess, x1, tmp, xw0, alp, der0, ln0, adl, t+(1.-g)*dt, g*dt, g/(1.-g)))
         return error("imex_integration: 3rd order, 2nd step");
      get_der_explicit(x1, tmp, alp, t + dt, d2e);
      get_der_implicit(x1, tmp, alp, xw0, ln0, t, (1. - g)*dt, d2i);
      // 3rd step
      for (int i= NG; i < N1+NG; i++)
         for (int l= 0; l < N0; l++)
            tmp[i][l]= var[i][l] + 0.25*dt*(d1e[i][l] + d2e[i][l] + (1. - 2.*g)*d1i[i][l]);
      if (cons2prim_bc_implicit(first_guess, x1, tmp, xw0, alp, der0, ln0, adl, t+0.5*dt, g*dt, 2.*g))
         return error("imex_integration: 3rd order, 3rd step");
      get_der_explicit(x1, tmp, alp, t + 0.5*dt, d3e);
      get_der_implicit(x1, tmp, alp, xw0, ln0, t, 0.5*dt, d3i);
      // last step
      for (int i= NG; i < N1+NG; i++)
         for (int l= 0; l < N0; l++)
            var[i][l]= var[i][l] + dt/6.*(d1e[i][l] + d1i[i][l] + d2e[i][l] + d2i[i][l] + 4.*(d3e[i][l] + d3i[i][l]));
      if (cons2prim_bc_final(first_guess, x1, var, t + dt))
         return error("imex_integration: 3rd order, last step");
   }
   return 0;
}

/// @brief compute the timestep with the Courant criterion
/// @param[in] x1 spatial grid
/// @param[in] var primitive variables
/// @return on success timestep > 0., on failure 0.
double timestep (const double x1[D1], const double var[D1][D0]) {
   if (PROBLEM == CLOUD && FLUID == DUST) {
      return CFL*fabs(x1[NG] - x1[NG-1]);
      /*
      int isurf= N1+NG;
      for (int i= NG; i < N1+NG; i++)
         if (var[i][0] < THR_RHO) {
            isurf= i-1;
            break;
         }
      const double wv= var[isurf][1];
      v= wv/sqrt(1. + wv*wv);
      const double rho= var[isurf][0];
      const double r= x1[isurf];
      const double m= var[isurf][N0];
      const double x= sqrt(r/(r - 2.*m));
      // formula from first term of Eq. (2) of Thorne (1966), ApJ 144:201.
      const double a= -x*x*x*(m/(r*r) + 4.*M_PI*r*rho*wv*wv); // negative
      const double dx= x1[isurf-1] - x1[isurf]; // negative
      const double dt= -(v + sqrt(v*v + 2.*a*dx))/a; // discriminant is positive 
      return CFL*dt;
      */
   } else if (PROBLEM == BJORKEN) {
      return 0.1*CFL;
   }
   double tmp= 0., ye= 0., ymu= 0.;
   double v, avm, avp, c;
   for (int i= NG; i < N1+NG; i++) {
      if (FLUID == PERFECT) {
         if (NSPECIES >= 2) ye= var[i][3];
         if (NSPECIES == 3) ymu= var[i][4];
      }
      c= eos2c(var[i][0], var[i][2], ye, ymu);
      v= var[i][1]/sqrt(1. + var[i][1]*var[i][1]); // v= Wv/W
      avp= fabs((v + c)/(1. + v*c));
      avm= fabs((v - c)/(1. - v*c));
      tmp= fmax(tmp, 2.*fmax(fabs(v), fmax(avm, avp))/fabs(x1[i+1] - x1[i-1]));
   }
   if (tmp <= 0.) return tmp; // propagate error
   else return CFL/tmp;
}

/// @brief initialize the simulation
/// @param[out] x1 spatial grid
/// @param[out] var primitive variables
void initialize(double x1[D1], double var[D1][D0]) {
   // init the grid
   for (int i= 0; i < D1; i++) x1[i]= (double) (i-NG)/(double) (N1 - 1);
   if (GRAVITATING) {
      // staggered grid, avoids singularity at r=0
      const double shift= 0.5*(x1[NG] - x1[NG-1]);
      for (int i= 0; i < D1; i++) x1[i]= (x1[i] + shift)*MAX_R;
   }

   if (PROBLEM == BJORKEN) {
      for (int i= 0; i < D1; i++) {
         var[i][0]= RHO0; // intial rest mass density
         var[i][1]= 0.; // null velocity
         var[i][2]= rho2eps(RHO0, S0)/RHO0 - 1.; // initial internal energy
         var[i][3]= PI0; // Π
      }
   } else if (PROBLEM == FLUXTUBE) {
      // Sod's shock tube
      for (int i= 0; i < D1; i++) {
         if (i < D1/2) {
            var[i][0]= RHOL; // density
            var[i][1]= VL/sqrt(1. - VL*VL); // velocity*Lorentz
            var[i][2]= PL/(GAM - 1.)/RHOL; // internal energy
         } else {
            var[i][0]= RHOR; // density
            var[i][1]= VR/sqrt(1. - VR*VR); // velocity*Lorentz
            var[i][2]= PR/(GAM - 1.)/RHOR; // internal energy
         }
      }
   } else if (PROBLEM == CLOUD) {
      double m= 0.;
      for (int i= NG; i < N1+NG; i++) {
         if (x1[i] <= R0) {
            var[i][0]= RHO0; // constant density
            m= 4./3.*M_PI*RHO0*x1[i]*x1[i]*x1[i]; // gravitational mass
            var[i][2]= rho2eps(RHO0, S0)/RHO0 - 1.; // internal energy
         } else {
            var[i][0]= MIN_RHO;
            var[i][2]= MIN_U; // internal energy
         }
         var[i][1]= 0.; // null velocity (Wv)
         var[i][N0]= m;
      }
   } else { // PROBLEM == STAR
      // we start from inside and integrate outward
      enum {INTERIOR, ATMOSPHERE} position= INTERIOR;
      // for each grid cell, we make 10 integration steps (this is a setting!)
      const int n= 10;
      // atmosphere pressure
      const double atm_p= rho2p(THR_RHO, S0);
      // I can start the integration from zero because I deal with it in the TOV equations
      double r= 0., p= rho2p(RHO0, S0), m= 0.;
      double dr, dp, dm, eps, rho;
      int isurf= NG+N1;
      // solve the TOV with the Euler method
      for (int i= NG; i < NG+N1; i++) {
         if (position == INTERIOR) {
            dr= (x1[i] - r) / (double) (n);
            for (int j= 0; j < n; j++) {
               // TOV equations for pressure and gravitational mass
               eps= rho2eps(p2rho(p, S0), S0);
               dp= (r <= 0.) ? 0. : -(eps + p)*(m/r + 4.*M_PI*r*r*p)/(r - 2.*m);
               dm= 4.*M_PI*eps*r*r;
               p= p + dr*dp;
               if (p <= atm_p) {
                  position= ATMOSPHERE;
                  isurf= i;
                  break;
               }
               m= m + dr*dm;
               r= r + dr;
            }
            r= x1[i];
         }
         // load the primitive variables
         if (position == INTERIOR) {
            // p is already sanitazed
            rho= p2rho(p, S0);
            var[i][0]= rho;
            var[i][2]= rho2eps(rho, S0)/rho - 1.;
         } else {
            var[i][0]= MIN_RHO;
            var[i][2]= MIN_U;
         }
         // the velocity is always 0
         var[i][1]= 0.; // Wv
         // the enclosed mass is the result of the last integration
         var[i][N0]= m;
      }
      // add a velocity perturbation to the initial condition
      if (V_PERTURBATION != 0.0) {
         const double fact= M_PI/x1[isurf];
         double v;
         for (int i= NG; i < isurf; i++) {
            v= V_PERTURBATION*sin(fact*x1[i]);
            var[i][1]= v/sqrt(1.0 - v*v);
         }
      }
   }
   // the fluid is at equilibrium
   for (int i= NG; i < NG+N1; i++) {
      if (FLUID == PERFECT && NSPECIES >= 2 && var[i][0] > MIN_RHO) {
         var[i][3]= fmax(0., fmin(1., var[i][0]/RHON*YE0));
         if (NSPECIES == 3) var[i][4]= fmax(0., fmin(1., var[i][0]/RHON*YMU0));
      } else // BULK and atmosphere
         for (int l= 3; l < N0; l++) var[i][l]= 0.;
   }
   // boundary conditions
   set_boundaries(var);
}

/// @brief check that the settings are consistent
/// @return error code (0 on success)
int check_settings(void) {
   if (N1 < 1) return error("must be N1 >= 1");
   else if (N0 < 1) return error("must be N0 > 0");
   else if (D0 < 1) return error("must be D0 > 0");
   else if (CFL <= 0. || CFL >= 1.) return error("must be 0 < CFL < 1");
   else if (END_TIME <= 0.) return error("must be END_TIME > 0");
   else if (GAM <= 1.) return error("must be GAM > 1");
   else if (GAMTH <= 1.) return error("must be GAMTH > 1");
   else if (YE0 < 0. || YE0 > 1.) return error("must be 0. <= YE0 <= 1.");
   else if (YMU0 < 0. || YMU0 > 1.) return error("must be 0. <= YMU0 <= 1.");
   else if (SLOPE_LIMITER < NONE || SLOPE_LIMITER > MCLIM) return error("unknown SLOPE_LIMITER");
   else if (ORDER_SPACE < 1 || ORDER_SPACE > 2) return error("must be 1 <= ORDER_SPACE <= 2");
   else if (PROBLEM < BJORKEN || PROBLEM > STAR) return error("unknown PROBLEM");
   else if (K0 < 0. && GAM == GAMTH && PROBLEM == FLUXTUBE)
      return error("must be K0 >= 0 for PROBLEM == FLUXTUBE and GAM == GAMTH");
   else if (K0 <= 0. && (GAM != GAMTH || PROBLEM != FLUXTUBE))
      return error("must be K0 > 0 for PROBLEM != FLUXTUBE or GAM != GAMTH");
   else if (KE < 0. && NSPECIES >= 2) return error("must be KE >= 0");
   else if (KMU < 0. && NSPECIES == 3) return error("must be KMU >= 0");
   else if (ORDER_TIME < 1 || ORDER_TIME > 3) return error("must be 1 <= ORDER_TIME <= 3");
   else if (MIN_RHO < 0.) return error("must be MIN_RHO >= 0");
   else if (MIN_RHO > THR_RHO) return error("must be MIN_RHO <= THR_RHO");
   else if (MAX_R <= 0.) return error("must be MAX_R > 0");
   else if (NSPECIES < 1 || NSPECIES > 3) return error("must be 1 <= NSPECIES <= 3");
   else if (TOL <= 0.) return error("must be TOL > 0");
   else if (TAU_BULK <= 0.) return error("must be TAU_BULK > 0.");
   else if (CHI_BULK <= 0.) return error("must be CHI_BULK > 0.");
   else if (FLUID < DUST || FLUID > BULK) return error("unknown FLUID");
   else if (R0 <= 0. || R0 >= MAX_R) return error("must be 0. < R0 < MAX_R");
   else if (RHO0 <= 0.0) return error("must be RHO0 > 0.0");
   else if (S0 < 0.0) return error("must be S0 >= 0.0");
   else if (PROBLEM == BJORKEN && FLUID != BULK) return error("PROBLEM = BJORKEN needs FLUID = BULK");
   else if (PROBLEM == FLUXTUBE && FLUID != PERFECT) return error("PROBLEM = FLUXTUBE needs FLUID = PERFECT");
   else if (FLUID == BULK && INCLUDE_LUMINOSITY && NSPECIES != 2)
      return error("INCLUDE_LUMINOSITY with FLUID = BULK needs NSPECIES = 2");
   else if (FLUID == PERFECT && NSPECIES >= 2 && PROBLEM == FLUXTUBE)
      return error("PROBLEM = FLUXTUBE not tested for FLUID = PERFECT with NSPECIES >= 2");
   else if (END_TIME <= 1. && PROBLEM == BJORKEN)
      return error("PROBLEM == BJORKEN needs END_TIME > 1");
   else return 0;
}

/// @brief print the variables on a text file in a gnuplot-friendly format
/// @param[in] time time of the snapshot
/// @param[in] x1 spatial grid
/// @param[in] var grid of primitive variables
/// @param[in] outfile output file
void output_profile(
   const double time,
   const double x1[D1],
   const double var[D1][D0],
   FILE *outfile)
{
   double ye= 0., ymu= 0.;
   if (VERBOSE) printf("VERBOSE: output_profile: writing outfile\n");
   fprintf(outfile, "# time= %f\n", time);
   for (int i= 0; i < D1; i++) {
      fprintf(outfile, "%.16e", x1[i]);
      for (int l= 0; l < D0; l++)
         fprintf(outfile, " %.16e", var[i][l]);
      if (FLUID == PERFECT) {
         if (NSPECIES >= 2) ye= var[i][3];
         if (NSPECIES == 3) ymu= var[i][4];
      }
      fprintf(outfile, " %.16e", eos2c(var[i][0],var[i][2],ye,ymu)); // cs
      fprintf(outfile, "\n");
   }
   fprintf(outfile, "\n\n");
   fflush(outfile);
}

/// @brief print central values on a text file at every timestep
/// @param[in] time time of the snapshot
/// @param[in] x1 spatial grid
/// @param[in] var primitive variables
/// @param[in] logfile output file
void output_log(
   const double time,
   const double x1[D1],
   const double var[D1][D0],
   FILE *logfile)
{
   fprintf(logfile, "%.16e", time);
   // beware: it does not print the auxiliary vars
   for (int l= 0; l < N0; l++)
      fprintf(logfile, " %.16e", var[NG][l]);
   const double ye= (FLUID == PERFECT && NSPECIES >= 2 ? var[NG][3] : 0.);
   const double ymu= (FLUID == PERFECT && NSPECIES == 3 ? var[NG][4] : 0.);
   fprintf(logfile, " %.16e", eos2c(var[NG][0],var[NG][2],ye,ymu)); // speed of sound
   if (GRAVITATING) {
      fprintf(logfile, " %.16e", var[NG+N1][N0]); // print the total gravitational mass
      fprintf(logfile, " %.16e", get_mb(x1, var)); // print the total baryonic mass
   }
   fprintf(logfile," %.16e", var[NG+N1][0]); // last density -- to check if there is mass shedding
   fprintf(logfile,"\n");
   fflush(logfile);
}

// MAIN ************************************************************************

/// @brief starting point of the program
/// @return error code (0 on success)
int main(void) {
   int err= check_settings();
   // if error, no need to finalize (I didn't open any file)
   if (err) return error("main: wrong settings");

   // settings are OK, initialize the grid and the variables
   double x1[D1], var[D1][D0];
   initialize(x1, var);
   
   // open the file with the output profile data
   FILE *outfile= fopen(OUTPUT_FILE, "w");
   // if error, no need to finalize (I didn't open any file)
   if (outfile == NULL) return error("main: outfile was not open");

   // open the logfile
   FILE *logfile= fopen(LOG_FILE, "w");
   if (logfile == NULL) {
      err= error("main: logfile was not open");
      goto finalize_2;
   }
   
   // evolution
   double time= (PROBLEM == BJORKEN ? 1. : 0.);
   double dt= 0.;
   double next_print_time= time;
   for (int i= 0; i < MAX_ITER; i++) {
      // output
      output_log(time, x1, var, logfile);
      if (time >= next_print_time) {
         output_profile(time, x1, var, outfile);
         next_print_time= next_print_time + OUTPUT_INTERVAL;
      }

      // determine the timestep
      dt= timestep(x1, var);
      if (dt <= 0.) {
         err= error("main: non-positive timestep: dt <= 0");
         goto finalize;
      }
      dt= fmin(dt, fmin(END_TIME, next_print_time) - time);

      if (VERBOSE)
         printf("VERBOSE: main: iter= %d, time= %f, timestep= %f\n", i, time, dt);

      // evolve
      err= imex_integration(time, dt, x1, var);
      if (err) goto finalize;
      time+= dt;
      if (time >= END_TIME) break;
   }
   if (time < END_TIME)
      printf("WARNING: main: max iterations reached, simulation didn't end at END_TIME\n");
   
   output_log(time, x1, var, logfile);
   output_profile(time, x1, var, outfile);

   // unwind and return
   finalize:
   fclose(logfile);
   finalize_2:
   fclose(outfile);
   return err;
}
