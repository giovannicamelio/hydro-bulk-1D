/// @file
/// @brief parameter file for the migration of Maxwell-Cattaneo with luminosity
/// @details the settings in this files determine the code behaviour
/// @author Giovanni Camelio
/// @date 2022

/// @brief slope limiter
/// @details used if @ref ORDER_SPACE == 2
static const enum {
   NONE,   ///< no slope limiter
   MINMOD, ///< use the `minmod' slope limiter
   MCLIM   ///< use the `MC' slope limiter
} SLOPE_LIMITER= MINMOD;

/// @brief type of fluid
static const enum {
   DUST,    ///< pressurless fluid
   PERFECT, ///< perfect fluid
   BULK     ///< bulk stress
} FLUID= BULK;

/// @brief problem to solve (including initial conditions)
static const enum {
   BJORKEN,   ///< Bjorken flow: Smirne coordinates and no fluxes
   FLUXTUBE,  ///< hydrodynamics in cartesian coordinates
   CLOUD,     ///< hydrodynamics with initial homogeneous pressurless cloud
   STAR       ///< hydrodynamics with initial solution of TOV eqs
} PROBLEM= STAR;

/// @brief verbosity of the output on stdout
/// @details true (verbose) / false (quiet)
#define VERBOSE true

/// @brief type of Riemann solver
/// @details true (LLF) / false (HLL),
///          used if @ref PROBLEM != @ref BJORKEN
#define USE_LLF_RIEMANN_SOLVER true

/// @brief control flag for the treatment of the bulk
/// @details true (Hiscock-Lindblom) / false (Maxwell-Cattaneo),
///          used if @ref FLUID == @ref BULK
#define HISCOCK_LINDBLOM false

/// @brief threshold for the temperature in the differentiation
/// @details double, @ref TINY_T > 0.0,
///          used if @ref HISCOCK_LINDBLOM == true
#define TINY_T 1e-61

/// @brief control flag for the inclusion of the energy luminosity
/// @details true (include luminosity) / false (no luminosity)
#define INCLUDE_LUMINOSITY true

/// @brief control flag for the type of reaction considered
/// @details true (dURCA) / false (mURCA),
///          used if @ref NSPECIES >= 2
#define USE_DURCA_REACTIONS true

/// @brief order of the spatial reconstruction
/// @details integer, 1 <= ORDER_SPACE <= 2,
///          useful if @ref PROBLEM != @ref BJORKEN
#define ORDER_SPACE 2

/// @brief order of the time integration
/// @details integer, 1 <= ORDER_TIME <= 3
#define ORDER_TIME 3

/// @brief number of independent species
/// @details integer, 1 <= NSPECIES <= 3
#define NSPECIES 2

/// @brief number of physical grid points (without ghosts)
/// @details integer, N1 > 0
#define N1 8001

/// @brief file name of the output profiles
/// @details is a string
#define OUTPUT_FILE "mig-mcq.dat"

/// @brief file name of output log
/// @details is a string
#define LOG_FILE "mig-mcq.log"

/// @brief Courant-Friedrichs-Lewy factor [#]
/// @details double, 0.0 < CFL < 1.0
#define CFL 0.5

/// @brief final time of the simulation [M☉]
/// @details double, END_TIME > 0.0 (for @ref BJORKEN, END_TIME > 1.0)
#define END_TIME (0.005*SEC)

/// @brief exponent of the gamma-law EOS [#]
/// @details double, GAM > 1.0,
///          used if @ref FLUID != @ref DUST
#define GAM (2.)

/// @brief thermal exponent of the gamma-law EOS [#]
/// @details double, GAMTH > 1.0
///          used if @ref FLUID != @ref DUST
#define GAMTH (1.52)

/// @brief minimal value of rest mass density [M☉⁻²]
/// @details double, MIN_RHO >= 0.0
#define MIN_RHO (1e-20)

/// @brief threeshold density of the atmosphere [M☉⁻²]
/// @details double, THR_RHO >= @ref MIN_RHO >= 0.0
#define THR_RHO (100.*MIN_RHO)

/// @brief maximal number of time iterations
/// @details integer, @ref MAX_ITER > 0
#define MAX_ITER 1000000

/// @brief output frequency of the full profile
/// @details is a double
#define OUTPUT_INTERVAL (0.05*END_TIME)

/// @brief EOS polytropic constant [c=G=M☉=kB=1]
/// @details double, K0 > 0.0; to increase the precision it can be K0 = 0.0
///          for @ref PROBLEM == @ref FLUXTUBE and @ref GAM == @ref GAMTH
#define K0 (200.)

/// @brief EOS polytropic constant for electrons [c=G=M☉=kB=1]
/// @details double, KE > 0.0
#define KE (0.603)

/// @brief EOS polytropic constant for muons [c=G=M☉=kB=1]
/// @details double, KMU > 0.0
///          @li for modified muons KMU = 0.147
///          @li for fitted muons KMU = 0.974
#define KMU (0.147)

/// @brief EOS polytropic constant for the temperature [c=G=M☉=kB=1]
/// @details double, KTH > 0.0
#define KTH (0.374)

/// @brief electron fraction at equilibrium [#]
/// @details double, 0.0 <= YE0 <= 1.0
#define YE0 (0.0570)

/// @brief muon fraction at equilibrium [#]
/// @details double, 0.0 <= YMU0 <= 1.0
///          @li for modified muons YMU0 = YE0 = 0.0570
///          @li for fitted muons YMU0 = 0.00860
#define YMU0 (YE0)

/// @brief central rest mass density [M☉⁻²]
/// @details double, RHO0 > 0.0
///          used if @ref PROBLEM != @ref SHOCKTUBE
#define RHO0 (4e-3)

/// @brief maximal radius of the physical (without ghost) grid [M☉]
/// @details double, MAX_R > 0.0,
///          used if @ref PROBLEM >= @ref CLOUD
#define MAX_R (80.)

/// @brief relative tolerance in the cons2prim inversion [#]
/// @details double, TOL > 0.0
#define TOL 1e-15

/// @brief bulk viscous parameter [M☉⁻²]
/// @details double, CHI_BULK > 0.0,
///          used if @ref FLUID == @ref BULK and @ref NSPECIES == 1
#define CHI_BULK (100./RHON)

/// @brief bulk viscous timescale [M☉]
/// @details double, CHI_TAU > 0.0,
///          used if @ref FLUID == @ref BULK and @ref NSPECIES == 1
#define TAU_BULK (1.)

/// @brief initial bulk viscous stress [M☉⁻²]
/// @details double, used if @ref PROBLEM == @ref BJORKEN and @ref NSPECIES == 1
#define PI0 (1e-8)

/// @brief left density in the shock tube
/// @details double, used if @ref PROBLEM == @ref FLUXTUBE
#define RHOL 0.9

/// @brief right density in the shock tube
/// @details double, used if @ref PROBLEM == @ref FLUXTUBE
#define RHOR 0.1

/// @brief left pressure in the shock tube
/// @details double, used if @ref PROBLEM == @ref FLUXTUBE
#define PL 1.0

/// @brief right pressure in the shock tube
/// @details double, used if @ref PROBLEM == @ref FLUXTUBE
#define PR 0.001

/// @brief left velocity in the shock tube
/// @details double, used if @ref PROBLEM == @ref FLUXTUBE
#define VL 0.

/// @brief right velocity in the shock tube
/// @details double, used if @ref PROBLEM == @ref FLUXTUBE
#define VR 0.

/// @brief initial radius for the dust collapse
/// @details double, R0 < MAX_R, used if @ref PROBLEM == @ref CLOUD,
#define R0 9.

/// @brief initial entropy per baryon
/// @details double, S0 >= 0.0, used if @ref PROBLEM == @ref BJORKEN
#define S0 0.

/// @brief amplitude of the initial velocity perturbation
/// @detail @li double,
///         @li used if @ref PROBLEM = @ref STAR,
///         @li v(r)= V_PERTURBATION * sin(π * r / R),
///         @li if V_PERTURBATION == 0.0, there is no initial perturbation
#define V_PERTURBATION (0.0)
