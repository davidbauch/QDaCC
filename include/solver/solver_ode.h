#pragma once

#include "global.h"
#include "system/system.h"
#include "misc/interpolant.h"

// Description: ODESolver class provides both Runge-Kutta functions of different orders and functions for different numerical operations
class ODESolver {
   private:
    // RK 4&5 coefficients
    double a2 = 1. / 5.;
    double a3 = 3. / 10.;
    double a4 = 4. / 5.;
    double a5 = 8. / 9.;
    double a6 = 1.0;

    double b11 = 1. / 5.;
    double b21 = 3. / 40.;
    double b31 = 44. / 45.;
    double b41 = 19372. / 6561.;
    double b51 = 9017. / 3168.;
    double b61 = 35. / 384.;
    double b22 = 9. / 40.;
    double b32 = -56. / 15.;
    double b42 = -25360. / 2187.;
    double b52 = -355. / 33.;
    double b33 = 32. / 9.;
    double b43 = 64448. / 6561.;
    double b53 = 46732. / 5247.;
    double b63 = 500. / 1113.;
    double b44 = -212. / 729.;
    double b54 = 49. / 176.;
    double b64 = 125. / 192.;
    double b55 = -5103 / 18656.;
    double b65 = -2187. / 6784.;
    double b66 = 11. / 84.;

    int track_gethamilton_read, track_gethamilton_write, track_gethamilton_calc, track_gethamilton_calcattempt;
    int dim;
    // Cache Matrices that allow for variable time steps.
    Dense cache1 = Dense::Zero( 1, 1 );
    Dense cache2 = Dense::Zero( 1, 1 );

    // Cached Entries
    std::vector<SaveState> savedStates;    // Vector for saved matrix-time tuples for densitymatrix
    std::vector<SaveState> savedHamiltons; // Vector for saved matrix-time tuples for hamilton operators

    // Description: Saves a tuple of a complex (density-)matrix and time, ensuring times and matrices don't get mixed up
    // Type: ODESolver private function
    // @param mat: [&Sparse] Matrix to save
    // @param t: [double] Corresponding time
    // @return: [void]
    void saveState( const Sparse &mat, const double t, std::vector<SaveState> &savedStates );

    // Description: Saves a tuple of a complex (Hamilton-)matrix and time, ensuring times and matrices don't get mixed up
    // Type: ODESolver private function
    // @param mat: [&Sparse] Matrix to save
    // @param t: [double] Corresponding time
    // @return [void]
    void saveHamilton( const Sparse &mat, const double t );

    // Description: Iterates Runge-Kutta of order 4 at time t onto rho using the systems hamilton operator.
    // Type: ODESolver private function
    // @param rho: [&Sparse] Input (density-) matrix
    // @param s: [&System] Class providing set of system functions
    // @param t: [double] Time to iterate at
    // @return: [Sparse] rho at time t+t_step
    Sparse iterateRungeKutta4( const Sparse &rho, System &s, const double t, std::vector<SaveState> &savedStates );
    // Description: Iterates Runge-Kutta of order 5 at time t onto rho using the systems hamilton operator.
    // Type: ODESolver private function
    // @param rho: [&Sparse] Input (density-) matrix
    // @param s: [&System] Class providing set of system functions
    // @param t: [double] Time to iterate at
    // @return: [Sparse] rho at time t+t_step
    Sparse iterateRungeKutta5( const Sparse &rho, System &s, const double t, std::vector<SaveState> &savedStates );

    // Desciption: Function to calculate the number of iterations used for tau direction calculations
    // Type: ODESolver private function
    // @param s: [&System] Class providing set of system functions
    // @return [int] Number of tau-direction iterations
    int getIterationNumberTau( System &s );

    int getIterationNumberSpectrum( System &s );
    double getTimeAt( int i );
    Sparse getRhoAt( int i );

    // Description: Gatheres a Hamiltonoperator using a systems get_hamilton() function. This workaround enables the saving of already calculated matrices for dublicate uses. Seems to have almost no influence on runtime.
    // Type: ODESolver private function
    // @param s: [&System] Class providing set of system functions
    // @param t: [double] Time to return Hamiltonoperator at
    // @return: [Sparse] hamilton matrix of type Sparse
    Sparse getHamilton( System &s, const double t, bool use_saved_hamiltons = false );

    // Description: Checks wether or not to save the current matrix for later calculations
    // Type: ODESolver private function
    // @param s: [&System] Class providing set of system functions
    // @param curIt: [&int] Iterator variable to use for iteration check
    // @return [bool] True if current matrix should be saved, else false. Saving depends on the iteration skip for tau direction and wether or not they are needed later.
    bool queueNow( System &s, int &curIt );

    // Description: Resets all temporary variables. Currently: out, akf_mat, track variables
    // Type: ODESolver private function
    // @param s: [&System] Class providing set of system functions
    // @return: [int] dimensions of temporary variables
    int reset( System &s );

    bool calculate_g1( System &s, const Sparse &op_creator, const Sparse &op_annihilator, Dense &cache, std::string purpose = "unknown" ); //std::vector<std::vector<SaveScalar>>
    bool calculate_g2( System &s, const Sparse &op_creator_1, const Sparse &op_annihilator_1, const Sparse &op_creator_2, const Sparse &op_annihilator_2, Dense &cache, std::string purpose = "unknown" );

    // Description: Stretches Data evaluated at various timesteps onto an equidistant grid, such that g1 and g2 can work with said grid.
    bool scale_grid( System &s, Dense &cache, std::vector<std::vector<SaveScalar>> &cache_noneq );

   public:
    ODESolver(){};
    ODESolver( System &s );

    // Description: Iterates Runge-Kutta with given order depending on the systems settings.
    // Type: ODESolver public function
    // @param rho: [&Sparse] Input (density-) matrix
    // @param s: [&System] Class providing set of system functions
    // @param t: [double] Time to iterate at
    // @param dir: [int] Time direction. Solver order can differ in both t and tau direction, depending on the systems settings. Default value is DIR_T
    // @return: [Sparse] rho at time t+t_step
    Sparse iterate( const Sparse &rho, System &s, const double t, std::vector<SaveState> &savedStates, const int dir = DIR_T );

    // Description: Calculates the normal t-direction via solving the von-Neumann equation for rho. May save some of the density matrices for later uses. Logs the calculation and outputs progress.
    // Type: ODESolver public function
    // @param s: [&System] Class providing set of system functions
    // @return: [bool] True if calculations are sucessfull, else false
    bool calculate_t_direction( System &s );

    //bool calculate_g2_0( System &s, const Sparse &op_creator, const Sparse &op_annihilator, std::string fileOutputName ); //moved to advancedPhotonStatistics
    bool calculate_spectrum( System &s, const Sparse &op_creator, const Sparse &op_annihilator, std::string fileOutputName, const int cache_index );
    bool calculate_advanced_photon_statistics( System &s, const Sparse &op_creator_1, const Sparse &op_annihilator_1, const Sparse &op_creator_2, const Sparse &op_annihilator_2, std::string fileOutputName );
    //template <typename T>
    //static Sparse iterate_definite_integral( const Sparse &rho, T rungefunction, const double t, const double step );
    //static std::vector<SaveState> calculate_definite_integral_vec( Sparse rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step );
    //static SaveState calculate_definite_integral( Sparse rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step );

    bool calculate_runge_kutta( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output = true );
};

/*
maddlab
double a2 = 1./5.;
double a3 = 3./10.;
double a4 = 4./5.;
double a5 = 8./9.;
double a6 = 1.0;

double b11 = 1./5.;
double b21 = 3./40.;
double b31 = 44./45.;
double b41 = 19372./6561.;
double b51 = 9017./3168.;
double b61 = 35./384.;
double b22 = 9./40.;
double b32 = -56./15.;
double b42 = -25360./2187.;
double b52 = -355./33.;
double b33 = 32./9.;
double b43 = 64448./6561.;
double b53 = 46732./5247.;
double b63 = 500./1113.;
double b44 = -212./729.;
double b54 = 49./176.;
double b64 = 125./192.;
double b55 = -5103/18656.;
double b65 = -2187./6784.;
double b66 = 11./84.;

a2=cast(1/5,dataType);
a3=cast(3/10,dataType);
a4=cast(4/5,dataType);
a5=cast(8/9,dataType);

b11=cast(1/5,dataType); 
b21=cast(3/40,dataType); 
b31=cast(44/45,dataType);
b41=cast(19372/6561,dataType);
b51=cast(9017/3168,dataType);
b61=cast(35/384,dataType);
b22=cast(9/40,dataType);
b32=cast(-56/15,dataType);
b42=cast(-25360/2187,dataType);
b52=cast(-355/33,dataType);
b33=cast(32/9,dataType);
b43=cast(64448/6561,dataType);
b53=cast(46732/5247,dataType);
b63=cast(500/1113,dataType);
b44=cast(-212/729,dataType);
b54=cast(49/176,dataType);
b64=cast(125/192,dataType);
b55=cast(-5103/18656,dataType);
b65=cast(-2187/6784,dataType);
b66=cast(11/84,dataType);
*/