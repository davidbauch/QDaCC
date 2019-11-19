#pragma once

#include "global.h"

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

    std::vector<SaveState> savedStates;    // Vector for saved matrix-time tuples for densitymatrix
    std::vector<SaveState> savedHamiltons; // Vector for saved matrix-time tuples for hamilton operators
    std::vector<dcomplex> out;
    void saveState( const MatType &mat, const double t, std::vector<SaveState> &savedStates );
    void saveHamilton( const MatType &mat, const double t );
    MatType iterateRungeKutta4( const MatType &rho, System_Parent &s, const double t, std::vector<SaveState> &savedStates );
    MatType iterateRungeKutta5( const MatType &rho, System_Parent &s, const double t, std::vector<SaveState> &savedStates );
    int getIterationNumberTau( System_Parent &s );
    int getIterationNumberSpectrum( System_Parent &s );
    double getTimeAt( int i );
    MatType getRhoAt( int i );
    MatType getHamilton( System_Parent &s, const double t, bool use_saved_hamiltons );
    bool queueNow( System_Parent &s, int &curIt );
    int reset( System_Parent &s );
    bool calculate_g1( System_Parent &s, const MatType &op_creator, const MatType &op_annihilator, DenseMat &cache, std::string purpose = "unknown" );
    bool calculate_g2( System_Parent &s, const MatType &op_creator_1, const MatType &op_annihilator_1, const MatType &op_creator_2, const MatType &op_annihilator_2, DenseMat &cache, std::string purpose = "unknown" );

        public : ODESolver(){};
    ODESolver( System_Parent &s );
    MatType iterate( const MatType &rho, System_Parent &s, const double t, std::vector<SaveState> &savedStates, const int dir );
    bool calculate_t_direction( System_Parent &s );
    //bool calculate_g2_0( System_Parent &s, const MatType &op_creator, const MatType &op_annihilator, std::string fileOutputName ); //moved to advancedPhotonStatistics
    bool calculate_spectrum( System_Parent &s, const MatType &op_creator, const MatType &op_annihilator, std::string fileOutputName );
    bool calculate_advanced_photon_statistics( System_Parent &s, const MatType &op_creator_1, const MatType &op_annihilator_1, const MatType &op_creator_2, const MatType &op_annihilator_2, std::string fileOutputName );
    template <typename T>
    static MatType iterate_definite_integral( const MatType &rho, T rungefunction, const double t, const double step );
    static std::vector<SaveState> calculate_definite_integral_vec( MatType rho, std::function<MatType( const MatType &, const double )> const &rungefunction, const double t0, const double t1, const double step );
    static SaveState calculate_definite_integral( MatType rho, std::function<MatType( const MatType &, const double )> const &rungefunction, const double t0, const double t1, const double step );
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