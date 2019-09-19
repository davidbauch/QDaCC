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

    int track_gethamilton_read, track_gethamilton_write, track_gethamilton_calc,track_gethamilton_calcattempt;

    // Vector for mat/time, tuple
    class SaveState {
       public:
        MatrixXcd mat;
        double t;
        SaveState( const MatrixXcd &mat, const double time ) : mat(mat), t(time) {};
    };

    std::vector<SaveState> savedStates; // Vector for saved matrix-time tuples for densitymatrix
    std::vector<SaveState> savedHamiltons; // Vector for saved matrix-time tuples for hamilton operators
    std::vector<std::complex<double>> out;
    MatrixXcd akf_mat;
    void saveState( const MatrixXcd &mat, const double t );
    void saveHamilton( const MatrixXcd &mat, const double t );
    bool queueNow( const System &s, int &curIt );
    MatrixXcd iterateRungeKutta4( const MatrixXcd &rho, const System &s, const double t );
    MatrixXcd iterateRungeKutta5( const MatrixXcd &rho, const System &s, const double t );
    int getIterationNumberTau( const System &s );
    int getIterationNumberSpectrum( const System &s );
    double getTimeAt( int i );
    MatrixXcd getRhoAt( int i );
    MatrixXcd getHamilton( const System &s, const double t, bool use_saved_hamiltons );
    int reset( const System &s );

   public:
    ODESolver(){};
    ODESolver( const System &s );
    MatrixXcd iterate( const MatrixXcd &rho, const System &s, const double t, const int dir );
    bool calculate_t_direction( System &s );
    bool calculate_g1( const System &s, const MatrixXcd &op_creator, const MatrixXcd &op_annihilator );
    bool calculate_g2_0( const System &s, const MatrixXcd &op_creator, const MatrixXcd &op_annihilator, std::string fileOutputName );
    bool calculate_spectrum( const System &s, std::string fileOutputName );
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