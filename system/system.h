#pragma once
// Dependencies
#include "operatormatrices.h"
#include "parameters.h"

class System_Parent {
   public:
    // Variables
    std::string name;
    std::string terminate_message;
    std::vector<std::string> arguments;
    Parameters parameters;
    OperatorMatrices operatorMatrices;
    // Constructor
    System_Parent(){};
    System_Parent( const std::vector<std::string> &input ){};
    // Functions:
    virtual double getTimeborderStart() { return parameters.t_start; };
    virtual double getTimeborderEnd() { return parameters.t_end; };
    virtual double getTimeStep() { return parameters.t_step; };
    virtual MatrixXcd getRho0() { return operatorMatrices.rho; }; // Rho is always complex
    // Matrix Commutator function
    template <typename T>
    T dgl_kommutator( const T &A, const T &B ) {
        return A * B - B * A;
    }
    // Matrix Anticommutator function
    template <typename T>
    T dgl_antikommutator( const T &A, const T &B ) {
        return A * B + B * A;
    }
    // Matrix Lindblad Decay Term
    template <typename T, typename T2, typename T3>
    T dgl_lindblad( const T &rho, const T2 &op, const T3 &opd ) {
        return 2.0 * op * rho * opd - opd * op * rho - rho * opd * op;
    }
    //virtual MatrixXcd dgl_timetrafo( const DenseMat &op, const double t ) { return op; };
    virtual DenseMat dgl_timetrafo( const DenseMat &op, const double t ) { return op; };
    // Expectationvalue
    template <typename T>
    std::complex<double> dgl_expectationvalue( const T &rho, const T &op, const double t ) {
        return ( rho * dgl_timetrafo( op, t ) ).trace();
    }
    // Time Transformation ( -> Interaction picture)
    // Transformed densitymatrix (q.r.t)
    template <typename T>
    T dgl_calc_rhotau( const T &rho, const T &op, const double t ) {
        return dgl_timetrafo( op, t ) * rho;
    }
    // Runge function
    virtual DenseMat dgl_rungeFunction( const DenseMat &rho, const MatrixXcd &H, const double t, std::vector<SaveState> &past_rhos ) {
        logs( "Warning; Returning virtual function!\n" );
        return rho;
    }
    // Output Expectation Values
    virtual void expectationValues( const DenseMat &rho, const double t ){};
    // Check for trace validity
    bool traceValid( DenseMat &rho, const double t, const bool forceOutput );
    // Return Hamiltonian for time t
    virtual DenseMat dgl_getHamilton( const double t ) {
        logs( "Warning; Returning virtual function\n" );
        return DenseMat( 1, 1 );
    };
    // Return solver order for time direction
    int getSolverOrder( const int dir );
    // To be called from child constructor
    void init();
    // Initial init. Should return true if initialization seemed valid.
    virtual bool init_system() { return false; };
    // Final exit. Should close any open files, should return true if seemed valid.
    virtual bool exit_system( const int failure ) { return false; };
    // Output predetermined parameter variables:
    bool calculate_spectrum();
    bool calculate_g2();
    bool use_interactionpicture();
    bool use_rwa();
    bool output_handlerstrings();
    bool output_operators();
    int getSolverRungeKuttaOrder( int dir );
    int getTimeTransformationAlg();
    int getIterationSkip( );
};