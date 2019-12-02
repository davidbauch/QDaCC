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
    virtual MatType getRho0() { return operatorMatrices.rho; };
    template <typename T>
    inline T getTrace( const DenseMat &mat ) const {
        return mat.trace();
    }
    template <typename T>
    inline T getTrace( const SparseMat &mat ) const {
        return getTrace<T>( DenseMat( mat ) );
    }
    // Matrix Commutator function
    template <typename T>
    inline T dgl_kommutator( const T &A, const T &B ) {
        return A * B - B * A;
    }
    // Matrix Anticommutator function
    template <typename T>
    inline T dgl_antikommutator( const T &A, const T &B ) {
        return A * B + B * A;
    }
    // Matrix Lindblad Decay Term
    template <typename T, typename T2, typename T3>
    inline T dgl_lindblad( const T &rho, const T2 &op, const T3 &opd ) {
        return 2.0 * op * rho * opd - opd * op * rho - rho * opd * op;
    }
    //virtual MatType dgl_timetrafo( const MatType &op, const double t ) { return op; };
    virtual MatType dgl_timetrafo( const MatType &op, const double t ) { return op; };
    // Expectationvalue
    template <typename T, typename T2>
    inline T2 dgl_expectationvalue( const T &rho, const T &op, const double t ) {
        return getTrace<T2>( ( rho * dgl_timetrafo( op, t ) ).eval() );
    }
    // Time Transformation ( -> Interaction picture)
    // Transformed densitymatrix (q.r.t)
    template <typename T>
    inline T dgl_calc_rhotau( const T &rho, const T &op, const double t ) {
        return dgl_timetrafo( op, t ) * rho;
    }
    // Order 2
    template <typename T>
    inline T dgl_calc_rhotau_2( const T &rho, const T &op1, const T &op2, const double t ) {
        return dgl_timetrafo( op1, t ) * rho * dgl_timetrafo( op2, t );
    }
    // Runge function
    virtual MatType dgl_rungeFunction( const MatType &rho, const MatType &H, const double t, std::vector<SaveState> &past_rhos ) {
        return rho;
    };
    // Output Expectation Values
    virtual void expectationValues( const MatType &rho, const double t, const std::vector<SaveState> &pastrhos ){};
    // Check for trace validity
    template <class T>
    bool traceValid( T &rho, double t_hit, bool force = false ) {
        double trace = std::real( getTrace<dcomplex>( rho ) );
        parameters.trace.emplace_back( trace );
        if ( trace < 0.99 || trace > 1.01 || force ) {
            if ( force )
                fmt::print( "{} {} -> trace check failed at t = {} with trace(rho) = {}\n", PREFIX_ERROR, global_message_error_divergent, t_hit, trace );
            terminate_message = global_message_error_divergent;
            parameters.numerics_calculate_spectrum = 0;
            FILE *fp_trace = std::fopen( ( parameters.subfolder + "trace.txt" ).c_str(), "w" );
            for ( int i = 0; i < (int)parameters.trace.size() && parameters.t_step * 1.0 * i < t_hit; i++ ) {
                fmt::print( fp_trace, "{:.10e} {:.15e}\n", parameters.t_step * 1.0 * ( i + 1 ), parameters.trace.at( i ) );
            }
            std::fclose( fp_trace );
            return false;
        } else {
            return true;
        }
    }
    // Return Hamiltonian for time t
    virtual MatType dgl_getHamilton( const double t ) {
        logs( "Warning; Returning virtual function\n" );
        return MatType( 1, 1 );
    };
    // Return solver order for time direction
    int getSolverOrder( const int dir );
    // To be called from child constructor
    void init();
    // Initial init. Should return true if initialization seemed valid.
    virtual bool init_system() { return false; };
    // Final exit. Should close any open files, should return true if seemed valid.
    virtual bool exit_system( const int failure ) { return false; };
    // Return specific operators:
    virtual MatType getOperator( std::string &argument ) { return MatType( 1, 1 ); }
    // Execute specific commands:
    virtual bool command( unsigned int index ) { return true; }
    // Output predetermined parameter variables:
    virtual bool calculate_spectrum() {
        return parameters.numerics_calculate_spectrum;
    }
    bool calculate_g2();
    bool use_interactionpicture();
    bool use_rwa();
    bool output_handlerstrings();
    bool output_operators();
    int getSolverRungeKuttaOrder( int dir );
    int getTimeTransformationAlg();
    int getIterationSkip();
};