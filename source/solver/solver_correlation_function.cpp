#include "solver/solver_ode.h"

double QDACC::Numerics::get_tdelta( const Dense &gmat_time, size_t fixed_index, size_t var_index ) {
    return var_index == 0 ? std::real( gmat_time( var_index + 1, fixed_index ) - gmat_time( var_index, fixed_index ) ) : std::real( gmat_time( var_index, fixed_index ) - gmat_time( var_index - 1, fixed_index ) );
}
double QDACC::Numerics::get_taudelta( const Dense &gmat_time, size_t fixed_index, size_t var_index ) {
    return var_index == 0 ? std::imag( gmat_time( fixed_index, var_index + 1 ) - gmat_time( fixed_index, var_index ) ) : std::imag( gmat_time( fixed_index, var_index ) - gmat_time( fixed_index, var_index - 1 ) );
}
double QDACC::Numerics::get_tdelta( const std::vector<SaveState> &savedStates, size_t var_index ) {
    return var_index == 0 ? savedStates[var_index + 1].t - savedStates[var_index].t : savedStates[var_index].t - savedStates[var_index - 1].t;
}

void QDACC::Numerics::ODESolver::calculate_g1( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[CorrelationFunction] G1(tau) for {} already exists.\n", purpose );
    }
    calculate_g1( s, std::vector<std::string>{ s_op_i }, s_op_j, { purpose } );
}

void QDACC::Numerics::ODESolver::calculate_g1( System &s, const std::vector<std::string> &s_op_i, const std::string &s_op_j, const std::vector<std::string> &purposes ) {
    const auto size = s_op_i.size();
    calculate_g2( s, "internal_identitymatrix", std::vector<std::string>{ size, "internal_identitymatrix" }, s_op_i, s_op_j, purposes );
}

void QDACC::Numerics::ODESolver::calculate_g2( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[CorrelationFunction] G2(tau) for {} already exists.\n", purpose );
    }
    calculate_g2( s, s_op_i, std::vector<std::string>{ s_op_j }, { s_op_k }, s_op_l, { purpose } );
}

void QDACC::Numerics::ODESolver::calculate_g2( System &s, const std::string &s_op_i, const std::vector<std::string> &s_op_j, const std::vector<std::string> &s_op_k, const std::string &s_op_l, const std::vector<std::string> &purposes ) {
    Log::L2( "[CorrelationFunction] Preparing to calculate Correlation function\n" );
    Log::L2( "[CorrelationFunction] Generating Sparse Operator Matrices from String input...\n" );
   
    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_l = get_operators_matrix( s, s_op_l );
    Log::L2( "[CorrelationFunction] Iterated Operators are op_i = {} and op_l = {}\n", s_op_i, s_op_l );
    
    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Generator super-purpose
    std::string super_purpose = std::accumulate( std::next( purposes.begin() ), purposes.end(), purposes.front(), []( const std::string &a, const std::string &b ) { return a + " and " + b; } );

    std::vector<std::pair<Sparse, std::string>> eval_operators;
    
    // TEST:
    MultidimensionalCacheMatrix cache_test( {(int)matdim,(int)matdim}, "test" );

    // Preconstruct
    for ( auto current = 0; current < s_op_k.size(); current++ ) {
        const auto &purpose = purposes[current];
        // Cancel if Purpose already exists
        if ( cache.contains( purpose ) ) {
            Log::L2( "[CorrelationFunction] Matrix for {} already exists! Skipping!\n", purpose );
            continue;
        }
        const auto &op_j = get_operators_matrix( s, s_op_j[current] );
        const auto &op_k = get_operators_matrix( s, s_op_k[current] );
        // Construct Evaluation Operators
        eval_operators.emplace_back( op_j * op_k, purpose );

        Log::L2( "[CorrelationFunction] Preparing Cache Matrices for {}. Using Eval Operators op_j = {}, op_k = {}\n", purpose, s_op_j[current], s_op_k[current] );
        cache[purpose] = MultidimensionalCacheMatrix( {(int)matdim, (int)matdim}, purpose );
        auto &mat = cache[purpose];
        // Fill Time Matrix
#pragma omp parallel for collapse( 2 ) schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
        for ( size_t i = 0; i < matdim; i++ ) {
            for ( size_t j = 0; j < matdim; j++ ) {
                double tau = i + j < matdim ? s.parameters.grid_values[i + j] : s.parameters.grid_values.back() + s.parameters.grid_steps.back() * ( i + j - matdim + 1 );
                mat.get_time( i, j ) = s.parameters.grid_values[i] + 1.0i * tau;
            }
        }
    }

    // Return if no Functions need to be evaluated
    if ( eval_operators.empty() )
        return;

    // Create Timer and Progresbar
    Timer &timer = Timers::create( "Correlation-Loop (" + super_purpose + ")" ).start();
    auto progressbar = ProgressBar();

    // Calculate G2 Function
    Log::L2( "[CorrelationFunction] Calculating G(tau)... purpose: {}, saving to matrix of size {}x{},  iterating over {} saved states...\n", super_purpose, matdim, matdim, std::min<size_t>( matdim, savedStates.size() ) );
    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t i = 0; i < std::min<size_t>( matdim, savedStates.size() ); i++ ) {
        // Create and reserve past rho's vector
        std::vector<QDACC::SaveState> savedRhos;
        // Get Time from saved State
        double t_t = get_time_at( i );
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau_2( get_rho_at( i ), op_l, op_i, t_t );
        // Calculate Runge Kutta or PI
        if ( s.parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
            calculate_path_integral_correlation( rho_tau, t_t, s.parameters.t_end, timer, s, savedRhos, op_l, op_i, 1 /*ADM Cores*/ );
        } else {
            calculate_runge_kutta( rho_tau, t_t, s.parameters.t_end, timer, progressbar, super_purpose, s, savedRhos, false /*Output*/ );
        }
        // Interpolate saved states to equidistant timestep
        savedRhos = Numerics::interpolate_curve( savedRhos, t_t, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        for ( const auto &[eval, purpose] : eval_operators ) {
            auto &gmat = cache[purpose];
            for ( size_t j = 0; j < savedRhos.size(); j++ ) {
                const double t_tau = savedRhos.at( j ).t;
                gmat.get( i, j ) = s.dgl_expectationvalue<Sparse, Scalar>( savedRhos.at( j ).mat, eval, t_tau );
            } 
        }
        Timers::outputProgress( timer, progressbar, i, savedStates.size(), super_purpose );
    }

 
    timer.end();
    Timers::outputProgress( timer, progressbar, savedStates.size(), savedStates.size(), super_purpose, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[CorrelationFunction] G ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", super_purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );

    // Manually Apply the detector function
    for ( const auto &[eval, purpose] : eval_operators ) {
        auto &gmat = cache[purpose];
        if (gmat.hasBeenFourierTransformed()){ 
            Log::L2( "[CorrelationFunction] Detector Function already applied to {}. Skipping!\n", purpose);
            continue;
        }
        apply_detector_function( s, gmat, purpose );
    }
}

// Future TODO: am besten sollte es nur eine funktion geben, die (a rho b)t->t+tau ausrechnet, und das speichert. dann können calculate_g_i functionen nur noch die summation übernehmen. but idk.
// Dafür: cachematrix_matrix klasser, die vector von rho matrices speichert.

void QDACC::Numerics::ODESolver::calculate_g3( System &s, const std::string &s_op_i, const std::vector<std::string> &s_op_j, const std::vector<std::string> &s_op_k, const std::vector<std::string> &s_op_l, const std::vector<std::string> &s_op_m, const std::string &s_op_n, const std::vector<std::string> &purposes ) {

/*
    für G3 auch relevant:
    time shifted operators.
    dafür soll es möglich sein, ein t_min und t_max der korrelationsfunktion zu definieren
    und dann soll die korrelationsfunktion nur für diese zeiten berechnet werden.
    ein operator kann dann einen timeshift haben, beispielsweise wenn die dynamiken ungefär periodisch sein sollen

*/

/*
    vorgehen:
    in verschachteltem loops: 
    für jedes rho in rho(t):
        a_n rho(t) a_i^d nehmen und zu R(tau_1) integrieren
        dabei entstehen tau_1/dt matritzen, die in einem vector gespeichert werden

        for jedes r in R(tau_1):
            r*a_j^d*a_l nehmen und zu Q(tau_2) integrieren
            dabei entstehen tau_2/dt matritzen, die wieder in einem vector gespeichert werden
            
            für jedes q in Q(tau_1):
                erwartungswert von q*a_k^d*a_l ausrechnen und in g(tau_1,tau_2) speichern
                hierfür benötigen wir eine cacheMatrix klasse, die t,tau_1,tau_2 triplets speichern kann.
    
    am besten speichern wir die werte nicht mehr in eigen matritzen, sondern in einem (indexvector->index),value array.
    da G^n dense ist, sind die indices einfach ordered nach t,tau,tau2,...
    so können wir später auch noch G^N implementieren, wenn wir das wollen, indem wir den obigen loop einfach rekusiv aufrufen.

    TODO: cacheMatrixVec implementieren. statt Eigen matrix als m_matrix, ein vector
            getter und setter anpassen. getter und setter. vorher planen wie das interface funktionieren soll!
            getter(i,j) sollen weiterhin existieren, so müssen indist, conc, etc. nicht geändert werden.
            vector soll dann i->(a,b,c,d,..) mappen (t,tau,tau2,tau3,...). für G1 und G2 2D, für den rest nD
    TODO: loop implmenentieren
    TODO: testen, gucken ob sinnvoll.
            was geändert werden muss:
            - hier wie Gn gespeichert wird
            - wie Gn ausgegeben wird. indist,conc etc. können ihre syntax dank getter und setter behalten, aber 
                die outputfunktionen für Gn müssen angepasst werden um die weiteren indices zu berücksichtigen.
                dabei können einfach alle elemente ordered nach index ausgegeben werden, da der index eines elements nach
                t,tau,tau2,... sortiert ist, wobei t,tau,tau2,... von links angefangen inkrementiert sein sollen.


    TODO after: loop rekursiv umschreiben; nur nohc calculate_g_n haben, die für alles nehmen
                G1 wird dann wie jetzt als G2 mit einheitsmatritzen berechnet. 
                G2 ist dann die erste richtige korrelationsfunktion
                G(n>2) dann rekursiv. bei G2 sollte die rekursion schon nach einer iteration beendet sein,
                und das finale summieren in einen erwartungswert passieren. die indices sind dann nur 2D, wobei die
                indices für G(n>2) dann nD sind.
                wenn das funktioniert, mit vorheriger G3 vergleichen.
            hier stellt sich mir die frage, ob das überhaupt richtig ist lmao. weil so G1 ein spezialfall von G2 ist, aber alle höheren Gn sind 
            nur verschachtelte G2
*/

}