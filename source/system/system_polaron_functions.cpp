#include "system/system.h"

/**
 * @brief Calculates the Phonon Correlation Kerbal Phi(tau) = int_0^inf dw J(w) / w^2 * [coth(hbar w/(2 k_b T)) * cos(w t) - i*sin(w t)]
 *
 * @param t Current Time
 * @return Scalar value Phi(t)
 */
Scalar System::dgl_phonons_phi( const double t ) {
    Scalar integral = 0;
    double stepsize = 0.01 * parameters.p_phonon_wcutoff;
    for ( double w = stepsize; w < 10 * parameters.p_phonon_wcutoff; w += stepsize ) {
        double J;
        if ( parameters.p_phonon_qd_ae == 0.0 )
            J = parameters.p_phonon_alpha * w * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff );
        else
            J = w * parameters.hbar * std::pow( parameters.p_phonon_qd_de * std::exp( -w * w * parameters.p_phonon_qd_ae * parameters.p_phonon_qd_ae / ( 4. * parameters.p_phonon_qd_cs * parameters.p_phonon_qd_cs ) ) - parameters.p_phonon_qd_dh * std::exp( -w * w * parameters.p_phonon_qd_ae / parameters.p_phonon_qd_ratio * parameters.p_phonon_qd_ae / parameters.p_phonon_qd_ratio / ( 4. * parameters.p_phonon_qd_cs * parameters.p_phonon_qd_cs ) ), 2. ) / ( 4. * 3.1415 * 3.1415 * parameters.p_phonon_qd_rho * std::pow( parameters.p_phonon_qd_cs, 5. ) );
        integral += stepsize * J * ( std::cos( w * t ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) - 1.0i * std::sin( w * t ) );
    }
    return integral;
}

/**
 * @brief Initializes the Polaron Frame Functions by precalculating the Phi(tau) function and the corresponding Green functions
 *
 */
void System::initialize_polaron_frame_functions() {
    if ( parameters.p_phonon_T >= 0 ) {
        Log::L2( "[System-Polaron-Frame] Initializing Polaron Frame Functions.\n" );
        // Initialize Phi(tau)
        double tau = 0.0;
        double last = 1.0;
        double first = std::abs( dgl_phonons_phi( 0.0 ) );
        while ( parameters.p_phonon_tcutoff < 0 ? true : ( tau < parameters.p_phonon_tcutoff ) ) {
            phi_vector[tau] = dgl_phonons_phi( tau );
            dgl_phonons_greenf_matrix( tau, 'g' );
            if ( parameters.p_phonon_tcutoff < 0 and phi_vector.size() > 1 ) {
                double current = std::abs( phi_vector[tau] );
                if ( std::abs( 1.0 - current / last ) < 1E-2 or std::abs( 1.0 - current / first ) < 1E-3 ) {
                    parameters.p_phonon_tcutoff = tau;
                    Log::L2( "[System-Polaron-Frame] Polaron t-cutoff was automatically determined to t_cutoff = {}\n", parameters.p_phonon_tcutoff );
                    break;
                }
                last = current;
            }
            tau += parameters.numerics_subiterator_stepsize; // Use iterator stepsize here to allow for dt = -1 (variable timestep)
        }

        // Output Phonon Functions
        FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
        fmt::print( fp_phonons, "t\treal(phi(t))\timag(phi(t))\treal(g_u(t))\timag(g_u(t))\treal(g_g(t))\timag(g_g(t))\n" );
        for ( double t = parameters.t_start; t <= parameters.p_phonon_tcutoff; t += parameters.t_step ) {
            auto greenu = dgl_phonons_greenf( t, 'u' );
            auto greeng = dgl_phonons_greenf( t, 'g' );
            fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", t, std::real( phi_vector[t] ), std::imag( phi_vector[t] ), std::real( greenu ), std::imag( greenu ), std::real( greeng ), std::imag( greeng ) );
        }
        std::fclose( fp_phonons );
        if ( parameters.output_coefficients ) {
            // fp_phonons = std::fopen( ( parameters.subfolder + "phonons_lb.txt" ).c_str(), "w" );
            // fmt::print( fp_phonons, "t\tL_a_+\tL_a_-\tL_c_+\tL_c_-\n" );
            // for ( double t = parameters.t_start; t < parameters.t_end; t += parameters.t_step ) {
            //     fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\n", t, dgl_phonons_lindblad_coefficients( t, 'L', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'L', -1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', -1.0 ) );
            // }
            // std::fclose( fp_phonons );
        }
        Log::L2( "[System-Polaron-Frame] Done.\n" );
    }
}

/**
 * @brief The Polaron Transformed Operators are calculated by evaluating the Time Transformed X_g and X_u. This Runge Function is used to calculate this time transformation using the von-Neumann euqation instead of
 * numerically evaluating the tim transformation by e.g. using the numerical matrix exponential. The noteworthy difference here is to also include the explicit time dependency of the transformed operators.
 *
 * @param chi Current Chi(t)
 * @param t Current Time
 * @return Sparse: Applied Runge function
 */
Sparse System::dgl_phonons_rungefunc( const Sparse &chi, const double t ) {
    // TODO Maybe? Cache explicit times?
    double chirpcorrection = chirp.size() > 0 ? ( chirp.back().get( t ) + t * ( chirp.back().get( t ) - parameters.scaleVariable( chirp.back().derivative( t ), parameters.scale_value ) ) ) : 0;
    Sparse explicit_time = Sparse( chi.rows(), chi.cols() );
    for ( auto &[mode, param] : operatorMatrices.el_transitions ) {
        if ( param.direction == -1 )
            continue;
        explicit_time += ( param.energy + chirpcorrection ) * param.projector;
    }
    for ( auto &[mode, param] : parameters.input_photonic ) {
        for ( auto transition : param.string_v["CoupledTo"] ) {
            std::reverse( transition.begin(), transition.end() );
            explicit_time += 1.0i * ( operatorMatrices.el_transitions[transition].energy + chirpcorrection - operatorMatrices.ph_states[mode].energy ) * operatorMatrices.el_transitions[transition].projector * operatorMatrices.ph_transitions[mode + "b"].projector;
        }
    }
    int p = 0;
    int pp = 0;
    for ( auto &[mode, param] : parameters.input_pulse ) {
        for ( auto transition : param.string_v["CoupledTo"] ) {
            std::reverse( transition.begin(), transition.end() );
            explicit_time += 1.0i * ( pulse[p].get( t ) + 1.0i * parameters.scaleVariable( pulse[p].derivative( t, parameters.t_step ), parameters.scale_value ) ) * operatorMatrices.polaron_pulse_factors_explicit_time[pp++];
        }
        p++;
    }
    explicit_time = parameters.scaleVariable( explicit_time, 1.0 / parameters.scale_value );

    Sparse hamilton = dgl_getHamilton( t );

    return -1.0i * dgl_kommutator( hamilton, chi ) + explicit_time.cwiseProduct( QDLC::Matrix::sparse_projector( chi ) );
}

/**
 * @brief Transforms the Polaron Chi(t) to either X_g = Chi(t) + H.c. or X_u = i*Chi(t) - H.c.
 *
 * @param chi Current Chi(t)
 * @param mode Either "g" or "u"
 * @return Sparse: X_u/g
 */
Sparse System::dgl_phonons_chiToX( const Sparse &chi, const double t, const char mode ) {
    Sparse adjoint = chi.adjoint();
    if ( mode == 'g' ) {
        return chi + adjoint;
    }
    return 1.0i * ( chi - adjoint );
    // Sparse ret;
    // if ( mode == 'g' ) {
    //     // QD-Cavity
    //     Sparse adjoint = operatorMatrices.polaron_factors[0].adjoint();
    //     ret = operatorMatrices.polaron_factors[0] + adjoint;
    //     // QD-Pulse
    //     for ( int p = 1; p < parameters.input_pulse.size() + 1; p++ ) {
    //         Sparse adjoint = ( operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t ) ).adjoint();
    //         ret -= -( operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t ) + adjoint );
    //     }
    // } else {
    //     // QD-Cavity
    //     Sparse adjoint = operatorMatrices.polaron_factors[0].adjoint();
    //     ret = -1.i * ( -operatorMatrices.polaron_factors[0] + adjoint );
    //     // QD-Pulse
    //     for ( int p = 1; p < parameters.input_pulse.size() + 1; p++ ) {
    //         Sparse adjoint = ( operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t ) ).adjoint();
    //         ret -= 1.i * ( -operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t ) + adjoint );
    //     }
    // }
    // return dgl_timetrafo( ret, t );
}

/**
 * @brief Calculates the Polaron Green Function. This function is currently not cached.
 *
 * @param t Current Time
 * @param mode Either "g" or "u"
 * @return Scalar: G_u/g(t)
 */
Scalar System::dgl_phonons_greenf( double t, const char mode ) {
    if ( phi_vector.count( t ) == 0 )
        phi_vector[t] = dgl_phonons_phi( t );
    auto phi = phi_vector[t];
    if ( mode == 'g' ) {
        return parameters.p_phonon_b * parameters.p_phonon_b * ( std::cosh( phi ) - 1.0 );
    }
    return parameters.p_phonon_b * parameters.p_phonon_b * std::sinh( phi );
}

/**
 * @brief Calculates the Polaron Green Function as a Matrix considering the different coupling Rates for the Electronic Levels. This function is cached.
 *
 * @param t Current Time
 * @param mode Either "g" or "u"
 * @return Sparse&: Reference to cached G_u/g(t) in Matrix form
 */
Sparse &System::dgl_phonons_greenf_matrix( double t, const char mode ) {
    if ( not phi_vector.contains( t ) ) {
        phi_vector[t] = dgl_phonons_phi( t ); // std::cout << "Time " << t << " is not in phi vecotr\n";
    }
    if ( operatorMatrices.phi_vector_matrix_cache_g.contains( t ) ) {
        if ( mode == 'g' )
            return operatorMatrices.phi_vector_matrix_cache_g[t];
        return operatorMatrices.phi_vector_matrix_cache_u[t];
    } else {
        operatorMatrices.phi_vector_matrix_cache_g[t] = Sparse( operatorMatrices.polaron_phonon_coupling_matrix.cols(), operatorMatrices.polaron_phonon_coupling_matrix.rows() );
        operatorMatrices.phi_vector_matrix_cache_u[t] = Sparse( operatorMatrices.polaron_phonon_coupling_matrix.cols(), operatorMatrices.polaron_phonon_coupling_matrix.rows() );
    }
    auto phi = dgl_phonons_phi( t ); // phi_vector[t];
    auto &g = operatorMatrices.phi_vector_matrix_cache_g[t];
    auto &u = operatorMatrices.phi_vector_matrix_cache_u[t];
    // IDK why coefficient wise (unaryExpr) doesnt work lol. EDIT: because im stoopid. This function doesnt get evaluated often tho so its fine.
    for ( int k = 0; k < operatorMatrices.polaron_phonon_coupling_matrix.outerSize(); ++k ) {
        for ( Sparse::InnerIterator it( operatorMatrices.polaron_phonon_coupling_matrix, k ); it; ++it ) {
            auto row = it.row();
            auto col = it.col();
            auto val = it.value();
            g.coeffRef( row, col ) = ( std::cosh( phi * val ) - 1.0 );
            u.coeffRef( row, col ) = std::sinh( phi * val );
        }
    }
    Log::L3( "[System-Polaron-Frame] Calculated Phi Phonon Matrix cache for tau = {} -> Phi(tau) = {}\n", t, phi );
    if ( mode == 'g' )
        return operatorMatrices.phi_vector_matrix_cache_g[t];
    return operatorMatrices.phi_vector_matrix_cache_u[t];
}

/**
 * @brief Calculates to Lindblad Coefficients for the analytical solution of the Polaron Contributions.
 *
 * @param t Current Time
 * @param energy Transition Energy
 * @param coupling Cavity Coupling
 * @param pulse Current Pulse values
 * @param mode Mode
 * @param sign +1 or -1
 * @return double: L(t)
 */
double System::dgl_phonons_lindblad_coefficients( double t, double energy, double coupling, Scalar pulse, const char mode, const double sign ) {
    if ( t == 0.0 )
        t = parameters.t_step * 0.01;
    if ( phi_vector.count( t ) == 0 ) // TODO: get function for phi?
        phi_vector[t] = dgl_phonons_phi( t );
    double ret = 0;
    if ( mode == 'L' ) {
        double bpulsesquared = std::pow( std::abs( parameters.p_phonon_b * pulse ), 2.0 );

        double nu = std::sqrt( bpulsesquared + energy * energy );
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += parameters.t_step ) {
            Scalar f = ( energy * energy * std::cos( nu * tau ) + bpulsesquared ) / std::pow( nu, 2.0 );
            ret += std::real( ( std::cosh( phi_vector[t] ) - 1.0 ) * f + std::sinh( phi_vector[t] ) * std::cos( nu * tau ) ) - sign * std::imag( ( std::exp( phi_vector[t] ) - 1.0 ) * energy * std::sin( nu * tau ) / nu );
            i++;
        }
        ret *= 2.0 * bpulsesquared * parameters.t_step;
    } else if ( mode == 'C' ) {
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += parameters.t_step ) {
            ret += std::real( std::exp( 1.0i * sign * energy * tau ) * ( std::exp( phi_vector[t] ) - 1.0 ) );
            i++;
        }
        ret *= parameters.p_phonon_b * parameters.p_phonon_b * coupling * coupling * parameters.t_step;
    }
    return ret;
}

/**
 * @brief Calculates the Polaron Chi(t)
 *
 * @param t Current Time
 * @return Sparse: Interaction Picture Chi(t) in Matrix Form
 */
Sparse System::dgl_phonons_chi( const double t ) {
    // QD-Cavity
    Sparse ret = operatorMatrices.polaron_factors[0];
    // QD-Pulse
    for ( int p = 1; p < parameters.input_pulse.size() + 1; p++ ) {
        ret += operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t );
    }
    return dgl_timetrafo( ret, t );
}

/**
 * @brief Evaluates the Transformation U(t,tau) Chi(t) U(t,tau)^+ = Chi(t,tau) from tau' = 0 to tau' = tau
 *
 * @param chi_tau Current Chi(t)
 * @param t Current Time
 * @param tau Current Delay
 * @return Sparse: Transformed Chi(t) = Chi(t,tau)
 */
std::vector<QDLC::SaveState> System::dgl_phonons_calculate_transformation( Sparse &chi_tau, double t, double tau ) {
    // Backwards Integral
    if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_BACKWARDS_INTEGRAL ) {
        return QDLC::Numerics::calculate_definite_integral_vec( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, parameters.numerics_rk_tol, parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() );
    }
    // Matrix Exponential
    else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_TRANSFORMATION_MATRIX ) {
        // auto H = dgl_getHamilton( t );
        Sparse H = operatorMatrices.H_used + dgl_pulse( t ); // + dgl_chirp( t );
        std::vector<QDLC::SaveState> ret( std::ceil( tau / parameters.numerics_subiterator_stepsize ), QDLC::SaveState() );
#pragma omp parallel for ordered schedule( dynamic ) num_threads( parameters.numerics_phonons_maximum_threads )
        for ( int i = 0; i < ret.size(); i++ ) {
            double dtau = 1.0 * i * parameters.numerics_subiterator_stepsize;
            Sparse U = ( Dense( -1.0i * H * dtau ).exp() ).sparseView();
            // ret[i] = { dgl_timetrafo( ( U * chi_tau * U.adjoint() ).eval(), t - dtau ), t };
            ret[i] = { ( U * chi_tau * U.adjoint() ).eval(), t };
        }
        return ret;
    }
    // Backwards Integral if Threshold is met, else neglect transformation
    else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_MIXED ) {
        double threshold = 0;
        for ( auto &p : pulse )
            threshold += std::abs( p.get( t ) / p.maximum );                                   // TODO: + photon zahl, wenn photon number > 0.1 oder so dann auch. für große kopplungen gibts sonst starke abweichungen. vil. number*g draufaddieren.
        if ( threshold > 1E-4 || ( chirp.size() > 0 && chirp.back().derivative( t ) != 0 ) ) { // TODO: threshold als parameter
            return QDLC::Numerics::calculate_definite_integral_vec( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, parameters.numerics_rk_tol, parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() );
        }
    }
    return { { chi_tau, t } };
}

// TODO: im RK45 fall chis einfach aus t-richtung cachen, für tau richtung dann interpolieren. d.h. wenn ich t einsetze und t zwischen t0,t1 liegt _(t0,t1 schon gespeichert) dann interpoliere einfach linear zwischen den beiden. mal testen.
//  t-richtung sollte eigentlich schon genug infos berechnen für chi. lineare interpolation sollte legit sein. NICHT neu integrieren für tau (auch im RK4 oder RK5 fall nicht, da hat man eh schon alle gecached)
//  im tau-fall eifnach direkt auslesen und ggf interpolieren. dazu interpolate_single funktion schreiben, in solver interpolator auch die nehmen. inlined.
// FIX: einfach immer mit tstep integrieren, auf schon für direction, dann immer interpolieren. lel. die rhos könnte der sich hier auch lokal linear interpolieren für ohne markov
// NEXT: das hier immer nur mit tStepIterator rechnen, dann interpolieren für returnvalue.
/**
 * @brief Calculates the Time Transformed Chi for the Polaron Master Equation. This function either evaluates the PME or returns an interpolated Chi if allowed.
 *
 * @param rho Current Rho Matrix
 * @param t Current Time
 * @param past_rhos All past Rhos. If the general Markov Approximation is disabled, this vector is used instead of the simple Rho
 * @return Sparse
 */
Sparse System::dgl_phonons_pmeq( const Sparse &rho, const double t, const std::vector<QDLC::SaveState> &past_rhos ) {
    track_getcoefficient_calcattempt++;
    // All Contributions will (finally) be reduced onto this return value matrix
    Sparse ret = Sparse( parameters.maxStates, parameters.maxStates );
    // Most precise approximation used. Calculate polaron fram Chi by integrating backwards from t to t-tau.
    if ( parameters.numerics_phonon_approximation_order != PHONON_APPROXIMATION_LINDBLAD_FULL ) {
        // Calculate the initial Chi(t)
        Sparse chi = dgl_phonons_chi( t );
        Sparse chi_adjoint = chi.adjoint();
        // Calculate generic X_u and X_g from the initial Chi(t)
        Sparse XGT = dgl_timetrafo( chi + chi_adjoint, t );            // dgl_phonons_chiToX( chi, t, 'u' );
        Sparse XUT = dgl_timetrafo( 1.0i * ( chi - chi_adjoint ), t ); // dgl_phonons_chiToX( chi, t, 'g' );
        //// QD-Cavity
        // Sparse adj = operatorMatrices.polaron_factors[0].adjoint();
        // Sparse XGT = operatorMatrices.polaron_factors[0] + adj;
        // Sparse XUT = 1.i * operatorMatrices.polaron_factors[0] - 1.i * adj;
        //// QD-Pulse
        // for ( int p = 1; p < parameters.input_pulse.size() + 1; p++ ) {
        //     Sparse adj = operatorMatrices.polaron_factors[p].adjoint();
        //     Scalar pls = pulse[p - 1].get( t );
        //     XGT += operatorMatrices.polaron_factors[p] * pls + adj * std::conj( pls );
        //     XUT += 1.i * operatorMatrices.polaron_factors[p] * pls - 1.i * std::conj( pls ) * adj;
        // }
        // XGT = dgl_timetrafo( XGT, t );
        // XUT = dgl_timetrafo( XUT, t );
        int _taumax = (int)std::min( parameters.p_phonon_tcutoff / parameters.numerics_subiterator_stepsize, t / parameters.numerics_subiterator_stepsize ); // TODO: remove and replace with while (tau < tcutoff)
        // Log::L2( "Tau max = {}, stepsize = {}\n", _taumax, parameters.numerics_subiterator_stepsize );
        //  int _taumax = (int)std::floor( parameters.p_phonon_tcutoff / parameters.t_step ); // TODO: remove and replace with while (tau < tcutoff)
        //   Use markov approximation
        if ( parameters.numerics_phonon_approximation_markov1 ) {
            // Temporary variables
            Sparse chi_tau_back_u, chi_tau_back_g;
            // Look if we already calculated coefficient sum for this specific t-value or any value for t > t-value, then interpolate to that specific matrix
            // std::pair<double, double> time_pair = std::make_pair( t, 0.0 );

            // Index was found, chi(t-tau) sum used from saved vector
            auto map_index_to = savedCoefficients.begin();
            if ( parameters.numerics_use_saved_coefficients and savedCoefficients.count( t ) > 0 ) {
                Log::L3( "[System-Polaron-Frame] Thread #{} - Found {}\n", omp_get_thread_num(), t );
                auto &coeff = savedCoefficients[t][0.0];
                chi_tau_back_u = coeff.mat1;
                chi_tau_back_g = coeff.mat2;
                track_getcoefficient_read++;
            } else if ( parameters.numerics_use_saved_coefficients and savedCoefficients.size() > 2 and savedCoefficients.rbegin()->first > t ) {
                Log::L3( "[System-Polaron-Frame] Element {} should be in here as max is {}, trying to find it...\n", t, savedCoefficients.rbegin()->first );
                QDLC::SaveStateTau *min;
                QDLC::SaveStateTau *max;
                for ( auto &[nt, other] : savedCoefficients ) {
                    max = &other.begin()->second;
                    if ( nt > t )
                        break;
                    min = max;
                }

                Log::L3( "[System-Polaron-Frame] Found Phonon index! Interpolating from {} - {} to {}\n", min->t, max->t, t );
                chi_tau_back_u = min->mat1 + ( t - min->t ) / ( max->t - min->t ) * max->mat1;
                chi_tau_back_g = min->mat2 + ( t - min->t ) / ( max->t - min->t ) * max->mat2;
                track_getcoefficient_read_interpolated++;
            } else {
                Log::L3( "[System-Polaron-Frame] {} (Re)Calculating {}\n", omp_get_thread_num(), t );
                // Index was not found, (re)calculate chi(t-tau) sum
                // Initialize temporary matrices to zero for threads to write to
                std::vector<Sparse> threadmap_u( parameters.numerics_phonons_maximum_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
                std::vector<Sparse> threadmap_g( parameters.numerics_phonons_maximum_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
                // Calculate backwards integral and sum it into threadmaps. Threadmaps will later be summed into one coefficient matrix.
                // for ( int cur_thread = 0; cur_thread < parameters.numerics_phonons_maximum_threads; cur_thread++ ) {
                //     for ( int cur = 0; cur < _taumax; cur += parameters.numerics_phonons_maximum_threads ) {
                // auto X_transformed_g = dgl_phonons_calculate_transformation( XGT, t, _taumax * parameters.numerics_subiterator_stepsize );
                // auto X_transformed_u = dgl_phonons_calculate_transformation( XUT, t, _taumax * parameters.numerics_subiterator_stepsize );
                auto chis_transformed = dgl_phonons_calculate_transformation( chi, t, parameters.numerics_subiterator_stepsize * _taumax );
                //  std::cout << "Starting XUT for t =" << t << " \n"
                //            << Dense( XUT ).format( operatorMatrices.output_format ) << std::endl;
#pragma omp parallel for ordered schedule( dynamic ) shared( savedCoefficients ) num_threads( parameters.numerics_phonons_maximum_threads )
                for ( int _tau = 0; _tau < _taumax; _tau++ ) {
                    double tau = ( 1.0 * _tau ) * parameters.numerics_subiterator_stepsize;
                    // Sparse chi_tau = dgl_phonons_chi( t ); // dgl_phonons_chi( t - tau );
                    //  Log::L3( "[System-Polaron-Frame] Thread #{} - _tau = {}, size of vec = {}, Chi Index = {}\n", omp_get_thread_num(), _tau, chis_transformed_u.size(), std::min<int>( chis_transformed_u.size() - 1, _tau ) );
                    //  auto &X_tau_back_u = X_transformed_u.at( std::min<int>( X_transformed_u.size() - 1, _tau ) ).mat;
                    // auto &X_tau_back_g = X_transformed_g.at( std::min<int>( X_transformed_g.size() - 1, _tau ) ).mat;
                    auto &chi_tau_back = chis_transformed.at( std::min<int>( chis_transformed.size() - 1, _tau ) ).mat;
                    Sparse chi_tau_back_adjoint = chi_tau_back.adjoint();
                    auto X_tau_back_g = chi_tau_back + chi_tau_back_adjoint;           // dgl_phonons_chiToX( chi_tau_back, 'g' );
                    auto X_tau_back_u = 1.i * ( chi_tau_back - chi_tau_back_adjoint ); // dgl_phonons_chiToX( chi_tau_back, 'u' );
                    auto thread = omp_get_thread_num();
                    threadmap_u[thread] += X_tau_back_u.cwiseProduct( dgl_phonons_greenf_matrix( tau, 'u' ) ); // dgl_phonons_greenf( tau, 'u' ) * X_tau_back_u;
                    threadmap_g[thread] += X_tau_back_g.cwiseProduct( dgl_phonons_greenf_matrix( tau, 'g' ) ); // dgl_phonons_greenf( tau, 'g' ) * X_tau_back_g;
                }
// Sum all contributions from threadmaps into one coefficient
#pragma omp parallel sections
                {
#pragma omp section
                    chi_tau_back_u = std::accumulate( threadmap_u.begin(), threadmap_u.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
#pragma omp section
                    chi_tau_back_g = std::accumulate( threadmap_g.begin(), threadmap_g.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
                }
                // Save coefficients
                if ( parameters.numerics_use_saved_coefficients or parameters.numerics_force_caching ) {
#pragma omp master
                    savedCoefficients[t][0.0] = QDLC::SaveStateTau( chi_tau_back_u, chi_tau_back_g, t, 0 );
                    track_getcoefficient_write++;
                }
                track_getcoefficient_calculate++;
            }
            // Calculate phonon contributions from (saved/calculated) coefficients and rho(t)
            Sparse integrant = parameters.numerics_subiterator_stepsize * ( dgl_kommutator( XUT, chi_tau_back_u * rho ) + dgl_kommutator( XGT, chi_tau_back_g * rho ) );
            Sparse adjoint = integrant.adjoint();
            Log::L3( "[System-Polaron-Frame] Thread #{} - Adding ret value for {}\n", omp_get_thread_num(), t );
            ret -= integrant + adjoint;
        } else {
            std::vector<Sparse> threadmap_u( parameters.numerics_phonons_maximum_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
            // Dont use markov approximation; Full integral
            auto chis_transformed = dgl_phonons_calculate_transformation( chi, t, parameters.t_step * _taumax );
#pragma omp parallel for ordered schedule( dynamic ) shared( savedCoefficients ) num_threads( parameters.numerics_phonons_maximum_threads )
            for ( int _tau = 0; _tau < _taumax; _tau++ ) {
                int rho_index = std::max( 0, (int)past_rhos.size() - 1 - _tau ); // TODO: Das hier ist mit var timesteps der rhos auch grütze!
                double tau = ( 1.0 * _tau ) * parameters.t_step;

                // Index not found, or no saving is used. Recalculate chi(t-tau).
                auto &chi_tau_back = chis_transformed.at( std::min<int>( chis_transformed.size() - 1, _tau ) ).mat;
                Sparse chi_tau_back_adjoint = chi_tau_back.adjoint();
                auto chi_tau_back_g = chi_tau_back + chi_tau_back_adjoint;           // dgl_phonons_chiToX( chi_tau_back, 'u' );
                auto chi_tau_back_u = 1.i * ( chi_tau_back - chi_tau_back_adjoint ); // dgl_phonons_chiToX( chi_tau_back, 'g' );
                // If saving is used, save current chi(t-tau).
                if ( parameters.numerics_use_saved_coefficients or parameters.numerics_force_caching ) {
#pragma omp critical
                    savedCoefficients[t][tau] = QDLC::SaveStateTau( chi_tau_back_u, chi_tau_back_g, t, tau );
                }
                Sparse integrant = dgl_phonons_greenf_matrix( tau, 'u' ).cwiseProduct( dgl_kommutator( XUT, chi_tau_back_u * rho ) ); // FIXME: mit rk45 oder var gitter passt das rho hier nicht mehr! --> rho interpolieren!
                integrant += dgl_phonons_greenf_matrix( tau, 'g' ).cwiseProduct( dgl_kommutator( XGT, chi_tau_back_g * rho ) );       // FIXME: mit rk45 oder var gitter passt das rho hier nicht mehr! --> rho interpolieren!
                Sparse adjoint = integrant.adjoint();
                auto thread = omp_get_thread_num();
                threadmap_u.at( thread ) += ( integrant + adjoint ) * parameters.t_step;
            }

            ret -= std::accumulate( threadmap_u.begin(), threadmap_u.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
        }
    } else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_LINDBLAD_FULL ) {
        double chirpcorrection = chirp.size() > 0 ? chirp.back().get( t ) : 0;
        for ( auto &[name, mat] : parameters.input_pulse ) {
            for ( int p = 0; p < mat.string_v["CoupledTo"].size(); p++ ) {
                auto &mode = mat.string_v["CoupledTo"][p];
                auto &transition = operatorMatrices.el_transitions[mode].hilbert;
                auto mode_transposed = mode;
                std::reverse( mode_transposed.begin(), mode_transposed.end() );
                auto &transition_transposed = operatorMatrices.el_transitions[mode_transposed].hilbert;
                // std::cout << p << " m1 = " << mode << ", m2 = " << mode_transposed << std::endl;
                auto delta_E = operatorMatrices.el_transitions[mode].energy - mat.numerical_v["Frequency"][p] + chirpcorrection;
                // std::cout << "Pulsed Transition " << mode << " -> coeeff = " << dgl_phonons_lindblad_coefficients( t, delta_E, 0.0, pulse[p].get( t ), 'L', operatorMatrices.el_transitions[mode].direction ) << std::endl;
                ret += dgl_phonons_lindblad_coefficients( t, delta_E, 0.0, pulse[p].get( t ), 'L', -1 ) * dgl_lindblad( rho, transition, transition_transposed );
                ret += dgl_phonons_lindblad_coefficients( t, delta_E, 0.0, pulse[p].get( t ), 'L', 1 ) * dgl_lindblad( rho, transition_transposed, transition );
            }
        }
        for ( auto &[name, mat] : parameters.input_photonic ) {
            for ( auto &mode : mat.string_v["CoupledTo"] ) {
                int c = 0;
                auto &transition = operatorMatrices.el_transitions[mode].hilbert;
                auto mode_transposed = mode;
                std::reverse( mode_transposed.begin(), mode_transposed.end() );
                auto &transition_transposed = operatorMatrices.el_transitions[mode_transposed].hilbert;
                auto &optical_transition = operatorMatrices.ph_transitions[name + "b"].hilbert;
                auto &optical_transition_transposed = operatorMatrices.ph_transitions[name + "bd"].hilbert;
                auto delta_E = operatorMatrices.el_transitions[mode].energy - mat.numerical["Frequency"] + chirpcorrection;
                // std::cout << "Transition " << mode << "-" << name + "bd (" << dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c++] * parameters.p_omega_coupling, 0.0, 'C', -1 ) << "),   " << mode_transposed << "-" << name + "b (" << dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c++] * parameters.p_omega_coupling, 0.0, 'C', 1 ) << ")" << std::endl;
                ret += dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c++] * parameters.p_omega_coupling, 0.0, 'C', -1 ) * dgl_lindblad( rho, transition * optical_transition_transposed, transition_transposed * optical_transition );
                ret += dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c++] * parameters.p_omega_coupling, 0.0, 'C', 1 ) * dgl_lindblad( rho, transition_transposed * optical_transition, transition * optical_transition_transposed );
            }
        }
    }
    Log::L3( "[System-Polaron-Frame] Thread #{} - Returning for {}\n", omp_get_thread_num(), t );
    return ret;
}