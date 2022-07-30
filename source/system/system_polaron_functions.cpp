#include "system/system.h"

Scalar System::dgl_phonons_J( const double w ) {
    if ( parameters.p_phonon_qd_ae == 0.0 )
        return parameters.p_phonon_alpha * w * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff );
    return w * parameters.hbar * std::pow( parameters.p_phonon_qd_de * std::exp( -w * w * parameters.p_phonon_qd_ae * parameters.p_phonon_qd_ae / ( 4. * parameters.p_phonon_qd_cs * parameters.p_phonon_qd_cs ) ) - parameters.p_phonon_qd_dh * std::exp( -w * w * parameters.p_phonon_qd_ae / parameters.p_phonon_qd_ratio * parameters.p_phonon_qd_ae / parameters.p_phonon_qd_ratio / ( 4. * parameters.p_phonon_qd_cs * parameters.p_phonon_qd_cs ) ), 2. ) / ( 4. * 3.1415 * 3.1415 * parameters.p_phonon_qd_rho * std::pow( parameters.p_phonon_qd_cs, 5. ) );
}

Scalar System::dgl_phonons_phi( const double tau ) {
    Scalar integral = 0;
    for ( double w = parameters.p_phonon_wcutoffdelta; w < 10.0 * parameters.p_phonon_wcutoff; w += parameters.p_phonon_wcutoffdelta ) {
        Scalar J = dgl_phonons_J( w );
        integral += parameters.p_phonon_wcutoffdelta * J * ( std::cos( w * tau ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) - 1.0i * std::sin( w * tau ) );
    }
    return integral;
}

void System::initialize_polaron_frame_functions() {
    if ( parameters.p_phonon_T >= 0 ) {
        Log::L2( "[System-PME] Initializing Polaron Frame Functions.\n" );
        // Initialize Phi(tau)
        double tau = 0.0;
        double last = 1.0;
        double first = std::abs( dgl_phonons_phi( 0.0 ) );
        while ( parameters.p_phonon_tcutoff < 0 ? true : ( tau < parameters.p_phonon_tcutoff ) ) {
            phi_vector[tau] = dgl_phonons_phi( tau );
            dgl_phonons_greenf_matrix( tau, 'g' );
            // FIXME: doesnt really work for stepsizes < 100fs
            if ( parameters.p_phonon_tcutoff < 0 and phi_vector.size() > 1 ) {
                double current = std::real( phi_vector[tau] );
                if ( std::abs( 1.0 - last / current ) < 1E-3 ) {
                    parameters.p_phonon_tcutoff = tau;
                    Log::L2( "[System-PME] Polaron t-cutoff was automatically determined to t_cutoff = {}\n", parameters.p_phonon_tcutoff );
                    break;
                }
                last = current;
            }
            tau += parameters.numerics_subiterator_stepsize; // Use iterator stepsize here to allow for dt = -1 (variable timestep)
        }

        // Output Phonon Functions.
        if ( parameters.output_dict.contains( "greenf" ) ) {
            auto &file = FileOutput::add_file( "phonon_greenfunctions" );
            file << fmt::format( "t\treal(phi(t))\timag(phi(t))\treal(g_u(t))\timag(g_u(t))\treal(g_g(t))\timag(g_g(t))\n" );
            for ( double t = parameters.t_start; t <= parameters.p_phonon_tcutoff; t += parameters.numerics_subiterator_stepsize ) {
                auto greenu = dgl_phonons_greenf( t, 'u' );
                auto greeng = dgl_phonons_greenf( t, 'g' );
                file << fmt::format( "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", t, std::real( phi_vector[t] ), std::imag( phi_vector[t] ), std::real( greenu ), std::imag( greenu ), std::real( greeng ), std::imag( greeng ) );
            }
        }
        if ( parameters.output_dict.contains( "phononJ" ) ) {
            auto &file = FileOutput::add_file( "phonon_spectral" );
            file << fmt::format( "omega\tJ(omega)\n" );
            for ( double w = parameters.p_phonon_wcutoffdelta; w < 10.0 * parameters.p_phonon_wcutoff; w += parameters.p_phonon_wcutoffdelta ) {
                file << fmt::format( "{}\t{}\n", w, std::real( dgl_phonons_J( w ) ) );
            }
        }
        if ( parameters.output_dict.contains( "phononcoefficients" ) ) {
            // fp_phonons = std::fopen( ( parameters.working_directory + "phonons_lb.txt" ).c_str(), "w" );
            // fmt::print( fp_phonons, "t\tL_a_+\tL_a_-\tL_c_+\tL_c_-\n" );
            // for ( double t = parameters.t_start; t < parameters.t_end; t += parameters.t_step ) {
            //     fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\n", t, dgl_phonons_lindblad_coefficients( t, 'L', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'L', -1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', -1.0 ) );
            // }
            // std::fclose( fp_phonons );
        }
        Log::L2( "[System-PME] Done.\n" );
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
    double chirpcorrection = not chirp.empty() ? ( chirp.back().get( t ) + t * ( chirp.back().get( t ) - parameters.scaleVariable( chirp.back().derivative( t ), parameters.scale_value ) ) ) : 0;
    auto explicit_time = Sparse( chi.rows(), chi.cols() );
    for ( const auto &[mode, param] : operatorMatrices.el_transitions ) {
        if ( param.direction == -1 )
            continue;
        explicit_time += ( param.energy + chirpcorrection ) * param.projector;
    }
    for ( auto &[mode, param] : parameters.input_photonic ) {
        for ( const auto &transition : param.string_v["CoupledTo"] ) {
            auto transition_transposed = operatorMatrices.el_transitions[transition].name_transposed;
            explicit_time += 1.0i * ( operatorMatrices.el_transitions[transition_transposed].energy + chirpcorrection - operatorMatrices.ph_states[mode].energy ) * operatorMatrices.el_transitions[transition_transposed].projector * operatorMatrices.ph_transitions[mode + "b"].projector;
        }
    }
    int p = 0;
    int pp = 0;
    for ( auto &[mode, param] : parameters.input_pulse ) {
        for ( const auto &transition : param.string_v["CoupledTo"] ) {
            explicit_time += 1.0i * ( pulse[p].get( t ) + 1.0i * parameters.scaleVariable( pulse[p].derivative( t, parameters.t_step ), parameters.scale_value ) ) * operatorMatrices.polaron_pulse_factors_explicit_time[pp++];
        }
        p++;
    }
    explicit_time = parameters.scaleVariable( explicit_time, 1.0 / parameters.scale_value );

    Sparse hamilton = dgl_get_hamilton( t );

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
    return chi;
    // Sparse adjoint = chi.adjoint();
    // if ( mode == 'g' ) {
    //     return chi + adjoint;
    // }
    // return 1.0i * ( chi - adjoint );
    //  Sparse ret;
    //  if ( mode == 'g' ) {
    //      // QD-Cavity
    //      Sparse adjoint = operatorMatrices.polaron_factors[0].adjoint();
    //      ret = operatorMatrices.polaron_factors[0] + adjoint;
    //      // QD-Pulse
    //      for ( int p = 1; p < parameters.input_pulse.size() + 1; p++ ) {
    //          Sparse adjoint = ( operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t ) ).adjoint();
    //          ret -= -( operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t ) + adjoint );
    //      }
    //  } else {
    //      // QD-Cavity
    //      Sparse adjoint = operatorMatrices.polaron_factors[0].adjoint();
    //      ret = -1.i * ( -operatorMatrices.polaron_factors[0] + adjoint );
    //      // QD-Pulse
    //      for ( int p = 1; p < parameters.input_pulse.size() + 1; p++ ) {
    //          Sparse adjoint = ( operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t ) ).adjoint();
    //          ret -= 1.i * ( -operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t ) + adjoint );
    //      }
    //  }
    //  return dgl_timetrafo( ret, t );
}

/**
 * @brief Calculates the Polaron Green Function. This function is currently not cached.
 *
 * @param t Current Time
 * @param mode Either "g" or "u"
 * @return Scalar: G_u/g(t)
 */
Scalar System::dgl_phonons_greenf( double tau, const char mode ) {
    if ( not phi_vector.contains( tau ) )
        phi_vector[tau] = dgl_phonons_phi( tau );
    auto phi = phi_vector[tau];
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
Sparse &System::dgl_phonons_greenf_matrix( double tau, const char mode ) {
    if ( operatorMatrices.pme_greenfunction_matrix_cache_g.contains( tau ) ) {
        if ( mode == 'g' )
            return operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
        return operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
    } else {
        operatorMatrices.pme_greenfunction_matrix_cache_g[tau] = Sparse( operatorMatrices.polaron_phonon_coupling_matrix.cols(), operatorMatrices.polaron_phonon_coupling_matrix.rows() );
        operatorMatrices.pme_greenfunction_matrix_cache_u[tau] = Sparse( operatorMatrices.polaron_phonon_coupling_matrix.cols(), operatorMatrices.polaron_phonon_coupling_matrix.rows() );
    }
    if ( not phi_vector.contains( tau ) ) {
        phi_vector[tau] = dgl_phonons_phi( tau ); // std::cout << "Time " << t << " is not in phi vecotr\n";
    }
    auto phi = phi_vector[tau];
    auto &g = operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
    auto &u = operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
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
    Log::L3( "[System-PME] Calculated Phi Phonon Matrix cache for tau = {} -> Phi(tau) = {}\n", tau, phi );
    if ( mode == 'g' )
        return operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
    return operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
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
double System::dgl_phonons_lindblad_coefficients( double t, double energy, double coupling, Scalar current_pulse, const char mode, const double scaling, const double sign ) {
    Log::L3( "[System-PME]     Calculating Lindblad Rates for t = {}, integral border = {}, deltaE = {}, g = {}, Omega = {}, mode = {}, sign = {}, coupling scaling = {}, integral stepsize = {}\n", t, parameters.p_phonon_tcutoff, energy, coupling, current_pulse, mode, sign, scaling, parameters.numerics_subiterator_stepsize );
    if ( t == 0.0 )
        t = parameters.numerics_subiterator_stepsize * 0.01;
    double ret = 0;
    if ( mode == 'L' ) {
        double bpulsesquared = std::pow( std::abs( parameters.p_phonon_b * scaling * current_pulse ), 2.0 );
        double nu = std::sqrt( bpulsesquared + energy * energy );
        Log::L3( "[System-PME]         Nu = {}, Omega*<B>^2 = {}\n", nu, bpulsesquared );
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += parameters.numerics_subiterator_stepsize ) {
            if ( not phi_vector.contains( tau ) )
                phi_vector[tau] = dgl_phonons_phi( tau );
            auto phi = phi_vector[tau];
            Scalar f = ( energy * energy * std::cos( nu * tau ) + bpulsesquared ) / std::pow( nu, 2.0 );
            ret += std::real( ( std::cosh( phi ) - 1.0 ) * f + std::sinh( phi ) * std::cos( nu * tau ) ) - sign * std::imag( ( std::exp( phi ) - 1.0 ) * energy * std::sin( nu * tau ) / nu );
            i++;
        }
        ret *= 2.0 * bpulsesquared * parameters.numerics_subiterator_stepsize;
    } else if ( mode == 'C' ) {
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += parameters.numerics_subiterator_stepsize ) {
            if ( not phi_vector.contains( tau ) )
                phi_vector[tau] = dgl_phonons_phi( tau );
            auto phi = phi_vector[tau];
            ret += std::real( std::exp( 1.0i * sign * energy * tau ) * ( std::exp( phi ) - 1.0 ) );
            i++;
        }
        ret *= parameters.p_phonon_b * parameters.p_phonon_b * scaling * scaling * coupling * coupling * parameters.numerics_subiterator_stepsize;
    }
    if ( std::isnan( ret ) )
        ret = 0.0;
    Log::L3( "[System-PME]         Return value: {}\n", ret );
    return ret;
}

/**
 * @brief Calculates the Polaron Chi(t)
 *
 * @param t Current Time
 * @return Sparse: Interaction Picture Chi(t) in Matrix Form
 */
Sparse System::dgl_phonons_chi( const double t ) {
    Log::L3( "[System-PME] Generating Chi({})...\n", t );
    // QD-Cavity
    Sparse ret = operatorMatrices.polaron_factors[0];
    // QD-Pulse
    for ( int i = 0; i < pulse.size(); i++ ) {
        auto p = pulse[i].get( t );
        if ( parameters.numerics_use_rwa )
            ret += operatorMatrices.pulse_mat[2 * i + 1] * p;
        else
            ret += operatorMatrices.pulse_mat[2 * i + 1] * ( p + std::conj( p ) );
        Log::L3( "[System-PME]     Added Pulse with Omega(t) = {}\n", pulse[i].get( t ) );
    }
    Log::L3( "[System-PME]     Done, final untransformed Chi:\n{}\n", Dense( ret ).format( operatorMatrices.output_format ) );
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
Sparse System::dgl_phonons_calculate_transformation( double t, double tau ) {
    // Backwards Integral
    if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_BACKWARDS_INTEGRAL ) {
        // TODO
        // return QDLC::Numerics::calculate_definite_integral_vec( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, std::get<1>( parameters.numerics_rk_tol.front() ), parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() );
    }
    // Matrix Exponential
    else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_TRANSFORMATION_MATRIX ) {
        auto chi_tau = dgl_phonons_chi( t - tau ); //TODO maybe cache these.
        Sparse U = ( Dense( -1.0i * dgl_get_hamilton( t ) * tau ).exp() ).sparseView();
        return ( U * chi_tau * U.adjoint() );
    }
    return dgl_phonons_chi( t - tau );
    // Backwards Integral if Threshold is met, else neglect transformation
    // else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_MIXED ) {
    //    double threshold = 0;
    //    for ( auto &p : pulse )
    //        threshold += std::abs( p.get( t ) / p.maximum );                                     // TODO: + photon zahl, wenn photon number > 0.1 oder so dann auch. für große kopplungen gibts sonst starke abweichungen. vil. number*g draufaddieren.
    //    if ( threshold > 1E-4 || ( not chirp.empty() and chirp.back().derivative( t ) != 0 ) ) { // TODO: threshold als parameter
    //        return QDLC::Numerics::calculate_definite_integral_vec( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, std::get<1>( parameters.numerics_rk_tol.front() ), parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() );
    //    }
    //}
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
    auto ret = Sparse( parameters.maxStates, parameters.maxStates );
    // Most precise approximation used. Calculate polaron fram Chi by integrating backwards from t to t-tau.
    if ( parameters.numerics_phonon_approximation_order != PHONON_APPROXIMATION_LINDBLAD_FULL ) {
        // Calculate the initial Chi(t) which is in the interaction frame.
        Sparse chi = dgl_phonons_chi( t ); // dgl_timetrafo( dgl_phonons_chi( t ), t );
        Sparse chi_adjoint = chi.adjoint();
        // Calculate generic X_u and X_g from the initial Chi(t)
        Sparse XGT = chi + chi_adjoint;           // dgl_phonons_chiToX( chi, t, 'u' ); // X_G = Chi + H.c.
        Sparse XUT = 1.i * ( chi - chi_adjoint ); // dgl_phonons_chiToX( chi, t, 'g' ); // X_U = i(Chi - H.c.)
        Log::L3( "[System-PME] Calculating the PME contribution using the Polaron Transformation for t = {}. Initial Time Transformed Chi(t):\n{}\nXG(t):\n{}\nXU(t):\n{}\n", t, Dense( chi ).format( operatorMatrices.output_format ), Dense( XGT ).format( operatorMatrices.output_format ), Dense( XUT ).format( operatorMatrices.output_format ) );
        // Log::L2( "Tau max = {}, stepsize = {}\n", tau_max, parameters.numerics_subiterator_stepsize );
        // Use markov approximation
        if ( parameters.numerics_phonon_approximation_markov1 ) {
            // Temporary variables
            Sparse chi_tau_back_g;
            Sparse chi_tau_back_u;
            // Look if we already calculated coefficient sum for this specific t-value or any value for t > t-value, then interpolate to that specific matrix
            // std::pair<double, double> time_pair = std::make_pair( t, 0.0 );

            // Index was found, chi(t-tau) sum used from saved vector
            if ( parameters.numerics_use_saved_coefficients and savedCoefficients.contains( t ) ) {
                Log::L3( "[System-PME]     Thread #{} - Found {}\n", omp_get_thread_num(), t );
                const auto &coeff = savedCoefficients[t][0.0];
                chi_tau_back_u = coeff.mat1;
                chi_tau_back_g = coeff.mat2;
                track_getcoefficient_read++;
            }
            // Index wasn't found, try to interpolate
            else if ( parameters.numerics_use_saved_coefficients and savedCoefficients.size() > 2 and savedCoefficients.rbegin()->first > t ) {
                Log::L3( "[System-PME]     Element {} should be in here as max is {}, trying to find it...\n", t, savedCoefficients.rbegin()->first );
                // Hacky way to find [min,<t>,max], should work because this is an ordered map.
                const auto &greater_or_equal_than_t = savedCoefficients.lower_bound( t );
                const auto &smaller_than_t = std::prev( greater_or_equal_than_t );
                //  Interpolate
                Log::L3( "[System-PME]     Found Phonon index! Interpolating from {} - {} to {}\n", smaller_than_t->second.begin()->second.t, greater_or_equal_than_t->second.begin()->second.t, t );
                chi_tau_back_u = QDLC::Math::lerp( smaller_than_t->second.begin()->second.mat1, greater_or_equal_than_t->second.begin()->second.mat1, ( t - smaller_than_t->first ) / ( greater_or_equal_than_t->first - smaller_than_t->first ) );
                chi_tau_back_g = QDLC::Math::lerp( smaller_than_t->second.begin()->second.mat2, greater_or_equal_than_t->second.begin()->second.mat2, ( t - smaller_than_t->first ) / ( greater_or_equal_than_t->first - smaller_than_t->first ) );
                track_getcoefficient_read_interpolated++;
            }
            // Not using saved_coefficients or index wasn't found and interpolation wasn't succesfull, recalculating.
            else {
                Log::L3( "[System-PME]     Thread #{} - (Re)Calculating {} using {} steps in the PME integral\n", omp_get_thread_num(), t, tau_max );
                // Index was not found, (re)calculate chi(t-tau) sum
                // Initialize temporary matrices to zero for threads to write to
                std::vector<Sparse> threadmap_u( parameters.numerics_maximum_secondary_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
                std::vector<Sparse> threadmap_g( parameters.numerics_maximum_secondary_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
                // Calculate backwards integral and sum it into threadmaps. Threadmaps will later be summed into one coefficient matrix.
                int tau_max = std::min<int>( parameters.p_phonon_tcutoff/parameters.numerics_subiterator_stepsize, t/parameters.numerics_subiterator_stepsize );
#pragma omp parallel for schedule( dynamic ) num_threads( parameters.numerics_maximum_secondary_threads )
                for ( int tau_index = 0; tau_index < tau_max; tau_index++ ) {
                    double tau = parameters.numerics_subiterator_stepsize * tau_index;
                    Sparse chi_tau_back = dgl_phonons_calculate_transformation(t,tau);
                    Sparse chi_tau_back_adjoint = chi_tau_back.adjoint();
                    auto X_tau_back_g = chi_tau_back + chi_tau_back_adjoint;           // dgl_phonons_chiToX( chi_tau_back, 'g' );
                    auto X_tau_back_u = 1.i * ( chi_tau_back - chi_tau_back_adjoint ); // dgl_phonons_chiToX( chi_tau_back, 'u' );
                    auto thread = omp_get_thread_num();
                    threadmap_g[thread] += dgl_phonons_greenf( tau, 'g' ) * X_tau_back_g; //.cwiseProduct( dgl_phonons_greenf_matrix( tau, 'g' ) );
                    threadmap_u[thread] += dgl_phonons_greenf( tau, 'u' ) * X_tau_back_u; //.cwiseProduct( dgl_phonons_greenf_matrix( tau, 'u' ) );
                    Log::L3( "[System-PME]         Added Contributions for tau = {}. The sum of the contributions is {}, the phi value is {}\n", tau, threadmap_u[thread].sum() + threadmap_g[thread].sum(), phi_vector[tau] );
                }
// Sum all contributions from threadmaps into one coefficient
#pragma omp parallel sections
                {
#pragma omp section
                    {std::cout << omp_get_thread_num() << std::endl; chi_tau_back_u = std::accumulate( threadmap_u.begin(), threadmap_u.end(), Sparse( parameters.maxStates, parameters.maxStates ) ); }
#pragma omp section
                    {std::cout << omp_get_thread_num() << std::endl; chi_tau_back_g = std::accumulate( threadmap_g.begin(), threadmap_g.end(), Sparse( parameters.maxStates, parameters.maxStates ) ); }
                }
                // Save coefficients
                if ( parameters.numerics_enable_saving_coefficients ) {
                    savedCoefficients[t][0.0] = QDLC::SaveStateTau( chi_tau_back_u, chi_tau_back_g, t, 0 );
                    track_getcoefficient_write++;
                }
                track_getcoefficient_calculate++;
            }
            // Calculate phonon contributions from (saved/calculated) coefficients and rho(t)
            Sparse integrant = parameters.numerics_subiterator_stepsize * ( dgl_kommutator( XUT, chi_tau_back_u * rho ) + dgl_kommutator( XGT, chi_tau_back_g * rho ) );
            Sparse adjoint = integrant.adjoint();
            ret -= integrant + adjoint;
        } else {
            //TODO: FIXME
            //std::vector<Sparse> threadmap_u( parameters.numerics_maximum_secondary_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
            //// Dont use markov approximation; Full integral
            //auto chis_transformed = dgl_phonons_calculate_transformation( chi, t, parameters.t_step * tau_max );
            //#pragma omp parallel for ordered schedule( dynamic ) shared( savedCoefficients ) num_threads( parameters.numerics_maximum_secondary_threads )
            //for ( int tau_index = 0; tau_index < tau_max; tau_index++ ) {
            //    int rho_index = std::max( 0, (int)past_rhos.size() - 1 - tau_index ); // TODO: Das hier ist mit var timesteps der rhos auch grütze!
            //    double tau = ( 1.0 * tau_index ) * parameters.t_step;
//
            //    // Index not found, or no saving is used. Recalculate chi(t-tau).
            //    auto &chi_tau_back = chis_transformed.at( std::min<int>( chis_transformed.size() - 1, tau_index ) ).mat;
            //    Sparse chi_tau_back_adjoint = chi_tau_back.adjoint();
            //    auto chi_tau_back_g = chi_tau_back + chi_tau_back_adjoint;           // dgl_phonons_chiToX( chi_tau_back, 'u' );
            //    auto chi_tau_back_u = 1.i * ( chi_tau_back - chi_tau_back_adjoint ); // dgl_phonons_chiToX( chi_tau_back, 'g' );
            //    // If saving is used, save current chi(t-tau).
            //    if ( parameters.numerics_use_saved_coefficients or parameters.numerics_enable_saving_coefficients ) {
            //#pragma omp critical
            //        savedCoefficients[t][tau] = QDLC::SaveStateTau( chi_tau_back_u, chi_tau_back_g, t, tau );
            //    }
            //    Sparse integrant = dgl_phonons_greenf_matrix( tau, 'u' ).cwiseProduct( dgl_kommutator( XUT, chi_tau_back_u * rho ) ); // FIXME: mit rk45 oder var gitter passt das rho hier nicht mehr! --> rho interpolieren!
            //    integrant += dgl_phonons_greenf_matrix( tau, 'g' ).cwiseProduct( dgl_kommutator( XGT, chi_tau_back_g * rho ) );       // FIXME: mit rk45 oder var gitter passt das rho hier nicht mehr! --> rho interpolieren!
            //    Sparse adjoint = integrant.adjoint();
            //    auto thread = omp_get_thread_num();
            //    threadmap_u[thread] += ( integrant + adjoint ) * parameters.t_step;
            //}
            //ret -= std::accumulate( threadmap_u.begin(), threadmap_u.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
        }
        Log::L3( "[System-PME]     Thread #{} - Adding ret value for {}. Matrix is: \n{}\n", omp_get_thread_num(), t, Dense( ret ).format( operatorMatrices.output_format ) );
    } else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_LINDBLAD_FULL ) {
        Log::L3( "[System-PME] Calculating Lindblad-Type Phonon Rates for t = {}\n", t );
        double chirpcorrection = chirp.empty() ? 0.0 : chirp.back().get( t );
        for ( auto &[name, mat] : parameters.input_pulse ) {
            int p = 0;
            for ( int m = 0; m < mat.string_v["CoupledTo"].size(); m++ ) {
                const auto &mode = mat.string_v["CoupledTo"][m];
                const auto &mode_transposed = operatorMatrices.el_transitions[mode].name_transposed;
                const auto &transition = operatorMatrices.el_transitions[mode].hilbert;
                const auto &transition_transposed = operatorMatrices.el_transitions[mode_transposed].hilbert;
                // std::cout << p << " m1 = " << mode << ", m2 = " << mode_transposed << std::endl;
                auto delta_E = mat.numerical_v["Frequency"][m] - operatorMatrices.el_transitions[mode].energy + chirpcorrection;
                auto r1 = dgl_phonons_lindblad_coefficients( t, delta_E, 0.0, pulse[p].get( t ), 'L', 1.0, -1 );
                auto r2 = dgl_phonons_lindblad_coefficients( t, delta_E, 0.0, pulse[p].get( t ), 'L', 1.0, 1 );
                Log::L3( "[System-PME]     Pulse induced Phonon Transition rates: {} ({}), {} ({}) for delta_E = {}\n", mode, r1, mode_transposed, r2, delta_E );
                ret += r1 * dgl_lindblad( rho, transition, transition_transposed );
                ret += r2 * dgl_lindblad( rho, transition_transposed, transition );
            }
            p++;
        }
        for ( auto &[name, mat] : parameters.input_photonic ) {
            for ( auto &mode : mat.string_v["CoupledTo"] ) {
                int c = 0;
                const std::string upper_level = operatorMatrices.el_transitions[mode].to; // QDLC::String::split( mode, parameters.transition_delimiter ).back();
                const double phonon_scaling = parameters.input_electronic[upper_level].numerical["PhononCoupling"];
                const auto &mode_transposed = operatorMatrices.el_transitions[mode].name_transposed;
                const auto &transition = operatorMatrices.el_transitions[mode].hilbert;
                const auto &transition_transposed = operatorMatrices.el_transitions[mode_transposed].hilbert;
                const auto &optical_transition = operatorMatrices.ph_transitions[name + "b"].hilbert;
                const auto &optical_transition_transposed = operatorMatrices.ph_transitions[name + "bd"].hilbert;
                auto delta_E = mat.numerical["Energy"] - operatorMatrices.el_transitions[mode].energy + chirpcorrection;
                auto r1 = dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c] * parameters.p_omega_coupling, 0.0, 'C', phonon_scaling, -1 );
                auto r2 = dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c] * parameters.p_omega_coupling, 0.0, 'C', phonon_scaling, 1 );
                Log::L3( "[System-PME]     Cavity induced Phonon Transition rates: {}-{}bd ({}), {}-{}b ({})", mode, name, r1, mode_transposed, name, r2 );
                ret += r1 * dgl_lindblad( rho, transition * optical_transition_transposed, transition_transposed * optical_transition );
                ret += r2 * dgl_lindblad( rho, transition_transposed * optical_transition, transition * optical_transition_transposed );
                c++;
            }
        }
    }
    Log::L3( "[System-PME] Thread #{} - Returning PME Contribution for {}:\n{}\n", omp_get_thread_num(), t, Dense( ret ).format( operatorMatrices.output_format ) );
    return ret;
}