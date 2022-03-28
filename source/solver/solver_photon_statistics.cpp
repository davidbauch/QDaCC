#include "solver/solver_ode.h"
#include <cmath>
#include <complex>
//#include <specfunc.h>

bool QDLC::Numerics::ODESolver::calculate_advanced_photon_statistics( System &s ) {
    // Calculate Spectra
    auto &spectrum_s = s.parameters.input_correlation["Spectrum"];
    for ( int i = 0; i < spectrum_s.string_v["Modes"].size(); i++ ) {
        const auto &[s_creator, s_annihilator] = get_operator_strings( spectrum_s.string_v["Modes"][i] );
        calculate_spectrum( s, s_creator, s_annihilator, spectrum_s.numerical_v["Center"][i], spectrum_s.numerical_v["Range"][i], spectrum_s.numerical_v["resW"][i], spectrum_s.string_v["Normalize"][i] == "True" );
    }
    // Calculate Indist
    auto &indist_s = s.parameters.input_correlation["Indist"];
    for ( int i = 0; i < indist_s.string_v["Modes"].size(); i++ ) {
        const auto &[s_creator, s_annihilator] = get_operator_strings( indist_s.string_v["Modes"][i] );
        calculate_indistinguishability( s, s_creator, s_annihilator );
    }
    // Calculate Conc
    auto &conc_s = s.parameters.input_correlation["Conc"];
    for ( auto &modes : conc_s.string_v["Modes"] ) {
        std::vector<std::string> s_creator, s_annihilator;
        for ( auto &mode : QDLC::String::splitline( modes, '-' ) ) {
            const auto &[ss_creator, ss_annihilator] = get_operator_strings( mode );
            s_creator.emplace_back( ss_creator );
            s_annihilator.emplace_back( ss_annihilator );
        }
        calculate_concurrence( s, s_creator[0], s_annihilator[0], s_creator[1], s_annihilator[1] );
    }
    // Calculate Raman
    auto &raman_s = s.parameters.input_correlation["Raman"];
    for ( int i = 0; i < raman_s.string_v["ElMode1"].size(); i++ ) {
        calculate_raman_population( s, raman_s.string_v["ElMode1"][i], raman_s.string_v["ElMode2"][i], raman_s.string_v["OpMode"][i], raman_s.string_v["PMode"][i] );
    }
    // Calculate G1/G2 functions
    auto &gs_s = s.parameters.input_correlation["GFunc"];
    for ( int i = 0; i < gs_s.string_v["Modes"].size(); i++ ) {
        double t_step = ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? s.parameters.t_step_pathint : s.parameters.t_step );
        /// TODO : in funktion
        auto modes = gs_s.string_v["Modes"][i];
        int order = std::abs( gs_s.numerical_v["Order"][i] );
        const auto &[s_creator, s_annihilator] = get_operator_strings( modes );
        std::string purpose = order == 1 ? get_operators_purpose( { s_creator, s_annihilator }, 1 ) : get_operators_purpose( { s_creator, s_annihilator, s_creator, s_annihilator }, 2 );
        Sparse creator, annihilator;
        if ( order == 1 ) {
            auto [a, b] = calculate_g1( s, s_creator, s_annihilator, purpose );
            creator = std::move( a );
            annihilator = std::move( b );
        } else {
            auto [a, b, _discard1, _discard2] = calculate_g2( s, s_creator, s_annihilator, s_creator, s_annihilator, purpose );
            creator = std::move( a );
            annihilator = std::move( b );
        }
        // Directly output corresponding matrix here so G1/2 functions calculated by other function calls are not output if they are not demanded.
        auto &gmat = cache[purpose];
        auto &gmat_time = cache[purpose + "_time"];
        // G2(t,tau)
        if ( gs_s.numerical_v["Integrated"][i] == 0 || gs_s.numerical_v["Integrated"][i] == 2 ) {
            Log::L2( "Saving G{} function matrix to {}_m.txt...\n", order, purpose );
            FILE *f_gfunc = std::fopen( ( s.parameters.subfolder + purpose + "_m.txt" ).c_str(), "w" );
            fmt::print( f_gfunc, "Time\tTau\tAbs\tReal\tImag\n" );
            for ( int k = 0; k < gmat.rows(); k++ ) {
                double t_t = std::real( gmat_time( k, 0 ) );
                for ( int l = 0; l < gmat.cols(); l++ ) {
                    double t_tau = std::imag( gmat_time( k, l ) ) - std::real( gmat_time( k, l ) );
                    fmt::print( f_gfunc, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t_t, t_tau, std::abs( gmat( k, l ) ), std::real( gmat( k, l ) ), std::imag( gmat( k, l ) ) );
                }
                fmt::print( f_gfunc, "\n" );
            }
            std::fclose( f_gfunc );
        }
        // G2(0)
        std::vector<Scalar> topv, g2ofzero;
        for ( int i = 0; i < std::min<int>( gmat.cols(), savedStates.size() ); i++ ) {
            g2ofzero.emplace_back( 0 );
            topv.emplace_back( 0 );
        }
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
        for ( int i = 0; i < std::min<int>( gmat.cols(), savedStates.size() ); i++ ) {
            for ( int j = 0; j < gmat.cols() - i; j++ ) {
                topv[i] += gmat( i, j );
            }
        }
        Scalar topsumv = 0;
        Scalar bottomsumv = 0;
        for ( int k = 0; k < topv.size(); k++ ) {
            double t_t = std::real( gmat_time( k, 0 ) );
            int t = rho_index_map[t_t];
            topsumv += topv[k];
            bottomsumv += s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( t ), creator * annihilator, t_t );
            g2ofzero[k] = 2.0 * topsumv / std::pow( bottomsumv, 2.0 );
        }
        // G2(t,0) and G2(tau)
        if ( gs_s.numerical_v["Integrated"][i] == 1 || gs_s.numerical_v["Integrated"][i] == 2 ) {
            Log::L2( "[PhotonStatistics] Saving G{} integrated function to {}.txt...\n", order, purpose );
            FILE *f_gfunc = std::fopen( ( s.parameters.subfolder + purpose + ".txt" ).c_str(), "w" );
            fmt::print( f_gfunc, "Time\tAbs(g{0}(tau))\tReal(g{0}(tau))\tImag(g{0}(tau))\tAbs(g{0}(t,0))\tReal(g{0}(t,0))\tImag(g{0}(t,0))\tAbs(g{0}(0))\tReal(g{0}(0))\tImag(g{0}(0))\n", order );
            for ( int l = 0; l < topv.size(); l++ ) { // gmat.cols()
                Scalar g2oftau = 0;
                for ( int k = 0; k < gmat.rows(); k++ ) {
                    g2oftau += gmat( k, l ) * Numerics::get_tdelta( gmat_time, l, k );
                }
                double t_tau = std::imag( gmat_time( 0, l ) );
                size_t tau_index = rho_index_map[t_tau];
                Scalar g2oft = s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( tau_index ), creator * creator * annihilator * annihilator, t_tau ); // / std::pow( s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( l ), creator * annihilator, getTimeAt( l ) ), 2.0 );
                fmt::print( f_gfunc, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t_tau, std::abs( g2oftau ), std::real( g2oftau ), std::imag( g2oftau ), std::abs( g2oft ), std::real( g2oft ), std::imag( g2oft ), std::abs( g2ofzero[l] ), std::real( g2ofzero[l] ), std::imag( g2ofzero[l] ) );
            }
            std::fclose( f_gfunc );
        }
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    // Calculate Conc
    auto &wigner_s = s.parameters.input_correlation["Wigner"];
    for ( int i = 0; i < wigner_s.string_v["Modes"].size(); i++ ) {
        calculate_wigner( s, wigner_s.string_v["Modes"][i], wigner_s.numerical_v["X"][i], wigner_s.numerical_v["Y"][i], wigner_s.numerical_v["Res"][i], wigner_s.numerical_v["Skip"][i] );
    }

    // Output Spectra and Rest in seperate Files
    for ( auto &[mode, data] : to_output["Spectrum"] ) {
        Log::L2( "[PhotonStatistics] Saving Emission Spectrum to spectrum_" + mode + ".txt...\n" );
        FILE *f_spectrum = std::fopen( ( s.parameters.subfolder + "spectrum_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_spectrum, "Omega\t{}\n", mode );
        for ( int i = 0; i < to_output["Spectrum"][mode].size(); i++ ) {
            fmt::print( f_spectrum, "{:.8e}\t{:.8e}\n", std::real( to_output["Spectrum_frequency"][mode][i] ), std::real( to_output["Spectrum"][mode][i] ) );
        }
        std::fclose( f_spectrum );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output["Indist"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Indistinguishability and Visibility to indist_" + mode + ".txt...\n" );
        FILE *f_indist = std::fopen( ( s.parameters.subfolder + "indist_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_indist, "Time\tIndist_{}\tVisibility_{}\n", mode, mode );
        for ( int i = 0; i < to_output["Indist"][mode].size(); i++ ) {
            fmt::print( f_indist, "{:.8e}\t{:.8e}\t{:.8e}\n", std::real( to_output["Indist"]["Time"][i] ), std::real( to_output["Indist"][mode][i] ), std::real( to_output["Visibility"][mode][i] ) );
        }
        std::fclose( f_indist );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output["Conc"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Concurrence to conc_" + mode + ".txt...\n" );
        FILE *f_indist = std::fopen( ( s.parameters.subfolder + "conc_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_indist, "Time\t{0}\t{0}_simple\t{0}(g2(0))\t{0}_simple(g2(0))\n", mode );
        // fmt::print( f_indist, "Time\t{0}\t{0}(g2(0))\n", mode );
        for ( int i = 0; i < to_output["Conc"][mode].size(); i++ ) {
            fmt::print( f_indist, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", std::real( to_output["Conc"]["Time"][i] ), std::real( to_output["Conc"][mode][i] ), std::real( to_output["Conc_simple"][mode][i] ), std::real( to_output["Conc_g2zero"][mode][i] ), std::real( to_output["Conc_g2zero_simple"][mode][i] ) );
            // fmt::print( f_indist, "{:.8e}\t{:.8e}\t{:.8e}\n", std::real( to_output["Conc"]["Time"][i] ), std::real( to_output["Conc"][mode][i] ), std::real( to_output["Conc_g2zero"][mode][i] ) );
        }
        std::fclose( f_indist );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output["Raman"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Raman Population to raman_" + mode + ".txt...\n" );
        FILE *f_raman = std::fopen( ( s.parameters.subfolder + "raman_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_raman, "Time\t{0}\tEM({0})\n", mode );
        for ( int i = 0; i < to_output["Raman"][mode].size(); i++ ) {
            fmt::print( f_raman, "{:.8e}\t{:.8e}\t{:.8e}\n", getTimeAt( i ), std::real( to_output["Raman"][mode][i] ), std::real( to_output["RamanEmProb"][mode][i] ) );
        }
        std::fclose( f_raman );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output_m["TwoPMat"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Two-Photon Matrix to twopmat_" + mode + ".txt...\n" );
        FILE *f_twophot = std::fopen( ( s.parameters.subfolder + "twopmat_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_twophot, "Time\t" );
        std::vector<std::string> modes = { "11", "12", "21", "22" };
        for ( int k = 0; k < 4; k++ ) {
            for ( int l = 0; l < 4; l++ ) {
                fmt::print( f_twophot, "Re({}{})\t", modes[k], modes[l] );
            }
        }
        for ( int k = 0; k < 4; k++ ) {
            for ( int l = 0; l < 4; l++ ) {
                fmt::print( f_twophot, "Im({}{})\t", modes[k], modes[l] );
            }
        }
        fmt::print( f_twophot, "\n" );
        for ( int i = 0; i < to_output_m["TwoPMat"][mode].size(); i++ ) {
            fmt::print( f_twophot, "{:.8e}\t", std::real( to_output["Conc"]["Time"][i] ) );
            auto &mat = to_output_m["TwoPMat"][mode][i];
            for ( int k = 0; k < 4; k++ ) {
                for ( int l = 0; l < 4; l++ ) {
                    fmt::print( f_twophot, "{:.8e}\t", std::real( mat( k, l ) ) );
                }
            }
            for ( int k = 0; k < 4; k++ ) {
                for ( int l = 0; l < 4; l++ ) {
                    fmt::print( f_twophot, "{:.8e}\t", std::imag( mat( k, l ) ) );
                }
            }

            fmt::print( f_twophot, "\n" );
        }
        std::fclose( f_twophot );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output_m["Wigner"] ) {
        auto &wigner_s = s.parameters.input_correlation["Wigner"];
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Wigner function to wigner_" + mode + ".txt...\n" );
        FILE *f_wigner = std::fopen( ( s.parameters.subfolder + "wigner_" + mode + ".txt" ).c_str(), "w" );
        if ( mode.starts_with( "rho_" ) ) {
            // Gather base:
            fmt::print( f_wigner, "Time\t" );
            std::string smode = mode.substr( 4 );
            int base = s.operatorMatrices.el_states.count( smode ) != 0 ? s.operatorMatrices.el_states[smode].base : s.operatorMatrices.ph_states[smode].base;
            if ( base == 0 ) {
                for ( auto &[name, dat] : s.parameters.input_electronic )
                    for ( auto &[name2, dat2] : s.parameters.input_electronic )
                        fmt::print( f_wigner, "Re(|{}><{}|)\t", name, name2 );
                for ( auto &[name, dat] : s.parameters.input_electronic )
                    for ( auto &[name2, dat2] : s.parameters.input_electronic )
                        fmt::print( f_wigner, "Im(|{}><{}|)\t", name, name2 );
            } else {
                for ( int i = 0; i < data[0].rows(); i++ )
                    for ( int j = 0; j < data[0].rows(); j++ )
                        fmt::print( f_wigner, "Re(|{}_{}><{}_{}|)\t", smode, i, smode, j );
                for ( int i = 0; i < data[0].rows(); i++ )
                    for ( int j = 0; j < data[0].rows(); j++ )
                        fmt::print( f_wigner, "Im(|{}_{}><{}_{}|)\t", smode, i, smode, j );
            }
            fmt::print( f_wigner, "\n" );
        } else {
            fmt::print( f_wigner, "Time\t{}\n", mode );
        }
        for ( int i = 0; i < data.size(); i++ ) {
            fmt::print( f_wigner, "{:.8e}\t", std::real( to_output["Wigner"]["Time"][i] ) );
            auto &currentwigner = data[i];
            for ( int k = 0; k < currentwigner.rows(); k++ ) {
                for ( int l = 0; l < currentwigner.cols(); l++ ) {
                    fmt::print( f_wigner, "{:.8e}\t", std::real( currentwigner( k, l ) ) );
                }
            }
            for ( int k = 0; k < currentwigner.rows(); k++ ) {
                for ( int l = 0; l < currentwigner.cols(); l++ ) {
                    fmt::print( f_wigner, "{:.8e}\t", std::imag( currentwigner( k, l ) ) );
                }
            }
            fmt::print( f_wigner, "\n" );
        }
        std::fclose( f_wigner );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    return true;
}

// electronic_transition1, electronic_transition2, optical_transition
bool QDLC::Numerics::ODESolver::calculate_raman_population( System &s, const std::string &source_transitions, const std::string &raman_transition, const std::string &optical_transition, const std::string &pulse_mode ) {
    // auto [m_electronic_transition1, m_electronic_transition2] = get_operators_matrices( s, source_transitions, );
    // if ( s.operatorMatrices.el_transitions.count( raman_transition ) == 0 ) {
    //     auto ket = s.operatorMatrices.el_states[raman_transition.substr( 0, 1 )].ket;
    //     auto bra = s.operatorMatrices.el_states[raman_transition.substr( 1, 1 )].bra;
    //     auto current = s.operatorMatrices.base_selfhilbert;
    //     current.front() = ket * bra;
    //     Sparse transition_hilbert = QDLC::Matrix::tensor( current ).sparseView();
    //     {
    //         auto &state = s.operatorMatrices.extra_transitions[raman_transition];
    //         state.ket = ket;
    //         state.bra = state.ket.transpose();
    //         state.self_hilbert = state.ket * state.bra;
    //         state.base = 0;
    //         state.hilbert = transition_hilbert;
    //     }
    // }
    // Sparse m_raman_transition = s.operatorMatrices.el_transitions.count( raman_transition ) == 0 ) ? s.operatorMatrices.extra_transitions[raman_transition].hilbert : s.operatorMatrices.el_transitions[raman_transition].hilbert;
    //// Sparse m_electronic_transition3 = m_electronic_transition1 * m_electronic_transition2;
    // auto [m_optical_transition, dump] = get_operators_matrices( s, optical_transition + "bd", optical_transition + "bd" );
    //// Energies
    // std::vector<double> e_source_transitions;
    //... = s.operatorMatrices.el_transitions[electronic_transition1].energy;
    // double e_electronic_transition2 = s.operatorMatrices.el_transitions[electronic_transition2].energy;
    // double e_optical_transition = s.operatorMatrices.ph_transitions[optical_transition + "bd"].energy;
    // int pulse_index = s.parameters.input_pulse[pulse_mode].numerical["PulseIndex"];
    //
    // Log::L2( "Using w_1 = {}, w_2 = {}, w_c = {}, pulse index = {}({})\n", e_electronic_transition1, e_electronic_transition2, e_optical_transition, pulse_mode, pulse_index );
    //
    // std::vector<Scalar> raman_pop;
    // std::vector<Scalar> raman_emission;
    // Scalar before = 0;
    // std::string fout = electronic_transition1 + electronic_transition2 + optical_transition + pulse_mode;
    //// Progress
    // ProgressBar progressbar = ProgressBar();
    // Timer &timer_r = Timers::create( "Raman (" + fout + ")" );
    // timer_r.start();
    //
    //// for ( size_t T = 0; T < savedStates.size(); T++ ) {
    ////     double t = savedStates.at( T ).t;
    ////     double dt = get_tdelta( savedStates, T );
    ////     double chirpcorrection = s.chirp.size() > 0 ? s.chirp.back().get( t ) : 0;
    ////     double w1 = e_electronic_transition1 + chirpcorrection;
    ////     double w2 = e_electronic_transition2 + chirpcorrection;
    ////     double wc = e_optical_transition;
    ////     Scalar next = 0;
    ////
    ////    double sigma1 = s.parameters.p_omega_pure_dephasing + s.parameters.p_omega_decay;
    ////    double sigma2 = s.parameters.p_omega_pure_dephasing + 3. * s.parameters.p_omega_decay;
    ////    Sparse op = m_electronic_transition3 * m_optical_transition;
    ////    double A = std::exp( -s.parameters.p_omega_cavity_loss * dt );
    ////    Scalar B, R;
    ////    for ( size_t i = 0; i < savedStates.size(); i++ ) {
    ////        double tdd = savedStates.at( i ).t;
    ////        double dtau = get_tdelta( savedStates, i );
    ////        B = std::exp( -1.0i * ( w2 - wc - 0.5i * ( s.parameters.p_omega_cavity_loss + sigma2 ) ) * ( t - tdd ) ) - std::exp( -1.0i * ( w1 - wc - 0.5i * ( s.parameters.p_omega_cavity_loss + sigma1 ) ) * ( t - tdd ) );
    ////        R = s.dgl_expectationvalue<Sparse, Scalar>( savedStates.at( i ).mat, op, tdd ) * std::conj( s.pulse.at( pulse_index ).get( tdd ) );
    ////        next += B * R * dtau;
    ////    }
    ////    before = A * before + dt * next;
    ////    raman_pop.emplace_back( before );
    ////    raman_emission.emplace_back( T == 0 ? s.parameters.p_omega_cavity_loss * before * dt : raman_emission.back() + s.parameters.p_omega_cavity_loss * before * dt );
    ////    timer_r.iterate();
    ////    Timers::outputProgress( timer_r, progressbar, timer_r.getTotalIterationNumber(), savedStates.size(), "Raman (" + fout + "): " );
    ////}
    //
    // for ( size_t T = 0; T < savedStates.size(); T++ ) {
    //    Scalar integral = 0;
    //    double t = savedStates.at( T ).t;
    //    double dT = get_tdelta( savedStates, T );
    //    for ( size_t i = 0; i < T; i++ ) {
    //        double td = savedStates.at( td ).t;
    //        double dt = get_tdelta( savedStates, td );
    //        double chirpcorrection = s.chirp.size() > 0 ? s.chirp.back().get( td ) : 0;
    //        double w1 = e_electronic_transition1 + chirpcorrection;
    //        double w2 = e_electronic_transition2 + chirpcorrection;
    //        double wc = e_optical_transition;
    //
    //        double sigma1 = s.parameters.p_omega_pure_dephasing + s.parameters.p_omega_decay;
    //        double sigma2 = s.parameters.p_omega_pure_dephasing + 3. * s.parameters.p_omega_decay;
    //        Sparse op = m_raman_transition * m_optical_transition;
    //        for ( size_t j = 0; j < i; j++ ) {
    //            double tdd = savedStates.at( j ).t;
    //            double dtau = get_tdelta( savedStates, j );
    //            Scalar B = std::exp( -1.0i * ( w1 - wc - 0.5i * ( s.parameters.p_omega_cavity_loss + sigma1 ) ) * ( td - tdd ) ) - std::exp( -1.0i * ( w2 - wc - 0.5i * ( s.parameters.p_omega_cavity_loss + sigma2 ) ) * ( td - tdd ) );
    //            Scalar R = s.dgl_expectationvalue<Sparse, Scalar>( savedStates.at( j ).mat, op, tdd ) * s.pulse.at( pulse_index ).get( tdd );
    //            integral += 2.0i * s.parameters.p_omega_coupling * std::exp( -s.parameters.p_omega_cavity_loss * ( t - td ) ) * B * R * dtau * dt;
    //            // std::cout << s.dgl_expectationvalue<Sparse, Scalar>( savedStates.at( j ).mat, op, tdd ) << ", " << std::conj( s.pulse.at( pulse_index ).get( tdd ) ) << "," << dt << "," << dtau << "," << B << "," << R << "," << s.parameters.p_omega_coupling << "," << s.parameters.p_omega_cavity_loss << "," << std::exp( -s.parameters.p_omega_cavity_loss * ( T - td ) ) << " - " << integral << std::endl;
    //        }
    //    }
    //    raman_pop.emplace_back( std::real( integral ) );
    //    raman_emission.emplace_back( T == 0 ? s.parameters.p_omega_cavity_loss * integral * dT : raman_emission.back() + s.parameters.p_omega_cavity_loss * integral * dT );
    //    timer_r.iterate();
    //    Timers::outputProgress( timer_r, progressbar, timer_r.getTotalIterationNumber(), savedStates.size(), "Raman (" + fout + "): " );
    //}
    // timer_r.end();
    // Timers::outputProgress( timer_r, progressbar, timer_r.getTotalIterationNumber(), savedStates.size(), "Raman (" + fout + ")", Timers::PROGRESS_FORCE_OUTPUT );
    // to_output["Raman"][fout] = raman_pop;
    // to_output["RamanEmProb"][fout] = raman_emission;
    return false;
}