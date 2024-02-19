#include "solver/solver_ode.h"
#include <cmath>
#include <complex>
#include <ranges>
// #include <specfunc.h>

void _output_gfunc_matrix( QDACC::System &s, auto &f_gfunc, auto &gmat, bool use_triangular = true ) {
    // Build G function header.
    const auto dim = gmat.dim();
    f_gfunc << "Time\tTau";
    for ( size_t tau = 0; tau < gmat.getOrder() - 2; tau++ ) f_gfunc << "\tTau_" << tau;
    f_gfunc << "\tAbs\tReal\tImag\n";

    // Output G function matrix.
    for ( size_t i = 0; i < gmat.size(); i++ ) {
        // Get index vector corresponding to the flat index i
        auto index_vector = gmat.get_index( i );
        Scalar el;
        if ( use_triangular and not gmat.insideTriangular( index_vector ) ) continue;
        el = gmat.get( i );
        // Calculate times
        auto times = std::views::iota( (size_t)0, index_vector.size() ) | std::views::transform( [&]( size_t a ) { return s.getTimeOf( a, index_vector ); } );
        // Output Times
        for ( const double time : times ) f_gfunc << std::format( "{:.8e}\t", time );
        // Output Element
        f_gfunc << std::format( "{:.8e}\t{:.8e}\t{:.8e}", std::abs( el ), std::real( el ), std::imag( el ) ) << std::endl;
        // Force linebreak after each row. This is only required when using gnuplot.
        if ( i % dim == dim - 1 ) f_gfunc << std::endl;
    }
}

void _output_integrated_matrix( QDACC::System &s, auto &savedStates, auto &f_gfunc, auto &gmat, auto &creator, auto &annihilator, int order, std::string purpose ) {
    // Integral buffers; index:buffer. Initially filled with all zeros.
    auto num_buffers = std::max<int>( 2, order );
    const auto dim = gmat.dim();
    std::vector<std::vector<Scalar>> result_buffer( num_buffers, std::vector<Scalar>( dim, { 0.0, 0.0 } ) );

    // Integrate the Gn matrix elements.
    // auto threads = s.parameters.numerics_maximum_primary_threads;
    // #pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( threads )
    for ( size_t i = 0; i < gmat.size(); i++ ) {
        // Thread index
        auto thread = omp_get_thread_num();
        // Get index vector corresponding to the flat index i. Don't calculate times, we don't need them.
        auto index_vector = gmat.get_index( i );
        // Instead, calculate dts
        auto dts = std::views::iota( 0, int( index_vector.size() ) ) | std::views::transform( [&]( int a ) { return s.getDeltaTimeOf( a, index_vector ); } );
        // Convert dts into a std::vector
        std::vector<double> dts_vec( dts.begin(), dts.end() );
        // Calculate the product of all dts.
        double dt = std::accumulate( dts_vec.begin(), dts_vec.end(), 1.0, std::multiplies<double>() );
        // Calculate the sum reduction of the Gn matrix element over all indices.
        for ( size_t tau_index = 0; tau_index < index_vector.size(); tau_index++ ) {
            auto time_index = index_vector[tau_index];
            result_buffer[tau_index][time_index] += gmat.get( i ) * dt / dts_vec[tau_index];
        }
    }

    // Create fixed creator and annihilator chains (because std::pow isnt working for Eigen matrices)
    MatrixMain super_creator = creator;
    MatrixMain super_annihilator = annihilator;
    for ( int i = 1; i < order; i++ ) {
        super_creator = ( super_creator * creator ).eval();
        super_annihilator = ( super_annihilator * annihilator ).eval();
    }
    MatrixMain creator_annihilator = creator * annihilator;
    MatrixMain super_creator_annihilator = super_creator * super_annihilator;

    // Output file headers
    f_gfunc << "Time";
    for ( size_t tau_i = 0; tau_i < result_buffer.size(); tau_i++ ) f_gfunc << std::format( "\tRe(G^{0}(tau_{1}))\tIm(G^{0}(tau_{1}))", order, tau_i );
    f_gfunc << std::format( "\tRe(G^{0}(0))\tIm(G^{0}(0))\tRe(G^{0}_pop)\tIm(G^{0}_pop)\n", order );

    // Output elements
    Log::L2( "[PhotonStatistics] Outputting G{} integrated function to {}.txt, iterating over {}({}) states...\n", order, purpose, savedStates.size(), result_buffer.front().size() );
    for ( size_t t = 0; t < std::min( savedStates.size(), result_buffer.front().size() ); t++ ) {
        const double t_t = s.parameters.grid_values.at( t );
        const auto &rho = savedStates.at( t ).mat;
        // Gi of zero
        Scalar g_n_of_order = s.dgl_expectationvalue<MatrixMain>( rho, creator_annihilator, t_t );
        Scalar g_n_pop = std::pow( s.dgl_expectationvalue<MatrixMain>( rho, creator_annihilator, t_t ), order );
        f_gfunc << std::format( "{:.8}\t", t_t );
        for ( size_t tau_i = 0; tau_i < result_buffer.size(); tau_i++ ) {
            f_gfunc << std::format( "{:.8}\t{:.8}\t", std::real( result_buffer[tau_i][t] ), std::imag( result_buffer[tau_i][t] ) );
        }
        f_gfunc << std::format( "{:.8}\t{:.8}\t{:.8}\t{:.8}\n", std::real( g_n_of_order ), std::imag( g_n_of_order ), std::real( g_n_pop ), std::imag( g_n_pop ) );
    }
}

bool QDACC::Numerics::ODESolver::calculate_advanced_photon_statistics( System &s ) {
    // Calculate Spectra
    auto &all_spectra = s.parameters.input_correlation["Spectrum"];
    for ( auto &spectrum_s : all_spectra )
        for ( size_t i = 0; i < spectrum_s.string_v["Modes"].size(); i++ ) {
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, spectrum_s.string_v["Modes"][i] );
            calculate_spectrum( s, s_creator, s_annihilator, spectrum_s.property_set["Center"][i], spectrum_s.property_set["Range"][i], spectrum_s.property_set["resW"][i], spectrum_s.property_set["Order"][i],
                                spectrum_s.string_v["Normalize"][i] == "True" );
        }
    // Calculate Indist
    auto &all_indist = s.parameters.input_correlation["Indist"];
    for ( auto &indist_s : all_indist )
        for ( size_t i = 0; i < indist_s.string_v["Modes"].size(); i++ ) {
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, indist_s.string_v["Modes"][i] );
            calculate_indistinguishability( s, s_creator, s_annihilator );
        }
    // Calculate Conc
    auto &all_conc = s.parameters.input_correlation["Conc"];
    for ( auto &conc_s : all_conc )
        for ( auto i = 0; i < conc_s.string_v["Modes"].size(); i++ ) {
            const auto &modes = conc_s.string_v["Modes"][i];
            const std::map<std::string, int> orders = { { "full", 3 }, { "outin", 2 }, { "outer", 1 } };
            const auto order = conc_s.string_v["Order"][i];
            const int matrix_priority_evaluation = orders.contains( order ) ? orders.at( order ) : 3;
            std::vector<std::string> s_creator, s_annihilator;
            for ( auto &mode : QDACC::String::splitline( modes, '-' ) ) {
                const auto &[ss_creator, ss_annihilator] = get_operator_strings( s, mode );
                s_creator.emplace_back( ss_creator );
                s_annihilator.emplace_back( ss_annihilator );
            }
            // Spectra
            bool use_spec = not conc_s.property_set["Center"].empty();
            const double spec_center = use_spec ? conc_s.property_set["Center"].front().get() : -1.;
            const double spec_range = use_spec ? conc_s.property_set["Range"].front().get() : -1.;
            const double spec_resW = use_spec ? conc_s.property_set["resW"].front().get() : -1.;
            // Calculate TPM for modes
            calculate_static_two_photon_matrix( s, s_creator[0], s_annihilator[0], s_creator[1], s_annihilator[1], matrix_priority_evaluation, spec_center, spec_range, spec_resW );
            // Calculate Concurrence
            calculate_concurrence( s, s_creator[0], s_annihilator[0], s_creator[1], s_annihilator[1], conc_s.string_v["Method"][i] );
        }
    // Calculate Raman
    auto &all_raman = s.parameters.input_correlation["Raman"];
    for ( auto &raman_s : all_raman )
        for ( size_t i = 0; i < raman_s.string_v["ElMode1"].size(); i++ ) {
            calculate_raman_population( s, raman_s.string_v["ElMode1"][i], raman_s.string_v["ElMode2"][i], raman_s.string_v["OpMode"][i], raman_s.string_v["PMode"][i] );
        }
    // Output Both masks to files
    auto &gmat = cache.begin()->second;
    for ( const auto &[mode, current_detector_temporal_mask] : detector_temporal_mask ) {
        Log::L2( "[PhotonStatistics] Saving Detector Matrix to detector_temporal_mask_{}.txt...\n", mode );
        auto &f_detector = FileOutput::add_file( "detector_temporal_mask_" + mode );
        f_detector << "Time_index\ttau_index\tD(t)*D(t+tau)\n";
        for ( size_t k = 0; k < current_detector_temporal_mask.dim(); k++ ) {
            const auto t = s.getTimeOf( 0, { k, 0 } );
            for ( size_t l = 0; l < current_detector_temporal_mask.dim(); l++ ) {
                const auto tau = s.getTimeOf( 1, { k, l } );
                const auto val = std::real( current_detector_temporal_mask.get( k, l ) );
                f_detector << std::format( "{}\t{}\t{:.8e}\n", t, tau, val );
            }
            f_detector << "\n";
        }
        f_detector.close();
    }
    for ( auto &[mode, current_detector_frequency_mask] : detector_frequency_mask ) {
        Log::L2( "[PhotonStatistics] Saving Spectral Detector to detector_spectral_mask_{}.txt...\n", mode );
        auto &f_detector = FileOutput::add_file( "detector_spectral_mask_" + mode );
        f_detector << "Omega\tD(omega)\n";
        for ( const auto &[omega, amplitude, delta_omega] : current_detector_frequency_mask ) {
            f_detector << std::format( "{:.8e}\t{:.8e}\n", omega, amplitude );
        }
        f_detector.close();
    }

    // Calculate Time Bin Coherences
    auto &all_timbinfuncs = s.parameters.input_correlation["GFuncTimeBins"];
    Log::L2( "[PhotonStatistics] Calculating G2/3 Time Bin coherences...\n" );
    std::map<std::string, int> mode_output;
    for ( auto &gs_s : all_timbinfuncs ) {
        for ( size_t i = 0; i < gs_s.string_v["Modes"].size(); i++ ) {
            std::string purpose;
            auto modes = gs_s.string_v["Modes"][i];
            int order = std::abs( gs_s.property_set["Order"][i] );
            double start = std::abs( gs_s.property_set["Start"][i] );
            double bin_length = std::abs( gs_s.property_set["BinLength"][i] );
            int deph = gs_s.property_set["Deph"][i];
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, modes );
            const auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            Log::L2( "[PhotonStatistics] Operators are {} and {}.\n", s_creator, s_annihilator );
            switch ( order ) {
                case 2:
                    purpose = get_operators_purpose( { s_creator, s_creator, s_annihilator, s_annihilator }, "TimeBin-G" );
                    calculate_timebin_g2_correlations( s, s_creator, s_creator, s_annihilator, s_annihilator, purpose, start, bin_length, deph );
                    break;
                case 3:
                    purpose = get_operators_purpose( { s_creator, s_creator, s_creator, s_annihilator, s_annihilator, s_annihilator }, "TimeBin-G" );
                    calculate_timebin_g3_correlations( s, s_creator, s_creator, s_creator, s_annihilator, s_annihilator, s_annihilator, purpose, start, bin_length, deph );
                    break;
                default:
                    Log::Warning( "[PhotonStatistics] G^({}) function order not implemented!\n", order );
                    continue;
            }

            // We count the number of times a specific purpose has been output. TODO: move this into a function so 
            // The other correlation functions can use it too.
            //if (not mode_output.contains(purpose))
            //    mode_output[purpose] = 0;
            //else
            //    mode_output[purpose]++;
            
            std::vector<std::string> output_modes = { "EEEE", "EEEL", "EELE", "EELL", "ELEE", "ELEL", "ELLE", "ELLL", "LEEE", "LEEL", "LELE", "LELL", "LLEE", "LLEL", "LLLE", "LLLL" };
            if ( order == 3 )
                output_modes = {
                    "EEEEEE", "EEEEEL", "EEEELE", "EEEELL", "EEELLE", "EEELEE", "EEELEL", "EEELLL", "EELEEE", "EELEEL", "EELELE", "EELELL", "EELLLE", "EELLEE", "EELLEL", "EELLLL",
                    "ELEEEE", "ELEEEL", "ELEELE", "ELEELL", "ELELLE", "ELELEE", "ELELEL", "ELELLL", "ELLEEE", "ELLEEL", "ELLELE", "ELLELL", "ELLLLE", "ELLLEE", "ELLLEL", "ELLLLL",
                    "LLEEEE", "LLEEEL", "LLEELE", "LLEELL", "LLELLE", "LLELEE", "LLELEL", "LLELLL", "LEEEEE", "LEEEEL", "LEEELE", "LEEELL", "LEELLE", "LEELEE", "LEELEL", "LEELLL",
                    "LELEEE", "LELEEL", "LELELE", "LELELL", "LELLLE", "LELLEE", "LELLEL", "LELLLL", "LLLEEE", "LLLEEL", "LLLELE", "LLLELL", "LLLLLE", "LLLLEE", "LLLLEL", "LLLLLL",
                };
//#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
            for ( const auto &mode : output_modes ) {
                const auto inner_purpose = purpose + "_" + mode;
                auto &gmat = cache[inner_purpose];
                if ( gs_s.string_v["Integrated"][i] == "matrix" || gs_s.string_v["Integrated"][i] == "both" ) {
                    Log::L2( "[PhotonStatistics] Saving G2/3 time bin function matrix to {}_m.txt...\n", inner_purpose );
                    auto &f_gfunc = FileOutput::add_file( inner_purpose + "_m" );
                    _output_gfunc_matrix( s, f_gfunc, gmat, false /* no triangular cutoff */ );
                }
                // G(t, integral tau, ...), G( integral t, tau, ...), ..., G(0), Gpop(t). g(0) then is G(0)/Gpop
                if ( gs_s.string_v["Integrated"][i] == "time" || gs_s.string_v["Integrated"][i] == "both" ) {
                    Log::L2( "[PhotonStatistics] Saving G integrated time bin function to {}.txt...\n", inner_purpose );
                    auto &f_gfunc = FileOutput::add_file( inner_purpose );
                    _output_integrated_matrix( s, savedStates, f_gfunc, gmat, creator, annihilator, order, inner_purpose );
                }
            }

            if (order == 2)
                calculate_timebin_two_photon_matrix(s, purpose);
            else
                calculate_timebin_three_photon_matrix(s, purpose);

            calculate_timebin_generator_expectation_values(s, purpose);
            auto &timebin_exp = FileOutput::add_file( "generator_"+purpose );
            timebin_exp << "Time\tRe(G)\tIm(G)\n";
            auto& mat = to_output["GeneratorExpectationG"][purpose];
            for ( size_t k = 0; k < mat.size(); k++ ) {
                const auto t = s.getTimeOf( 0, { k, 0 } );
                const auto val = mat[k];
                timebin_exp << std::format( "{:.8e}\t{:.8e}\t{:.8e}\n", t, std::real( val ), std::imag( val ) );
            }
            timebin_exp.close();
        }
    }

    // Calculate G1/G2/G3 functions
    auto &all_gfuncs = s.parameters.input_correlation["GFunc"];
    Log::L2( "[PhotonStatistics] Calculating G1/G2/G3 functions...\n" );
    for ( auto &gs_s : all_gfuncs )
        for ( size_t i = 0; i < gs_s.string_v["Modes"].size(); i++ ) {
            auto modes = gs_s.string_v["Modes"][i];
            int order = std::abs( gs_s.property_set["Order"][i] );
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, modes );
            std::string purpose;
            const auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            switch ( order ) {
                case 1:
                    purpose = get_operators_purpose( { s_creator, s_annihilator } );
                    calculate_g1( s, s_creator, s_annihilator, purpose );
                    break;
                case 2:
                    purpose = get_operators_purpose( { s_creator, s_creator, s_annihilator, s_annihilator } );
                    calculate_g2( s, s_creator, s_creator, s_annihilator, s_annihilator, purpose );
                    break;
                case 3:
                    purpose = get_operators_purpose( { s_creator, s_creator, s_creator, s_annihilator, s_annihilator, s_annihilator } );
                    calculate_g3( s, s_creator, s_creator, s_creator, s_annihilator, s_annihilator, s_annihilator, purpose );
                    break;
                default:
                    Log::Warning( "[PhotonStatistics] G^({}) function order not implemented!\n", order );
                    continue;
            }

            // Directly output corresponding matrix here so G1/2/3 functions calculated by other function calls are not output if they are not demanded.
            auto &gmat = cache[purpose];
            // G^(n)(t,tau)
            if ( gs_s.string_v["Integrated"][i] == "matrix" || gs_s.string_v["Integrated"][i] == "both" ) {
                Log::L2( "[PhotonStatistics] Saving G{} function matrix to {}_m.txt...\n", order, purpose );
                auto &f_gfunc = FileOutput::add_file( purpose + "_m" );
                _output_gfunc_matrix( s, f_gfunc, gmat );
            }
            // G(t, integral tau, ...), G( integral t, tau, ...), ..., G(0), Gpop(t). g(0) then is G(0)/Gpop
            if ( gs_s.string_v["Integrated"][i] == "time" || gs_s.string_v["Integrated"][i] == "both" ) {
                Log::L2( "[PhotonStatistics] Saving G{} integrated function to {}.txt...\n", order, purpose );
                auto &f_gfunc = FileOutput::add_file( purpose );
                _output_integrated_matrix( s, savedStates, f_gfunc, gmat, creator, annihilator, order, purpose );
            }
        }
    // Calculate Wigner
    auto &all_wigner = s.parameters.input_correlation["Wigner"];
    for ( auto &wigner_s : all_wigner )
        for ( size_t i = 0; i < wigner_s.string_v["Modes"].size(); i++ ) {
            calculate_wigner( s, wigner_s.string_v["Modes"][i], wigner_s.property_set["X"][i], wigner_s.property_set["Y"][i], wigner_s.property_set["Res"][i], wigner_s.property_set["Skip"][i] );
        }

    // Output Spectra and Rest in seperate Files
    for ( auto &[mode, data] : to_output["Spectrum"] ) {
        Log::L2( "[PhotonStatistics] Saving Emission Spectrum to spectrum_" + mode + ".txt...\n" );
        auto &f_spectrum = FileOutput::add_file( "spectrum_" + mode );
        f_spectrum << std::format( "Omega\t{}\n", mode );
        for ( size_t i = 0; i < to_output["Spectrum"][mode].size(); i++ ) {
            f_spectrum << std::format( "{:.8e}\t{:.8e}\n", std::real( to_output["Spectrum_frequency"][mode][i] ), std::real( to_output["Spectrum"][mode][i] ) );
        }
    }
    for ( auto &[mode, data] : to_output["Indist"] ) {
        if ( mode.compare( "Time" ) == 0 ) continue;
        Log::L2( "[PhotonStatistics] Saving Indistinguishability and Visibility to indist_" + mode + ".txt...\n" );
        auto &f_indist = FileOutput::add_file( "indist_" + mode );
        f_indist << std::format( "Time\tIndist_{}\tVisibility_{}\n", mode, mode );
        for ( size_t i = 0; i < to_output["Indist"][mode].size(); i++ ) {
            f_indist << std::format( "{:.8e}\t{:.8e}\t{:.8e}\n", std::real( to_output["Indist"]["Time"][i] ), std::real( to_output["Indist"][mode][i] ), std::real( to_output["Visibility"][mode][i] ) );
        }
    }
    for ( auto &[mode, data] : to_output["Conc"] ) {
        if ( mode.compare( "Time" ) == 0 ) continue;
        Log::L2( "[PhotonStatistics] Saving Concurrence to conc_" + mode + ".txt...\n" );
        auto &f_conc = FileOutput::add_file( "conc_" + mode );
        auto &time = to_output["TwoPMat"]["Time"];
        f_conc << std::format( "Time\t{0}\t{0}_simple\t{0}_analytical\t{0}_fidelity\t{0}(g2(0))\t{0}_simple(g2(0))\t{0}_fidelity(g2(0))\n", mode );
        for ( size_t i = 0; i < to_output["Conc"][mode].size(); i++ ) {
            f_conc << std::format( "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", std::real( time[i] ), std::real( to_output["Conc"][mode][i] ), std::real( to_output["Conc_simple"][mode][i] ),
                                   std::real( to_output["Conc_analytical"][mode][i] ), std::real( to_output["Conc_fidelity"][mode][i] ), std::real( to_output["Conc_g2zero"][mode][i] ),
                                   std::real( to_output["Conc_g2zero_simple"][mode][i] ), std::real( to_output["Conc_g2zero_fidelity"][mode][i] ) );
        }
        if ( s.parameters.output_dict.contains( "conc" ) ) {
            auto &f_ev = FileOutput::add_file( "conc_eigenvalues_" + mode );
            f_ev << std::format( "Time\tRe({0})_0\tRe({0})_1\tRe({0})_2\tRe({0})_3\tIm({0})_0\tIm({0})_1\tIm({0})_2\tIm({0})_3\n", mode );
            for ( size_t i = 0; i < to_output["Conc"][mode].size(); i++ ) {
                f_ev << std::format( "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", std::real( time[i] ), std::real( to_output_m["ConcEV"][mode + "_EV"][i]( 0 ) ),
                                     std::real( to_output_m["ConcEV"][mode + "_EV"][i]( 1 ) ), std::real( to_output_m["ConcEV"][mode + "_EV"][i]( 2 ) ), std::real( to_output_m["ConcEV"][mode + "_EV"][i]( 3 ) ),
                                     std::imag( to_output_m["ConcEV"][mode + "_EV"][i]( 0 ) ), std::imag( to_output_m["ConcEV"][mode + "_EV"][i]( 1 ) ), std::imag( to_output_m["ConcEV"][mode + "_EV"][i]( 2 ) ),
                                     std::imag( to_output_m["ConcEV"][mode + "_EV"][i]( 3 ) ) );
            }
        }
    }
    for ( auto &[mode, data] : to_output["Raman"] ) {
        if ( mode.compare( "Time" ) == 0 ) continue;
        Log::L2( "[PhotonStatistics] Saving Raman Population to raman_" + mode + ".txt...\n" );
        auto &f_raman = FileOutput::add_file( "raman_" + mode );
        f_raman << std::format( "Time\t{0}\tEM({0})\n", mode );
        for ( size_t i = 0; i < to_output["Raman"][mode].size(); i++ ) {
            const auto &current_state = savedStates.at( i );
            f_raman << std::format( "{:.8e}\t{:.8e}\t{:.8e}\n", current_state.t, std::real( to_output["Raman"][mode][i] ), std::real( to_output["RamanEmProb"][mode][i] ) );
        }
    }
    if ( s.parameters.output_dict.contains( "tpm" ) )
        for ( auto &[mode, data] : to_output_m["TwoPMat"] ) {
            if ( mode.compare( "Time" ) == 0 ) continue;
            Log::L2( "[PhotonStatistics] Saving Two-Photon Matrix to twopmat_" + mode + ".txt...\n" );
            auto &mat = to_output_m["TwoPMat"][mode][0];
            auto dim = mat.rows();
            auto &f_twophot = FileOutput::add_file( "twopmat_" + mode );
            f_twophot << "Time\t";
            std::vector<std::string> modes = { "11", "12", "21", "22" };
            if (dim == 8)
                modes = {"111", "112", "121", "122", "211", "212", "221", "222"};
            for ( size_t k = 0; k < dim; k++ ) {
                for ( size_t l = 0; l < dim; l++ ) {
                    f_twophot << std::format( "Re({}{})\t", modes[l], modes[k] );
                }
            }
            for ( size_t k = 0; k < dim; k++ ) {
                for ( size_t l = 0; l < dim; l++ ) {
                    f_twophot << std::format( "Im({}{})\t", modes[l], modes[k] );
                }
            }
            f_twophot << "\n";
            for ( size_t i = 0; i < to_output_m["TwoPMat"][mode].size(); i++ ) {
                f_twophot << std::format( "{:.8e}\t", std::real( to_output["TwoPMat"]["Time"][i] ) );
                auto &mat = to_output_m["TwoPMat"][mode][i];
                for ( size_t k = 0; k < dim; k++ ) {
                    for ( size_t l = 0; l < dim; l++ ) {
                        f_twophot << std::format( "{:.8e}\t", std::real( mat( k, l ) ) );
                    }
                }
                for ( size_t k = 0; k < dim; k++ ) {
                    for ( size_t l = 0; l < dim; l++ ) {
                        f_twophot << std::format( "{:.8e}", std::imag( mat( k, l ) ) );
                        if ( k == dim-1 && l == dim-1 )
                            f_twophot << "\n";
                        else
                            f_twophot << "\t";
                    }
                }
            }
        }
    for ( auto &[mode, data] : to_output_m["Wigner"] ) {
        if ( mode.compare( "Time" ) == 0 ) continue;
        Log::L2( "[PhotonStatistics] Saving Wigner function to wigner_" + mode + ".txt...\n" );
        auto &f_wigner = FileOutput::add_file( "wigner_" + mode );
        if ( mode.starts_with( "rho_" ) ) {
            // Gather base:
            f_wigner << "Time\t";
            std::string smode = mode.substr( 4 );
            int base = s.operatorMatrices.el_states.contains( smode ) ? s.operatorMatrices.el_states[smode].base : s.operatorMatrices.ph_states[smode].base;
            if ( base == 0 ) {
                for ( auto &[name, dat] : s.parameters.input_electronic )
                    for ( auto &[name2, dat2] : s.parameters.input_electronic ) f_wigner << std::format( "Re(|{}><{}|)\t", name, name2 );
                for ( auto &[name, dat] : s.parameters.input_electronic )
                    for ( auto &[name2, dat2] : s.parameters.input_electronic ) f_wigner << std::format( "Im(|{}><{}|)\t", name, name2 );
                f_wigner << "\n";
            } else {
                for ( size_t i = 0; i < data[0].rows(); i++ )
                    for ( size_t j = 0; j < data[0].rows(); j++ ) f_wigner << std::format( "Re(|{}_{}><{}_{}|)\t", smode, i, smode, j );
                for ( size_t i = 0; i < data[0].rows(); i++ )
                    for ( size_t j = 0; j < data[0].rows(); j++ ) {
                        f_wigner << std::format( "Im(|{}_{}><{}_{}|)", smode, i, smode, j );
                        if ( i == data[0].rows() - 1 && j == data[0].rows() - 1 )
                            f_wigner << "\n";
                        else
                            f_wigner << "\t";
                    }
            }
        } else {
            f_wigner << std::format( "Time\t{}\n", mode );
        }
        for ( size_t i = 0; i < data.size(); i++ ) {
            f_wigner << std::format( "{:.8e}\t", std::real( to_output["Wigner"]["Time"][i] ) );
            auto &currentwigner = data[i];
            for ( size_t k = 0; k < currentwigner.rows(); k++ ) {
                for ( size_t l = 0; l < currentwigner.cols(); l++ ) {
                    f_wigner << std::format( "{:.8e}\t", std::real( currentwigner( k, l ) ) );
                }
            }
            for ( size_t k = 0; k < currentwigner.rows(); k++ ) {
                for ( size_t l = 0; l < currentwigner.cols(); l++ ) {
                    f_wigner << std::format( "{:.8e}", std::imag( currentwigner( k, l ) ) );
                    if ( k == currentwigner.rows() - 1 && l == currentwigner.cols() - 1 )
                        f_wigner << "\n";
                    else
                        f_wigner << "\t";
                }
            }
        }
    }
    return true;
}

// electronic_transition1, electronic_transition2, optical_transition
bool QDACC::Numerics::ODESolver::calculate_raman_population( System &s, const std::string &source_transitions, const std::string &raman_transition, const std::string &optical_transition, const std::string &pulse_mode ) {
    // auto [m_electronic_transition1, m_electronic_transition2] = get_operators_matrices( s, source_transitions, );
    // if ( s.operatorMatrices.el_transitions.count( raman_transition ) == 0 ) {
    //     auto ket = s.operatorMatrices.el_states[raman_transition.substr( 0, 1 )].ket;
    //     auto bra = s.operatorMatrices.el_states[raman_transition.substr( 1, 1 )].bra;
    //     auto current = s.operatorMatrices.base_selfhilbert;
    //     current.front() = ket * bra;
    //     Sparse transition_hilbert = QDACC::Matrix::tensor( current ).sparseView();
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
    // int pulse_index = s.parameters.input_pulse[pulse_mode].property["PulseIndex"];
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
    ////        R = s.dgl_expectationvalue<Sparse>( savedStates.at( i ).mat, op, tdd ) * std::conj( s.pulse.at( pulse_index ).get( tdd ) );
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
    //            Scalar R = s.dgl_expectationvalue<Sparse>( savedStates.at( j ).mat, op, tdd ) * s.pulse.at( pulse_index ).get( tdd );
    //            integral += 2.0i * s.parameters.p_omega_coupling * std::exp( -s.parameters.p_omega_cavity_loss * ( t - td ) ) * B * R * dtau * dt;
    //            // std::cout << s.dgl_expectationvalue<Sparse>( savedStates.at( j ).mat, op, tdd ) << ", " << std::conj( s.pulse.at( pulse_index ).get( tdd ) ) << "," << dt << "," << dtau << "," << B << "," << R << "," << s.parameters.p_omega_coupling << "," << s.parameters.p_omega_cavity_loss << "," << std::exp( -s.parameters.p_omega_cavity_loss * ( T - td ) ) << " - " << integral << std::endl;
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