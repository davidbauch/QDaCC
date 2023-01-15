#include "system/system.h"

using namespace QDLC;

double System::dgl_phonons_lindblad_coefficients( const double energy, const double coupling, const Scalar current_pulse, const char mode, const double scaling, const double sign ) {
    Log::L3( "[System-PME]     Calculating Lindblad Rates for integral border = {}, deltaE = {}, g = {}, Omega = {}, mode = {}, sign = {}, coupling scaling = {}, integral stepsize = {}\n", parameters.p_phonon_tcutoff, energy, coupling, current_pulse, mode, sign, scaling, parameters.numerics_subiterator_stepsize );
    double ret = 0;
    // Light/Pulse Mode
    if ( mode == 'L' ) {
        const double bpulsesquared = std::pow( std::abs( parameters.p_phonon_b * scaling * current_pulse ), 2.0 );
        const double nu = std::sqrt( bpulsesquared + energy * energy );

        Log::L3( "[System-PME]         Nu = {}, Omega*<B>^2 = {}\n", nu, bpulsesquared );

        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += parameters.numerics_subiterator_stepsize ) {
            if ( not phi_vector.contains( tau ) )
                phi_vector[tau] = dgl_phonons_phi( tau );
            const auto phi = phi_vector[tau];
            const Scalar f = ( energy * energy * std::cos( nu * tau ) + bpulsesquared ) / std::pow( nu, 2.0 );
            ret += std::real( ( std::cosh( phi ) - 1.0 ) * f + std::sinh( phi ) * std::cos( nu * tau ) ) - sign * std::imag( ( std::exp( phi ) - 1.0 ) * energy * std::sin( nu * tau ) / nu );
            i++;
        }
        ret *= 2.0 * bpulsesquared * parameters.numerics_subiterator_stepsize;
    }
    // Cavity Mode
    else if ( mode == 'C' ) {
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += parameters.numerics_subiterator_stepsize ) {
            if ( not phi_vector.contains( tau ) )
                phi_vector[tau] = dgl_phonons_phi( tau );
            const auto phi = phi_vector[tau];
            ret += std::real( std::exp( 1.0i * sign * energy * tau ) * ( std::exp( phi ) - 1.0 ) );
            i++;
        }
        ret *= std::pow( parameters.p_phonon_b * scaling * coupling, 2.0 ) * parameters.numerics_subiterator_stepsize;
    }
    // Check if Rate is NaN. Only possibility for NaN is when the rates should actually result in zero.
    if ( std::isnan( ret ) )
        ret = 0.0;
    Log::L3( "[System-PME]         Return value: {}\n", ret );
    return ret;
}

Sparse System::dgl_phonons_lindblad_contribution( const double t, const Sparse &rho ) {
    Sparse ret( rho.rows(), rho.cols() );
    Log::L3( "[System-PME] Calculating Lindblad-Type Phonon Rates for t = {}\n", t );
    const double chirpcorrection = chirp.empty() ? 0.0 : std::real( chirp.back().get( t ) );
    for ( auto &[name, mat] : parameters.input_pulse ) {
        int p = 0;
        for ( int m = 0; m < mat.string_v["CoupledTo"].size(); m++ ) {
            const auto &mode = mat.string_v["CoupledTo"][m];
            const auto &mode_transposed = operatorMatrices.el_transitions[mode].name_transposed;
            const auto &transition = operatorMatrices.el_transitions[mode].hilbert;
            const auto &transition_transposed = operatorMatrices.el_transitions[mode_transposed].hilbert;
            // std::cout << p << " m1 = " << mode << ", m2 = " << mode_transposed << std::endl;
            auto delta_E = mat.property_set["Frequency"][m] - operatorMatrices.el_transitions[mode].energy + chirpcorrection;
            auto r1 = dgl_phonons_lindblad_coefficients( delta_E, 0.0, pulse[p].get( t ), 'L', 1.0, -1 );
            auto r2 = dgl_phonons_lindblad_coefficients( delta_E, 0.0, pulse[p].get( t ), 'L', 1.0, 1 );
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
            const double phonon_scaling = parameters.input_electronic[upper_level].property["PhononCoupling"];
            const auto &mode_transposed = operatorMatrices.el_transitions[mode].name_transposed;
            const auto &transition = operatorMatrices.el_transitions[mode].hilbert;
            const auto &transition_transposed = operatorMatrices.el_transitions[mode_transposed].hilbert;
            const auto &optical_transition = operatorMatrices.ph_transitions[name + "b"].hilbert;
            const auto &optical_transition_transposed = operatorMatrices.ph_transitions[name + "bd"].hilbert;
            auto delta_E = mat.property["Energy"] - operatorMatrices.el_transitions[mode].energy + chirpcorrection;
            auto r1 = dgl_phonons_lindblad_coefficients( delta_E, mat.property_set["CouplingScaling"][c] * parameters.p_omega_coupling, 0.0, 'C', phonon_scaling, -1 );
            auto r2 = dgl_phonons_lindblad_coefficients( delta_E, mat.property_set["CouplingScaling"][c] * parameters.p_omega_coupling, 0.0, 'C', phonon_scaling, 1 );
            Log::L3( "[System-PME]     Cavity induced Phonon Transition rates: {}-{}bd ({}), {}-{}b ({})", mode, name, r1, mode_transposed, name, r2 );
            ret += r1 * dgl_lindblad( rho, transition * optical_transition_transposed, transition_transposed * optical_transition );
            ret += r2 * dgl_lindblad( rho, transition_transposed * optical_transition, transition * optical_transition_transposed );
            c++;
        }
    }
    return ret;
}