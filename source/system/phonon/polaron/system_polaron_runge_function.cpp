#include "system/system.h"

using namespace QDLC;

Sparse System::dgl_phonons_rungefunc( const Sparse &chi, const double t ) {
    // TODO Maybe? Cache explicit times?
    Scalar chirpcorrection = not chirp.empty() ? ( chirp.back().get( t ) + t * ( chirp.back().get( t ) - parameters.scaleVariable( chirp.back().derivative( t ), parameters.scale_value ) ) ) : 0;
    auto explicit_time = Sparse( chi.rows(), chi.cols() );

    // Electronic Transitions
    for ( const auto &[mode, param] : operatorMatrices.el_transitions ) {
        if ( param.direction == -1 )
            continue;
        explicit_time += ( param.energy + chirpcorrection ) * param.projector;
    }

    // Photonic Transitions
    for ( auto &[mode, param] : parameters.input_photonic ) {
        for ( const auto &transition : param.string_v["CoupledTo"] ) {
            auto transition_transposed = operatorMatrices.el_transitions[transition].name_transposed;
            explicit_time += 1.0i * ( operatorMatrices.el_transitions[transition_transposed].energy + chirpcorrection - operatorMatrices.ph_states[mode].energy ) * operatorMatrices.el_transitions[transition_transposed].projector * operatorMatrices.ph_transitions[mode + "b"].projector;
        }
    }

    // Pulse Driving
    int p = 0;
    int pp = 0;
    for ( auto &[mode, param] : parameters.input_pulse ) {
        for ( const auto &transition : param.string_v["CoupledTo"] ) {
            explicit_time += 1.0i * ( pulse[p].get( t ) + 1.0i * parameters.scaleVariable( pulse[p].derivative( t, parameters.t_step ), parameters.scale_value ) ) * operatorMatrices.polaron_pulse_factors_explicit_time[pp++];
        }
        p++;
    }
    // explicit_time = parameters.scaleVariable( explicit_time, 1.0 / parameters.scale_value ); // Deprecated

    Sparse hamilton = dgl_get_hamilton( t );
    return -1.0i * dgl_kommutator( hamilton, chi ) + explicit_time.cwiseProduct( QDLC::Matrix::sparse_projector( chi ) );
}