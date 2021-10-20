#include "system/fileoutput.h"

FileOutput::FileOutput( Parameters &p, OperatorMatrices &op ) {
    Log::L2( "Creating FileOutputs... " );
    output_no_dm = p.output_no_dm;
    if ( !output_no_dm ) {
        fp_densitymatrix = std::fopen( ( p.subfolder + "densitymatrix.txt" ).c_str(), "w" );
        if ( !fp_densitymatrix )
            Log::L2( "\nCould not open file for densitymatrix!\n" );
        else {
            fmt::print( fp_densitymatrix, "t\t" );
            if ( p.output_full_dm ) {
                for ( int i = 0; i < op.base.size(); i++ )
                    for ( int j = 0; j < op.base.size(); j++ ) {
                        fmt::print( fp_densitymatrix, "|{}><{}|\t", op.base.at( i ).substr(1,op.base.at( i ).size()-2 ), op.base.at( j ).substr(1,op.base.at( j ).size()-2 ) );
                    }
            } else
                for ( int i = 0; i < op.base.size(); i++ )
                    fmt::print( fp_densitymatrix, "|{0}><{0}|\t", op.base.at( i ).substr(1,op.base.at( i ).size()-2 ) );
            fmt::print( fp_densitymatrix, "\n" );
        }
    }
    fp_atomicinversion = std::fopen( ( p.subfolder + "electronic.txt" ).c_str(), "w" );
    if ( !fp_atomicinversion )
        Log::L2( "\nCould not open file for atomic inversion!\n" );
    else {
        fmt::print( fp_atomicinversion, "t\t" ); //|G><G|\t|X_H><X_H|\t|X_V><X_V|\t|B><B|\n" );
        for ( auto &[name, rem] : p.input_electronic )
            fmt::print( fp_atomicinversion, "|{}><{}|\t", name, name );
        if ( p.p_omega_decay > 0.0 )
            for ( auto &[name, rem] : p.input_electronic )
                if ( rem.numerical["DecayScaling"] != 0.0 )
                    fmt::print( fp_atomicinversion, "EM(|{}><{}|)\t", name, name );
        fmt::print( fp_atomicinversion, "\n" );
    }
    fp_photonpopulation = std::fopen( ( p.subfolder + "photonic.txt" ).c_str(), "w" );
    if ( !fp_photonpopulation )
        Log::L2( "\nCould not open file for photonpopulation!\n" );
    else {
        fmt::print( fp_photonpopulation, "t\t" ); //|G><G|\t|X_H><X_H|\t|X_V><X_V|\t|B><B|\n" );
        for ( auto &[name, rem] : p.input_photonic )
            fmt::print( fp_photonpopulation, "|{}><{}|\t", name, name );
        if ( p.p_omega_cavity_loss > 0.0 )
            for ( auto &[name, rem] : p.input_photonic )
                fmt::print( fp_photonpopulation, "EM(|{}><{}|)\t", name, name );
        fmt::print( fp_photonpopulation, "\n" );
        //fmt::print( fp_photonpopulation, "t\tHorizontal\tVertical\tEmission-Probability-H\tEmission-Probability-V{}{}\n", ( p.numerics_output_raman_population ? "\tRaman-Population-H\tRaman-Poppulation-V\tRaman-Emission-Probability-H\tRaman-Emission-Probability-V" : "" ), ( p.numerics_output_electronic_emission ? "\tElectronic-H\tElectronic-V" : "" ) );
    }
    Log::L2( "done!\n" );
}

void FileOutput::close() {
    Log::L2( "Closing file outputs..." );
    if ( !output_no_dm ) {
        std::fclose( fp_densitymatrix );
    }
    std::fclose( fp_atomicinversion );
    std::fclose( fp_photonpopulation );
    Log::L2( "Done!\n" );
}