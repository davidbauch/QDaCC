#include "system/fileoutput.h"

FileOutput::FileOutput( Parameters &p, OperatorMatrices &op ) {
    LOG2( "[System-Fileoutput] Creating FileOutputs...\n" );
    if ( p.input_conf["DMconfig"].string["output_mode"] != "none" ) {
        fp_densitymatrix = std::fopen( ( p.working_directory + "densitymatrix.txt" ).c_str(), "w" );
        if ( !fp_densitymatrix ) {
            LOG2( "[System-Fileoutput] Could not open file for densitymatrix!\n" );
        } else {
            fmt::print( fp_densitymatrix, "t\t" );
            if ( p.input_conf["DMconfig"].string["output_mode"] == "full" ) {
                for ( int i = 0; i < op.base.size(); i++ )
                    for ( int j = 0; j < op.base.size(); j++ ) {
                        fmt::print( fp_densitymatrix, "Re(|{}><{}|)\t", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ), op.base.at( j ).substr( 1, op.base.at( j ).size() - 2 ) );
                    }
                for ( int i = 0; i < op.base.size(); i++ )
                    for ( int j = 0; j < op.base.size(); j++ ) {
                        fmt::print( fp_densitymatrix, "Im(|{}><{}|)\t", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ), op.base.at( j ).substr( 1, op.base.at( j ).size() - 2 ) );
                    }
            } else {
                for ( int i = 0; i < op.base.size(); i++ )
                    fmt::print( fp_densitymatrix, "|{0}><{0}|\t", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ) );
            }
            fmt::print( fp_densitymatrix, "\n" );
        }
    }
    fp_electronic = std::fopen( ( p.working_directory + "electronic.txt" ).c_str(), "w" );
    if ( !fp_electronic ) {
        LOG2( "[System-Fileoutput] Could not open file for atomic inversion!\n" );
    } else {
        fmt::print( fp_electronic, "t\t" ); //|G><G|\t|X_H><X_H|\t|X_V><X_V|\t|B><B|\n" );
        for ( auto &[name, rem] : p.input_electronic )
            fmt::print( fp_electronic, "|{}><{}|\t", name, name );
        if ( p.p_omega_decay > 0.0 )
            for ( auto &[name, rem] : p.input_electronic )
                if ( rem.numerical["DecayScaling"] != 0.0 )
                    fmt::print( fp_electronic, "EM(|{}><{}|)\t", name, name );
        fmt::print( fp_electronic, "\n" );
    }
    fp_photonic = std::fopen( ( p.working_directory + "photonic.txt" ).c_str(), "w" );
    if ( !fp_photonic ) {
        LOG2( "[System-Fileoutput] Could not open file for photonpopulation!\n" );
    } else {
        fmt::print( fp_photonic, "t\t" ); //|G><G|\t|X_H><X_H|\t|X_V><X_V|\t|B><B|\n" );
        for ( auto &[name, rem] : p.input_photonic )
            fmt::print( fp_photonic, "|{}><{}|\t", name, name );
        if ( p.p_omega_cavity_loss > 0.0 )
            for ( auto &[name, rem] : p.input_photonic )
                fmt::print( fp_photonic, "EM(|{}><{}|)\t", name, name );
        fmt::print( fp_photonic, "\n" );
    }
    if ( p.numerics_output_rkerror ) {
        fp_numerical = std::fopen( ( p.working_directory + "numerical.txt" ).c_str(), "w" );
        if ( !fp_numerical ) {
            LOG2( "[System-Fileoutput] Could not open file for numerical data!\n" );
        }
    }
    if ( p.output_eigenvalues ) {
        fp_eigenvalues = std::fopen( ( p.working_directory + "eigenvalues.txt" ).c_str(), "w" );
        if ( !fp_eigenvalues ) {
            LOG2( "[System-Fileoutput] Could not open file for eigenvalue data!\n" );
        }
    }
}

void FileOutput::close( Parameters &p ) {
    LOG2( "[System-Fileoutput] Closing file outputs...\n" );
    LOG2( "[System-Fileoutput] Closing Density Matrix Output\n" );
    if ( p.input_conf["DMconfig"].string["output_mode"] != "none" ) {
        std::fclose( fp_densitymatrix );
    }
    LOG2( "[System-Fileoutput] Closing Electronic States Output\n" );
    std::fclose( fp_electronic );
    LOG2( "[System-Fileoutput] Closing Photonic States Output\n" );
    std::fclose( fp_photonic );
    if ( p.numerics_output_rkerror ) {
        LOG2( "[System-Fileoutput] Closing Numerical Output\n" );
        std::fclose( fp_numerical );
    }
    if ( p.output_eigenvalues ) {
        LOG2( "[System-Fileoutput] Closing Eigenvalue Output\n" );
        std::fclose( fp_eigenvalues );
    }
}