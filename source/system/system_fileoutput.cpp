#include "system/fileoutput.h"

using namespace QDACC;

std::ifstream FileOutput::load_file( const std::string &name, const std::string &file_ending ) { return Get().Iload_file( name, file_ending ); }
std::fstream &FileOutput::add_file( const std::string &name, const std::string &file_ending, const std::ios_base::openmode mode ) { return Get().Iadd_file( name, file_ending ); }
std::fstream &FileOutput::get_file( const std::string &name ) { return Get().Iget_file( name ); }
bool FileOutput::close_file( const std::string &name ) { return Get().Iclose_file( name ); }
bool FileOutput::close_all() { return Get().Iclose_all(); }
void FileOutput::init( Parameters &p, OperatorMatrices &op ) { return Get().Iinit( p, op ); }

std::fstream &FileOutput::Iget_file( const std::string &name ) {
    if ( not files.contains( name ) ) {
#pragma omp critical
        files[name] = std::fstream();
    }
    return files[name];
}

std::ifstream FileOutput::Iload_file( const std::string &name, const std::string &file_ending ) {
    auto file = std::ifstream( path + name + "." + file_ending );
    if ( file.is_open() )
        Log::L2( "[System-Fileoutput] Loaded file '{}.{}'\n", name, file_ending );
    else
        Log::Error( "[System-Fileoutput] Could not load file '{}.{}'!\n", name, file_ending );
    return file;
}

std::fstream &FileOutput::Iadd_file( const std::string &name, const std::string &file_ending, const std::ios_base::openmode mode ) {
    if ( files.contains( name ) and get_file( name ).is_open() ) {
        return get_file( name );
    }
    get_file( name ).open( path + name + "." + file_ending, mode );
    if ( get_file( name ).is_open() )
        Log::L2( "[System-Fileoutput] Added file '{}.{}'\n", name, file_ending );
    else
        Log::Error( "[System-Fileoutput] Could not add file '{}.{}'!\n", name, file_ending );
    return get_file( name );
}

bool FileOutput::Iclose_file( const std::string &name ) {
    if ( not files.contains( name ) ) return false;
    if ( not files[name].is_open() ) return false;
    Log::L2( "[System-Fileoutput] Closing '{}.txt'\n", name );
    files[name].close();
    return true;
}

bool FileOutput::Iclose_all() {
    Log::L2( "[System-Fileoutput] Closing file outputs...\n" );
    bool success = true;
    for ( auto &[name, file] : files ) {
        success = success and close_file( name );
    }
    files.clear();
    return success;
}

void FileOutput::Iinit( Parameters &p, OperatorMatrices &op ) {
    path = p.working_directory;
    Log::L2( "[System-Fileoutput] Creating FileOutputs...\n" );
    // Density Matrix
    if ( p.input_conf["DMconfig"].string["output_mode"] != "none" ) {
        auto &fp_densitymatrix = add_file( "densitymatrix" );
        fp_densitymatrix << "t";
        if ( p.input_conf["DMconfig"].string["output_mode"] == "full" ) {
            for ( int i = 0; i < op.base.size(); i++ )
                for ( int j = 0; j < op.base.size(); j++ ) {
                    fp_densitymatrix << std::format( "\tRe(|{}><{}|)", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ), op.base.at( j ).substr( 1, op.base.at( j ).size() - 2 ) );
                }
            for ( int i = 0; i < op.base.size(); i++ )
                for ( int j = 0; j < op.base.size(); j++ ) {
                    fp_densitymatrix << std::format( "\tIm(|{}><{}|)", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ), op.base.at( j ).substr( 1, op.base.at( j ).size() - 2 ) );
                }
        } else {
            for ( int i = 0; i < op.base.size(); i++ ) fp_densitymatrix << std::format( "\t|{0}><{0}|", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ) );
        }
        fp_densitymatrix << "\n";
    }

    // Electronic
    if ( not p.input_electronic.empty() ) {
        auto &fp_electronic = add_file( "electronic" );
        fp_electronic << "t";
        for ( auto &[name, rem] : p.input_electronic ) fp_electronic << std::format( "\t|{}><{}|", name, name );
        if ( p.p_omega_decay > 0.0 )
            for ( auto &[name, rem] : p.input_electronic )
                if ( rem.property["DecayScaling"] != 0.0 ) fp_electronic << std::format( "\tEM(|{}><{}|)\t", name, name );
        fp_electronic << "\n";
    }

    // Photonic
    if ( not p.input_photonic.empty() ) {
        auto &fp_photonic = add_file( "photonic" );
        fp_photonic << "t";
        for ( auto &[name, rem] : p.input_photonic ) fp_photonic << std::format( "\t|{}><{}|", name, name );
        if ( p.p_omega_cavity_loss > 0.0 )
            for ( auto &[name, rem] : p.input_photonic ) fp_photonic << std::format( "\tEM(|{}><{}|)", name, name );
        fp_photonic << "\n";
    }
    // Photon Expv
    if ( p.output_dict.contains( "photons" ) ) {
        for ( auto &[mode, state] : op.ph_states ) {
            auto& fp_photons = add_file( "photons_" + mode );
            fp_photons << "t";
            for ( int i = 0; i < state.self_hilbert.rows(); i++ )
                for ( int j = 0; j < state.self_hilbert.cols(); j++ ) fp_photons << std::format( "\t|{}><{}|", i, j );
            if ( p.p_omega_cavity_loss > 0.0 )
                for ( int i = 1; i < state.self_hilbert.rows(); i++ ) fp_photons << std::format( "\tEM(|{}><{}|)", i, i ); // Start at 1 because 0 is the vacuum state and can never decay
            fp_photons << "\n";
        }
    }
    // Custom Expv
    if ( not p.numerics_custom_expectation_values.empty() ) {
        add_file( "custom_expectation_values" ) << "t";
        for ( auto i = 0; i < p.numerics_custom_expectation_values.size(); i++ ) get_file( "custom_expectation_values" ) << std::format( "\tCustom_{}", i );
        get_file( "custom_expectation_values" ) << "\n";
    }
}