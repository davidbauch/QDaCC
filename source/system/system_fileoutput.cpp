#include "system/fileoutput.h"

std::ofstream &FileOutput::add_file( const std::string &name, const std::string &file_ending ) {
    return Get().Iadd_file( name, file_ending );
}
std::ofstream &FileOutput::get_file( const std::string &name ) {
    return Get().Iget_file( name );
}
bool FileOutput::close_file( const std::string &name ) {
    return Get().Iclose_file( name );
}
bool FileOutput::close_all() {
    return Get().Iclose_all();
}
void FileOutput::init( Parameters &p, OperatorMatrices &op ) {
    return Get().Iinit( p, op );
}

std::ofstream &FileOutput::Iget_file( const std::string &name ) {
    return files[name];
}

std::ofstream &FileOutput::Iadd_file( const std::string &name, const std::string &file_ending ) {
    if ( files.contains( name ) )
        return get_file( name );
    files[name].open( path + name + "." + file_ending );
    if ( get_file( name ).is_open() )
        Log::L2( "[System-Fileoutput] Added file '{}'\n", name );
    else
        Log::Error( "[System-Fileoutput] Could not add file '{}'!\n", name );
    return get_file( name );
}

bool FileOutput::Iclose_file( const std::string &name ) {
    if ( not files.contains( name ) )
        return false;
    if ( not files[name].is_open() )
        return false;
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
        fp_densitymatrix << "t\t";
        if ( p.input_conf["DMconfig"].string["output_mode"] == "full" ) {
            for ( int i = 0; i < op.base.size(); i++ )
                for ( int j = 0; j < op.base.size(); j++ ) {
                    fp_densitymatrix << fmt::format( "Re(|{}><{}|)\t", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ), op.base.at( j ).substr( 1, op.base.at( j ).size() - 2 ) );
                }
            for ( int i = 0; i < op.base.size(); i++ )
                for ( int j = 0; j < op.base.size(); j++ ) {
                    fp_densitymatrix << fmt::format( "Im(|{}><{}|)\t", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ), op.base.at( j ).substr( 1, op.base.at( j ).size() - 2 ) );
                }
        } else {
            for ( int i = 0; i < op.base.size(); i++ )
                fp_densitymatrix << fmt::format( "|{0}><{0}|\t", op.base.at( i ).substr( 1, op.base.at( i ).size() - 2 ) );
        }
        fp_densitymatrix << "\n";
    }

    // Electronic
    if ( not p.input_electronic.empty() ) {
        auto &fp_electronic = add_file( "electronic" );
        fp_electronic << "t\t";
        for ( auto &[name, rem] : p.input_electronic )
            fp_electronic << fmt::format( "|{}><{}|\t", name, name );
        if ( p.p_omega_decay > 0.0 )
            for ( auto &[name, rem] : p.input_electronic )
                if ( rem.numerical["DecayScaling"] != 0.0 )
                    fp_electronic << fmt::format( "EM(|{}><{}|)\t", name, name );
        fp_electronic << "\n";
    }

    // Photonic
    if ( not p.input_photonic.empty() ) {
        auto &fp_photonic = add_file( "photonic" );
        fp_photonic << "t\t";
        for ( auto &[name, rem] : p.input_photonic )
            fp_photonic << fmt::format( "|{}><{}|\t", name, name );
        if ( p.p_omega_cavity_loss > 0.0 )
            for ( auto &[name, rem] : p.input_photonic )
                fp_photonic << fmt::format( "EM(|{}><{}|)\t", name, name );
        fp_photonic << "\n";
    }
}