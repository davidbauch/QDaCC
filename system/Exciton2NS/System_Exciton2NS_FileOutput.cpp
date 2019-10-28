class FileOutput {
   public:
    FILE *fp_densitymatrix;
    FILE *fp_atomicinversion;
    FILE *fp_photonpopulation;
    FileOutput(){};
    FileOutput( const std::vector<std::string> filenames, const Parameters &p, const OperatorMatrices &op );
    void close();
};

FileOutput::FileOutput( const std::vector<std::string> filenames, const Parameters &p, const OperatorMatrices &op ) {
    logs.level2( "Creating FileOutputs... " );
    long unsigned int config_maximumNumberOfFiles = 3;
    if ( filenames.size() != config_maximumNumberOfFiles ) {
        logs.level2( "WARNING: Number of input elements mismatches maximum files... " );
    }
    fp_densitymatrix = std::fopen( ( p.subfolder + filenames.at( 0 ) ).c_str(), "w" );
    if ( !fp_densitymatrix )
        logs.level2( "\nCould not open file for densitymatrix!\n" );
    else {
        fmt::print( fp_densitymatrix, "States: " );
        if ( p.output_full_dm ) {
            for ( int i = 0; i < p.maxStates; i++ )
                for ( int j = 0; j < p.maxStates; j++ ) {
                    fmt::print( fp_densitymatrix, "|{}><{}|\t", op.base.at( i ), op.base.at( j ) );
                }
        } else
            for ( int i = 0; i < p.maxStates; i++ )
                fmt::print( fp_densitymatrix, "|{}><{}|\t", op.base.at( i ), op.base.at( i ) );
        fmt::print( fp_densitymatrix, "\n" );
    }
    fp_atomicinversion = std::fopen( ( p.subfolder + filenames.at( 1 ) ).c_str(), "w" );
    if ( !fp_atomicinversion )
        logs.level2( "\nCould not open file for atomic inversion!\n" );
    fp_photonpopulation = std::fopen( ( p.subfolder + filenames.at( 2 ) ).c_str(), "w" );
    if ( !fp_photonpopulation )
        logs.level2( "\nCould not open file for photonpopulation!\n" );
    logs.level2( "done!\n" );
}

void FileOutput::close() {
    logs.level2( "Closing file outputs..." );
    std::fclose( fp_densitymatrix );
    std::fclose( fp_atomicinversion );
    std::fclose( fp_photonpopulation );
    logs.level2( "Done!\n" );
}