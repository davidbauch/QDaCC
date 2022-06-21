#pragma once
#include "global.h"
#include "system/operatormatrices.h"
#include "system/parameters.h"

class FileOutput {
   private:
    std::unordered_map<std::string, std::ofstream> files;
    void Iinit( Parameters &p, OperatorMatrices &op );
    std::ofstream &Iadd_file( const std::string &name );
    std::ofstream &Iget_file( const std::string &name );
    bool Iclose_file( const std::string &name );
    bool Iclose_all();
    std::string path;

   public:
    FileOutput() = default;
    FileOutput( const FileOutput & ) = delete;
    static FileOutput &Get() {
        static FileOutput instance;
        return instance;
    }
    static std::ofstream &add_file( const std::string &name );
    static std::ofstream &get_file( const std::string &name );
    static bool close_file( const std::string &name );
    static bool close_all();
    static void init( Parameters &p, OperatorMatrices &op );
};