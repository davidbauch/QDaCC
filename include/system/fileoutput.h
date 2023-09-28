#pragma once
#include "global.h"
#include "system/operatormatrices.h"
#include "system/parameters.h"

namespace QDACC {

class FileOutput {
   private:
    std::unordered_map<std::string, std::fstream> files;
    void Iinit( Parameters &p, OperatorMatrices &op );
    std::ifstream Iload_file( const std::string &name, const std::string &file_ending );
    std::fstream &Iadd_file( const std::string &name, const std::string &file_ending, const std::ios_base::openmode mode = std::ios::out );
    std::fstream &Iget_file( const std::string &name );
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
    static std::ifstream load_file( const std::string &name, const std::string &file_ending );
    static std::fstream &add_file( const std::string &name, const std::string &file_ending = "txt", const std::ios_base::openmode mode = std::ios::out );
    static std::fstream &get_file( const std::string &name );
    static bool close_file( const std::string &name );
    static bool close_all();
    static void init( Parameters &p, OperatorMatrices &op );
};

} // namespace QDACC