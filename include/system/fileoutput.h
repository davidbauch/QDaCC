#pragma once
#include "global.h"
#include "system/operatormatrices.h"
#include "system/parameters.h"

class FileOutput {
   public:
    FILE *fp_densitymatrix;
    FILE *fp_electronic;
    FILE *fp_photonic;
    FILE *fp_numerical;
    FILE *fp_eigenvalues;
    FileOutput(){};
    FileOutput( Parameters &p, OperatorMatrices &op );
    void close( Parameters &p );
};