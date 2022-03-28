#pragma once
#include "global.h"
#include "system/operatormatrices.h"
#include "system/parameters.h"

class FileOutput {
   public:
    FILE *fp_densitymatrix;
    bool output_no_dm = false;
    FILE *fp_atomicinversion;
    FILE *fp_photonpopulation;
    FILE *fp_numerical;
    FileOutput(){};
    FileOutput( Parameters &p, OperatorMatrices &op );
    void close();
};