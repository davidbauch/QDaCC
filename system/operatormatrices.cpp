#include "operatormatrices.h"

void OperatorMatrices_Parent::init(const Parameters &p) {
    Timer &timer_operatormatrices = createTimer("Operator Matrices");
    timer_operatormatrices.start();
    logs.level2("Generating operator matrices... ");
    if (!generateOperators(p)) {
        logs.level2("Generating operator matrices failed! Exitting program...\n");
        logs.close();
        exit(EXIT_FAILURE);
    }
    timer_operatormatrices.end();
    logs.level2("successful. Elapsed time is {}ms\n",timer_operatormatrices.getWallTime(TIMER_MILLISECONDS));
}

MatrixXcd OperatorMatrices_Parent::create_photonic_operator(const int type, const int maxPhotons) {
    MatrixXcd ret = MatrixXcd::Zero(maxPhotons+1,maxPhotons+1);
    for (int i = 0; i < maxPhotons; i++) {
        if (type == OPERATOR_PHOTONIC_CREATE)
            ret(i+1,i) = sqrt(i+1);
        if (type == OPERATOR_PHOTONIC_ANNIHILATE)
            ret(i,i+1) = sqrt(i+1);
    }
    return ret;
}