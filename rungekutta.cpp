#pragma once
#include "global.h"

class ODESolver{
    private:
        // RK parameters
        double a2 = 1./5.;
        double a3 = 3./10.;
        double a4 = 4./5.;
        double a5 = 8./9.;
        double a6 = 1.0;

        double b11 = 1./5.;
        double b21 = 3./40.;
        double b31 = 44./45.;
        double b41 = 19372./6561.;
        double b51 = 9017./3168.;
        double b61 = 35./384.;
        double b22 = 9./40.;
        double b32 = -56./15.;
        double b42 = -25360./2187.;
        double b52 = -355./33.;
        double b33 = 32./9.;
        double b43 = 64448./6561.;
        double b53 = 46732./5247.;
        double b63 = 500./1113.;
        double b44 = -212./729.;
        double b54 = 49./176.;
        double b64 = 125./192.;
        double b55 = -5103/18656.;
        double b65 = -2187./6784.;
        double b66 = 11./84.;
        
        // Vector of Hamiltonmatrices for timesteps. Save time and matrix in extra class so mapping is easy
        // Disable to save memory
        class savedMatrix {
            public:
                MatrixXcd matrix;
                double t;
                savedMatrix(MatrixXcd _matrix, const double _t) {
                    t = _t;
                    matrix = _matrix;
                }
        };
        bool saveMatrices;
        //Timer test = createTimer("placeholder",false);
        std::vector<savedMatrix> savedMatrices;
        //Timer &rk_Hamiltons = createTimer("RungeKutta-Hamiltons",false); // die sind nicht threadsafe lol.
        //Timer &rk_kMatrices = createTimer("RungeKutta-K-Matrices",false);
        //allTimers.erase(allTimers.back()-2);
    public:
        ODESolver(System &s, bool _saveMatrices = false) {
            logs.level2("Creating ODESolver Class... ");
            saveMatrices = _saveMatrices;
            if (saveMatrices) {
                int matricesToSave = (s.parameters.t_end-s.parameters.t_start)/s.parameters.t_step * ( s.getSolverOrder(DIR_T)+s.getSolverOrder(DIR_TAU) > 8 ? 5 : 3 );
                savedMatrices.reserve(matricesToSave);
                logs.level2("Saving matrices... ");
            }
            logs.level2("Done!\n");
        }

        MatrixXcd getHamilton(const System &s, const double t) {
            if (!saveMatrices) {
                return s.dgl_getHamilton(t);
            }
            int i = std::floor(t/s.parameters.t_step);
            int start = i;
            // FIXME: für savedMatrices erforderlich, das RK(T) = 5, ODER das RK(T) = RK(Tau), lazy method aber ansosnten
            // schreibt RK(T)=4 t werte in den vector, und RK(Tau)=5 appendet neue hinten, für viel frühere zeiten, quasi nicht zu mappen
            // dafür müssen für tau alle steps bekannt sien (kein 4->5, nur 5->4, 4->4, 5->5), d.h. fehlberg auch für 4 benutzen!
            while (t < savedMatrices.at(i).t && i+20 < savedMatrices.size())
                i+=20;
            bool found = false;
            while (i >= start-20 && !found) {
                if (std::abs(t - savedMatrices.at(i).t) > s.parameters.t_step/20.0) {
                    found = true;
                }
                i--;
            }
            if (!found) {
                savedMatrix saved = savedMatrix(s.dgl_getHamilton(t),t);
                savedMatrices.emplace_back(saved);
                return saved.matrix;
            }
            return savedMatrices.at(i).matrix;
        }
        
        MatrixXcd iterateRungeKutta4(const MatrixXcd &rho, const System &s, const double t){
            /* Verschiedene H's fuer k1-4 ausrechnen */
            //double wallt = omp_get_wtime();
            MatrixXcd H_calc_k1  = getHamilton(s,t);
            MatrixXcd H_calc_k23 = getHamilton(s,t+s.parameters.t_step*0.5);
            MatrixXcd H_calc_k4  = getHamilton(s,t+s.parameters.t_step);
            //rk_Hamiltons.add(0, (omp_get_wtime()-wallt)/numThreads);
            /* k1-4 ausrechnen */
            //wallt = omp_get_wtime();
            MatrixXcd rk1 = s.dgl_rungeFunction( rho,                                  H_calc_k1,  t            );
            MatrixXcd rk2 = s.dgl_rungeFunction( rho + s.parameters.t_step*0.5*rk1,    H_calc_k23, t+s.parameters.t_step*0.5 );
            MatrixXcd rk3 = s.dgl_rungeFunction( rho + s.parameters.t_step*0.5*rk2,    H_calc_k23, t+s.parameters.t_step*0.5 );
            MatrixXcd rk4 = s.dgl_rungeFunction( rho + s.parameters.t_step*rk3,        H_calc_k4,  t+s.parameters.t_step     );
            //rk_kMatrices.add(0, (omp_get_wtime()-wallt)/numThreads);
            /* Dichtematrix */
            return rho + s.parameters.t_step/6.0*( rk1 + 2.*rk2 + 2.*rk3 + rk4 );
        }

        MatrixXcd iterateRungeKutta5(const MatrixXcd &rho, const System &s, const double t) {
            // Verschiedene H's fuer k1-6 ausrechnen
            //double wallt = omp_get_wtime();
            MatrixXcd H_calc_k1  = getHamilton(s,t);
            MatrixXcd H_calc_k2  = getHamilton(s,t+a2*s.parameters.t_step);
            MatrixXcd H_calc_k3  = getHamilton(s,t+a3*s.parameters.t_step);
            MatrixXcd H_calc_k4  = getHamilton(s,t+a4*s.parameters.t_step);
            MatrixXcd H_calc_k5  = getHamilton(s,t+a5*s.parameters.t_step);
            MatrixXcd H_calc_k6  = getHamilton(s,t+a6*s.parameters.t_step);
            //rk_Hamiltons.add(0, (omp_get_wtime()-wallt)/numThreads);
            // k1-4 ausrechnen
            //wallt = omp_get_wtime();
            MatrixXcd k1 = s.dgl_rungeFunction( rho,                                                                    H_calc_k1, t            );
            MatrixXcd k2 = s.dgl_rungeFunction( rho + s.parameters.t_step*b11*k1,                                       H_calc_k2, t+a2*s.parameters.t_step );
            MatrixXcd k3 = s.dgl_rungeFunction( rho + s.parameters.t_step*(b21*k1 + b22*k2),                            H_calc_k3, t+a3*s.parameters.t_step );
            MatrixXcd k4 = s.dgl_rungeFunction( rho + s.parameters.t_step*(b31*k1 + b32*k2 + b33*k3),                   H_calc_k4, t+a4*s.parameters.t_step );
            MatrixXcd k5 = s.dgl_rungeFunction( rho + s.parameters.t_step*(b41*k1 + b42*k2 + b43*k3 + b44*k4),          H_calc_k5, t+a5*s.parameters.t_step );
            MatrixXcd k6 = s.dgl_rungeFunction( rho + s.parameters.t_step*(b51*k1 + b52*k2 + b53*k3 + b54*k4 + b55*k5), H_calc_k6, t+a6*s.parameters.t_step );
            //rk_kMatrices.add(0, (omp_get_wtime()-wallt)/numThreads);
            // Dichtematrix
            return rho + s.parameters.t_step*( b61*k1 + b63*k3 + b64*k4 + b65*k5 + b66*k6 );
        }

        MatrixXcd iterate(const MatrixXcd &rho, const System &s, const double t, const int dir = DIR_T) {
            int order = s.getSolverOrder(dir);
            if (order == 4) {
                return iterateRungeKutta4(rho,s,t);
            } //else if (order == 5) {
            return iterateRungeKutta5(rho,s,t);
            //}
        }
};

/*
maddlab
double a2 = 1./5.;
double a3 = 3./10.;
double a4 = 4./5.;
double a5 = 8./9.;
double a6 = 1.0;

double b11 = 1./5.;
double b21 = 3./40.;
double b31 = 44./45.;
double b41 = 19372./6561.;
double b51 = 9017./3168.;
double b61 = 35./384.;
double b22 = 9./40.;
double b32 = -56./15.;
double b42 = -25360./2187.;
double b52 = -355./33.;
double b33 = 32./9.;
double b43 = 64448./6561.;
double b53 = 46732./5247.;
double b63 = 500./1113.;
double b44 = -212./729.;
double b54 = 49./176.;
double b64 = 125./192.;
double b55 = -5103/18656.;
double b65 = -2187./6784.;
double b66 = 11./84.;

a2=cast(1/5,dataType);
a3=cast(3/10,dataType);
a4=cast(4/5,dataType);
a5=cast(8/9,dataType);

b11=cast(1/5,dataType); 
b21=cast(3/40,dataType); 
b31=cast(44/45,dataType);
b41=cast(19372/6561,dataType);
b51=cast(9017/3168,dataType);
b61=cast(35/384,dataType);
b22=cast(9/40,dataType);
b32=cast(-56/15,dataType);
b42=cast(-25360/2187,dataType);
b52=cast(-355/33,dataType);
b33=cast(32/9,dataType);
b43=cast(64448/6561,dataType);
b53=cast(46732/5247,dataType);
b63=cast(500/1113,dataType);
b44=cast(-212/729,dataType);
b54=cast(49/176,dataType);
b64=cast(125/192,dataType);
b55=cast(-5103/18656,dataType);
b65=cast(-2187/6784,dataType);
b66=cast(11/84,dataType);
*/