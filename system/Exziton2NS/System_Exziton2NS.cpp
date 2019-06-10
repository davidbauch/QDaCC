#include "System_Exziton2NS_Parameters.cpp"
#include "System_Exziton2NS_OperatorMatrices.cpp"
#include "System_Exziton2NS_FileOutput.cpp"

// FIXME: Header für System 
#include "../../chirp.h" 
#include "../../pulse.h"

class System {
    public:
        // Name
        std::string name;
        // Vector of input arguments
        std::vector<std::string> arguments;
        std::string terminateMessage;

        Chirp chirp;
        Pulse pulse;
        Parameters parameters;
        OperatorMatrices operatorMatrices;
        FileOutput fileoutput;

        System() {};
        System(int argc, char **argv) {
            terminateMessage = global_message_normaltermination; /* Successfull termination */
            //arguments(argv, argv+argc);
            for (int i = 0; i < argc; i++)
                arguments.push_back(std::string(argv[i]));
            name = "Exziton (2NS)";
            logs.level2("Creating System Class for '{}'\n",name);
            parameters = Parameters( arguments );
            parameters.log();
            operatorMatrices = OperatorMatrices( parameters );
            fileoutput = FileOutput( {"densitymatrix.txt", "atomicinversion.txt", "photonpopulation.txt"}, parameters );
        }
        
        void init() {
            Timer &timer = createTimer("System init");
            timer.start();
            logs.level2("Initializing system...\n");
            // Adding Timers
            
            // Rest
            chirp = Chirp(parameters); 
            /*if (parameters.chirp_file_type.compare("spline") == 0)
                chirp.generateFromSpline(parameters);
            else if (parameters.chirp_file_type.compare("chain") == 0)
                chirp.generateFromInputFile(parameters);
            else
                chirp.generateFromParam(parameters); */
            if (parameters.chirp_total != 0)
                chirp.fileOutput(parameters.subfolder + "chirp.txt", parameters);
            
            pulse = Pulse(parameters); 
            pulse.generateFromParam(parameters);
            if (parameters.pulse_amp != 0)
                pulse.fileOutput(parameters.subfolder + "pulse.txt", parameters);
            // Check time trafo
            MatrixXcd temp = dgl_timetrafo(MatrixXcd::Ones(parameters.maxStates,parameters.maxStates),1.0) - (-1i*operatorMatrices.H_0).exp()*MatrixXcd::Ones(parameters.maxStates,parameters.maxStates)*(1i*operatorMatrices.H_0).exp();
            if (std::abs(temp.sum()) >= 1E-15){
                logs("Unitary timetransformation does not match matrix exponential! MaxCoeff = {}\n",std::abs(temp.sum()));
            }
            timer.end();
            logs.level2("Done! Elapsed time is {}ms\n",timer.getWallTime()*1E3);
        }

        MatrixXcd getRho0() const {
            return operatorMatrices.rho;
        }

        MatrixXcd dgl_kommutator(const MatrixXcd &A, const MatrixXcd &B) const {
            return A*B - B*A;
        }

        MatrixXcd dgl_antikommutator(const MatrixXcd &A, const MatrixXcd &B) const {
            return A*B + B*A;
        }

        MatrixXcd dgl_lindblad(const MatrixXcd &rho, const MatrixXcd &op, const MatrixXcd &opd) const {
            return 2.0*op*rho*opd - opd*op*rho - rho*opd*op;
        }

        std::complex<double> dgl_expectationvalue(const MatrixXcd &rho, const MatrixXcd &op, double t) const {
            return (rho*dgl_timetrafo(op,t)).trace();
        }

        MatrixXcd dgl_calc_rhotau(const MatrixXcd &rho, const MatrixXcd &op, double t) const {
            return dgl_timetrafo(op,t)*rho; 
        }

        MatrixXcd dgl_rungeFunction(const MatrixXcd &rho, const MatrixXcd &H, const double t) const {
            MatrixXcd ret = -1i*dgl_kommutator(H,rho);
            // Photone losses
            if (parameters.k != 0.0) {
                ret += parameters.k/2.0*( 2*operatorMatrices.photon_annihilate*rho*operatorMatrices.photon_create - dgl_antikommutator(operatorMatrices.photon_n,rho) ); /* ret + k/2*(2*b*rho*b^+ - [b^+b,rho]) */
                //ret += parameters.k/2.0*dgl_lindblad(rho,operatorMatrices.photon_annihilate,operatorMatrices.photon_create);
            } 
            // Pure Dephasing
            if (parameters.gamma_pure != 0.0) {
                ret -= parameters.gamma_pure/2.0*( operatorMatrices.atom_ground*rho*operatorMatrices.atom_exited + operatorMatrices.atom_exited*rho*operatorMatrices.atom_ground ); /* -gamma_pure/2*(|g><g|rho|e><e| + |e><e|rho|g><g|) */
                //ret += parameters.gamma_pure/2.0*dgl_lindblad(rho,operatorMatrices.atom_exited,operatorMatrices.atom_ground);
            }
            if (parameters.gamma != 0.0) {
                ret += parameters.gamma*dgl_lindblad(rho,operatorMatrices.atom_sigmaminus,operatorMatrices.atom_sigmaplus);
            }
            return ret;
        }

        MatrixXcd dgl_timetrafo(const MatrixXcd &A, double t) const {
            MatrixXcd ret = A;
            if (parameters.doInteractionPicture == 1) {
                int i,j,pn,pm;
                for (int n = 0; n < A.rows(); n++) {
                    i = n%2;
                    pn = (int)(n/2);
                    for (int m = 0; m < A.cols(); m++) {
                        j = m%2;
                        pm = (int)(m/2);
                        ret(n,m) = A(n,m)*std::exp(1i*t*( (parameters.omegaE)/2.*( delta(i,1)-delta(i,0)-delta(j,1)+delta(j,0) ) + parameters.omegaC*(pn-pm) ) );            
                    }
                }
            } 
            return ret;
        }

        MatrixXcd dgl_chirp(const double t) const {
            if (parameters.chirp_total == 0)
                return MatrixXcd::Zero(parameters.maxStates, parameters.maxStates);
            return 0.5*operatorMatrices.atom_inversion*chirp.get(t);
        }
        MatrixXcd dgl_pulse(const double t) const {
            if (parameters.pulse_amp == 0)
                return MatrixXcd::Zero(parameters.maxStates, parameters.maxStates);
            return 0.5*(operatorMatrices.atom_sigmaplus*pulse.get(t) + operatorMatrices.atom_sigmaminus*std::conj(pulse.get(t)));
        }

        MatrixXcd dgl_phonons(const double t) const {
            return MatrixXcd::Zero(1,1);
        }
        
        void expectationValues(MatrixXcd &rho, const double t) const {
            std::fprintf(fileoutput.fp_atomicinversion,"%.15e\t%.15e\n",t,std::real(dgl_expectationvalue(rho,operatorMatrices.atom_inversion,t))); 
            std::fprintf(fileoutput.fp_photonpopulation,"%.15e\t%.15e\n",t,std::real(dgl_expectationvalue(rho,operatorMatrices.photon_n,t))); 
            std::fprintf(fileoutput.fp_densitymatrix,"%e\t",t);
            for (int j = 0; j < parameters.maxStates; j++)
                std::fprintf(fileoutput.fp_densitymatrix,"%e\t",std::real(rho(j,j)));
            std::fprintf(fileoutput.fp_densitymatrix,"\n");
        }

        bool traceValid(MatrixXcd &rho, double t_hit, bool force = false) {
            double trace = std::real(rho.trace());
            parameters.trace.emplace_back(trace);
            if (trace < 0.99 || trace > 1.01 || force) {
                if (force)
                    fmt::print("{} {} -> trace check failed at t = {} with trace(rho) = {}\n",PREFIX_ERROR,global_message_error_divergent,t_hit,trace);
                terminateMessage = global_message_error_divergent;
                parameters.doSpectrum = 0;
                FILE *fp_trace = std::fopen((parameters.subfolder + "trace.txt").c_str(),"w");
                for (int i = 0; i < parameters.trace.size() && parameters.t_step*1.0*i<t_hit ; i++) {
                    fmt::print(fp_trace,"{:.10e} {:.15e}\n",parameters.t_step*1.0*(i+1),parameters.trace.at(i));
                }
                std::fclose(fp_trace);
                return false;
            } else {
                return true;
            }           
        }

        MatrixXcd dgl_getHamilton(double t) const {
            return dgl_timetrafo( operatorMatrices.H_used + dgl_chirp(t) + dgl_pulse(t), t );
        }

        // Returns order for time direction to use
        int getSolverOrder(int dir = DIR_T) const {
            return (dir == DIR_T ? parameters.orderRKT : parameters.orderRKTau);
        }
};