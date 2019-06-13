static void help() {
    fmt::print("--time [start] [end] [step]\n\t--tstart [start]\n\t--tend [end]\n\t--tstep [step]\n");
    fmt::print("--system [omegaE] [omegaC] [g] [kappa] [gamma_pure] [gamma] else standard values (1.52eV,1.52eV,66mueV,66mueV,3mueV,1mueV) are used\n\t--we [omegaE]\n\t--wc [omegaC]\n\t--coupling [g]\n\t--kappa [kappa]\n\t--gamma [gamma]\n\t--gammapure [gamma_pure]\n");
    fmt::print("--chirp ['[Array Time]'] ['[Array Y]'] ['[Array d/dt]'] [type]\n\t--chirpT ['[Array Time]']\n\t--chirpY ['[Array Y]']\n\t--chirpDDT ['[Array d/dt]']\n\t--chirpType [type] where type = monotone, hermite, linear, spline\n");
    fmt::print("--pulse [Center] [Amplitude] [Frequency] [Sigma] [Type]\n\t-pulse for standard pulse\n\t--pulseCenter [Center]\n\t--pulseAmp [Amplitude]\n\t--pulseFreq [Frequency]\n\t--pulseSigma [Sigma]\n\t--pulseType [Type] where Type = cw, gauss, gauss_pi\n");
    fmt::print("--dimensions [maximum Photons] [Initial state]\n\t--maxPhotons [maximum Photons]\n\t--initState [Initial state], has to be smaller than (2*n+1)\n");
    fmt::print("--spectrum [Iteration Skips] [Center] [Range] [Omega Skips] enables spectrum\n\t-spectrum enables spectrum centered at cavity\n\t--specXIt [Iterations skips (int)]\n\t--specCenter [Center]\n\t--specRange [Range]\n\t--specSkip [0/1] True/False\n");
    fmt::print("-RK5 enables Runge Kutta of order 5 for T and Tau direction\n\t-RK5T enables Runge Kutta of order 5 for T direction\n\t-RK5Tau enables Runge Kutta of order 5 for Tau direction\n");
    fmt::print("-noInteractionpic disables Interaction picture - enabled by default\n-noRWA disables rotating wave approximation - enabled by default\n-noHandler disables handerl strings and enables loadbar output (for console)\n");
    fmt::print("--Threads [number] number of threads to use for both AKF and Spectrum integral calculation\n\n");
}

class Parameters {
    public:
        // Mandatory parameters, move to parent class:
        // Time variables
        double t_start, t_end, t_step;
        // System Dimensions
        int maxStates, maxPhotonNumber;
        // Starting state:
        int rho_0;
        // Non Numeric parameters:
        std::string subfolder;
        int doSpectrum;
        int advancedLogging;
        int outputHandlerStrings;

        // Non mandatory parameters, dependant on system chosen:
        // System Parameters
        double omegaC;
        double omegaE;
        double g;
        double k;
        double gamma_pure, gamma;
        int doRWA, doInteractionPicture;
        // Calculated System properties:
        double init_detuning, max_detuning, init_rabifrequenz, max_rabifrequenz;
        // Chirp and Pulse properties:
        double pulse_center, pulse_amp, pulse_freq, pulse_sigma, chirp_onofftime;
        double chirp_total;
        std::vector<double> chirp_t,chirp_y,chirp_ddt;

        std::string pulsetype, chirp_type;

        // Runtime parameters and other stuff
        int maxItMainRK; 
        int maxItTotal; 
        int akf_maxThreads;
        int orderRKT, orderRKTau;
        std::vector<double> trace;

        int akf_everyXIt,akf_skip_omega;
        int akf_spec_max_w;
        double akf_deltaWmax;
        std::vector<std::complex<double>> akf_spec_out;
        double akf_omega_center, akf_omega_range;
        std::vector<double> akf_spec_possiblew;

        double scaling;

        Parameters() {};
        Parameters(const std::vector<std::string> &arguments) {
            Timer &timer = createTimer("Parameters");
            timer.start();

            // Parse parameters and convert to SI, Input parameters 
            logs.wrapInBar("Conversion of input variables",LOG_SIZE_FULL,LOG_LEVEL_2,LOG_BAR_0); logs.level2("\n");
            logs.level2("Parsing input variables... ");
            int index = 1;
            bool legacy = false;
            // Looking for legacy input, usefull because the python handler doesnt understand named inputs :(
            if ( (index = vec_find_str("-legacy",arguments)) != -1) {
                t_start             = 0.0;
                t_end               = getNextInput<double>(arguments,"t_end",++index);
                t_step              = getNextInput<double>(arguments,"t_step",index);
                omegaE              = getNextInput<double>(arguments,"omegaE",index);
                omegaC              = getNextInput<double>(arguments,"omegaC",index);
                g                   = getNextInput<double>(arguments,"g",index);
                k                   = getNextInput<double>(arguments,"k",index);
                gamma_pure          = getNextInput<double>(arguments,"gamma_pure",index);
                gamma               = getNextInput<double>(arguments,"gamma",index);
                
                chirp_t             = getNextInputVector<double>(arguments,"chirp_time",index);
                chirp_y             = getNextInputVector<double>(arguments,"chirp_yval",index);
                chirp_ddt           = getNextInputVector<double>(arguments,"chirp_diff",index);
                chirp_type          = getNextInputString(arguments,"chirp_file_type",index);
                
                pulse_center        = getNextInput<double>(arguments,"pulse_center",index);
                pulse_amp           = getNextInput<double>(arguments,"pulse_amp",index);
                pulse_freq          = getNextInput<double>(arguments,"pulse_freq",index);
                pulse_sigma         = getNextInput<double>(arguments,"pulse_sigma",index);
                pulsetype           = getNextInputString(arguments,"pulsetype",index);
                
                maxPhotonNumber     = getNextInput<int>(arguments,"maxPhotonNumber",index);
                rho_0               = getNextInput<int>(arguments,"rho_0",index);
                
                doSpectrum          = getNextInput<double>(arguments,"doSpectrum",index);
                akf_everyXIt        = getNextInput<int>(arguments,"akf_everyXIt",index);
                akf_omega_center    = getNextInput<double>(arguments,"akf_omega_center",index);
                akf_omega_range     = getNextInput<double>(arguments,"akf_omega_range",index);
                akf_skip_omega      = getNextInput<int>(arguments,"akf_skip_omega",index);
                
                doInteractionPicture= getNextInput<int>(arguments,"doInteractionPicture",index);
                doRWA               = getNextInput<int>(arguments,"doRWA",index);
                akf_maxThreads      = getNextInput<int>(arguments,"akf_maxThreads",index);
                double rkType       = getNextInput<double>(arguments,"rungekuttatype",index);
                advancedLogging     = getNextInput<int>(arguments,"advancedLogging",index);
                subfolder           = getNextInputString(arguments,"subfolder",index);
                orderRKT            = (int)(rkType/10); // 4 or 5
                orderRKTau          = ((int)rkType)%10; // 4 or 5
                legacy              = true;
            }
            
            
            // Look for --time, if not found, standard values are used (t0 = 0, t1 = 1ns, deltaT = auto)
            if ( (index = vec_find_str("--time", arguments)) != -1) {
                t_start         = 0.0;
                t_end           = getNextInput<double>(arguments,"t_end",++index);
                t_step          = getNextInput<double>(arguments,"t_step",index);
            } else if(!legacy) {
                t_start         = 0.0;
                t_end           = convertParam<double>("1.0ns");
                t_step          = -1;
            }
             // Look for single parameter corrections
            if ( (index = vec_find_str("--tstart", arguments)) != -1) {
                t_start         = getNextInput<double>(arguments,"t_start single",++index);
            }
            if ( (index = vec_find_str("--tend", arguments)) != -1) {
                t_end           = getNextInput<double>(arguments,"t_end single",++index);
            }
            if ( (index = vec_find_str("--tstep", arguments)) != -1) {
                t_step          = getNextInput<double>(arguments,"t_step single",++index);
            }

            // Look for --system, if not found, standard system is used (g=66mueV, k=66mueV, gamma_pure = 3mueV, gamma = 1mueV)
            if ( (index = vec_find_str("--system", arguments)) != -1) {
                omegaE          = getNextInput<double>(arguments,"omegaE",++index);
                omegaC          = getNextInput<double>(arguments,"omegaC",index);
                g               = getNextInput<double>(arguments,"g",index);
                k               = getNextInput<double>(arguments,"k",index);
                gamma_pure      = getNextInput<double>(arguments,"gamma_pure",index);
                gamma           = getNextInput<double>(arguments,"gamma",index);
            } else if(!legacy) {
                omegaE          = convertParam<double>("1.52eV");
                omegaC          = convertParam<double>("1.52eV");
                g               = convertParam<double>("66mueV");
                k               = convertParam<double>("66mueV");
                gamma_pure      = convertParam<double>("3mueV");
                gamma           = convertParam<double>("1mueV");
            }
            // Look for single parameter corrections
            if ( (index = vec_find_str("--we", arguments)) != -1) {
                omegaE          = getNextInput<double>(arguments,"omegaE single",++index);
            }
            if ( (index = vec_find_str("--wc", arguments)) != -1) {
                omegaC          = getNextInput<double>(arguments,"omegaC single",++index);
            }
            if ( (index = vec_find_str("--coupling", arguments)) != -1) {
                g               = getNextInput<double>(arguments,"g single",++index);
            }
            if ( (index = vec_find_str("--kappa", arguments)) != -1) {
                k               = getNextInput<double>(arguments,"k single",++index);
            }
            if ( (index = vec_find_str("--gammapure", arguments)) != -1) {
                gamma_pure      = getNextInput<double>(arguments,"gamma_pure single",++index);
            }
            if ( (index = vec_find_str("--gamma", arguments)) != -1) {
                gamma           = getNextInput<double>(arguments,"gamma single",++index);
            }

            // Look for --chirp, if not found, standard system is used (no chirp, everything zero)
            if ( (index = vec_find_str("--chirp", arguments)) != -1) {
                chirp_t         = getNextInputVector<double>(arguments,"chirp_time",++index);
                chirp_y         = getNextInputVector<double>(arguments,"chirp_yval",index);
                chirp_ddt       = getNextInputVector<double>(arguments,"chirp_diff",index);
                chirp_type      = getNextInputString(arguments,"chirp_file_type",index);    
            } else if(!legacy) {
                chirp_t         = {t_start,t_end};
                chirp_y         = {0.0,0.0};
                chirp_ddt       = {0.0,0.0};
                chirp_type      = "monotone";
            }
            // Look for single parameter corrections
            if ( (index = vec_find_str("--chirpT", arguments)) != -1) {
                chirp_t         = getNextInputVector<double>(arguments,"chirp_time single",++index);
            }
            if ( (index = vec_find_str("--chirpY", arguments)) != -1) {
                chirp_y         = getNextInputVector<double>(arguments,"chirp_yval single",++index);
            }
            if ( (index = vec_find_str("--chirpDDT", arguments)) != -1) {
                chirp_ddt       = getNextInputVector<double>(arguments,"chirp_ddt single",++index);
            }
            if ( (index = vec_find_str("--chirpType", arguments)) != -1) {
                chirp_type      = getNextInputString(arguments,"chirp_type single",++index);
            }
            
            // Look for --pulse, if not found, standard system is used (no pulse, everything zero)
            if ( (index = vec_find_str("--pulse", arguments)) != -1) {
                pulse_center    = getNextInput<double>(arguments,"pulse_center",++index);
                pulse_amp       = getNextInput<double>(arguments,"pulse_amp",index);
                pulse_freq      = getNextInput<double>(arguments,"pulse_freq",index);
                pulse_sigma     = getNextInput<double>(arguments,"pulse_sigma",index);
                pulsetype       = getNextInputString(arguments,"pulsetype",index);
            } else if(!legacy) {
                pulse_center    = 0.0;
                pulse_amp       = 0.0;
                pulse_freq      = 0.0;
                pulse_sigma     = 1.0;
                pulsetype       = "cw";
            }
            if ( (index = vec_find_str("-pulse", arguments)) != -1) {
                pulse_amp       = 1.0;
                pulse_freq      = omegaE;
                pulse_sigma     = convertParam<double>("20ps");
                pulse_center    = 10.0*pulse_sigma;
                pulsetype       = "gauss_pi";
            }
            // Look for single parameter corrections
            if ( (index = vec_find_str("--pulseCenter", arguments)) != -1) {
                pulse_center    = getNextInput<double>(arguments,"pulse_center single",++index);
            }
            if ( (index = vec_find_str("--pulseAmp", arguments)) != -1) {
                pulse_amp       = getNextInput<double>(arguments,"pulse_amp single",++index);
            }
            if ( (index = vec_find_str("--pulseFreq", arguments)) != -1) {
                pulse_freq      = getNextInput<double>(arguments,"pulse_freq single",++index);
            }
            if ( (index = vec_find_str("--pulseSigma", arguments)) != -1) {
                pulse_sigma     = getNextInput<double>(arguments,"pulse_sigma single",++index);
            }
            if ( (index = vec_find_str("--pulseType", arguments)) != -1) {
                pulsetype       = getNextInputString(arguments,"pulsetype single",++index);
            }

            // Look for --dimensions, if not found, standard system is used (maxphotons = 0, starting state = |g,0>)
            if ( (index = vec_find_str("--dimensions", arguments)) != -1) {
                maxPhotonNumber = getNextInput<int>(arguments,"maxPhotonNumber",++index);
                rho_0           = getNextInput<int>(arguments,"rho_0",index);
            } else {
                maxPhotonNumber = 1;
                rho_0           = 0;
            }
            // Look for single parameter corrections
            if ( (index = vec_find_str("--maxPhotons", arguments)) != -1) {
                maxPhotonNumber = getNextInput<int>(arguments,"maxPhotonNumber single",++index);
            }
            if ( (index = vec_find_str("--initState", arguments)) != -1) {
                rho_0           = getNextInput<int>(arguments,"rho_0 single",++index);
            }

            // Look for (-RK4), -RK5, (-RK4T), (-RK4Tau), -RK5T, -RK5Tau
            orderRKT        = 4;
            orderRKTau      = 4;
            if ( (index = vec_find_str("-RK5", arguments)) != -1) {
                orderRKT        = 5;
                orderRKTau      = 5;
            }
            if ( (index = vec_find_str("-RK5T", arguments)) != -1) {
                orderRKT        = 5;
            }   
            if ( (index = vec_find_str("-RK5Tau", arguments)) != -1) {
                orderRKTau      = 5;
            }

            // Look for --spectrum, if not found, no spectrum is evaluated
            if ( (index = vec_find_str("--spectrum", arguments)) != -1) {
                doSpectrum      = 1;
                akf_everyXIt    = getNextInput<int>(arguments,"akf_everyXIt",++index);
                akf_omega_center= getNextInput<double>(arguments,"akf_omega_center",index);
                akf_omega_range = getNextInput<double>(arguments,"akf_omega_range",index);
                akf_skip_omega  = getNextInput<int>(arguments,"akf_skip_omega",index);
            } else if(!legacy) {
                doSpectrum      = 0;
                akf_everyXIt    = 1;
                akf_omega_center= omegaC;
                akf_omega_range = 0.0;
                akf_skip_omega  = 0;
            }
            if ( (index = vec_find_str("-spectrum", arguments)) != -1) {
                doSpectrum      = 1;
                akf_everyXIt    = 1;
                akf_omega_center= omegaC;
                akf_omega_range = (g+k+gamma+gamma_pure)*10.0;
                akf_skip_omega  = 0;
            }
            // Look for single parameter corrections
            if ( (index = vec_find_str("--specXIt", arguments)) != -1) {
                akf_everyXIt    = getNextInput<int>(arguments,"akf_everyXIt single",++index);
            }
            if ( (index = vec_find_str("--specCenter", arguments)) != -1) {
                akf_omega_center= getNextInput<double>(arguments,"akf_omega_center single",++index);
            }
            if ( (index = vec_find_str("--specRange", arguments)) != -1) {
                akf_omega_range = getNextInput<double>(arguments,"akf_omega_range single",++index);
            }
            if ( (index = vec_find_str("--specSkip", arguments)) != -1) {
                akf_skip_omega  = getNextInput<int>(arguments,"akf_skip_omega single",++index);
            }

            // Look for other parameters
            if ( (index = vec_find_str("-nointeractionpic", arguments)) != -1) {
                doInteractionPicture = 0;
            } else if(!legacy) {
                doInteractionPicture = 1;
            }
            if ( (index = vec_find_str("-noRWA", arguments)) != -1) {
                doRWA = 0;
            } else if(!legacy) {
                doRWA = 1;
            }
            if ( (index = vec_find_str("--Threads", arguments)) != -1) {
                akf_maxThreads = getNextInput<int>(arguments,"akf_maxThreads",++index);
            } else if(!legacy) {
                akf_maxThreads = 1;
            }
            if ( (index = vec_find_str("-noHandler", arguments)) != -1) {
                outputHandlerStrings = 0;
            } else {
                outputHandlerStrings = 1;
            }

            subfolder           = arguments.back();

            logs.level2("Done! Elapsed time is {}ms.\nProcessing input variables...",timer.getWallTime()*1E3);

            // Calculate/Recalculate some parameters:
            // Adjust pulse area if pulsetype is "gauss_pi"
            if (pulsetype.compare("gauss_pi") == 0) {
                pulse_amp = M_PI/(std::sqrt(2.0*M_PI)*pulse_sigma);
            }

            // Calculate Rabi frequencies:
            init_rabifrequenz   = rabiFrequency(omegaE-omegaC, g, (int)(rho_0/2.0));
            max_rabifrequenz    = rabiFrequency(omegaE-omegaC, g, (int)(rho_0/2.0) + 1);
            
            // Calculate minimum step necessary to resolve Rabi-oscillation if step=-1
            if (t_step == -1) {
                if (doRWA)
                    t_step = std::min(2./8.*M_PI/std::max(std::max(init_rabifrequenz,max_rabifrequenz),g+k+gamma+gamma_pure),1E-12);
                if (!doRWA)
                    t_step = std::min(2./8.*M_PI/std::max(std::max(init_rabifrequenz,max_rabifrequenz), omegaE+omegaC),1E-12);
            }

            // Calculate the maximum dimensions for operator matrices (max states)
            maxStates           = 2*(maxPhotonNumber+1);

            // Calculate stuff for RK
            maxItMainRK         = (int)std::ceil((t_end-t_start)/t_step);

            // Mandatory: rescale chirp ddt into chirp/ps
            for (long unsigned i = 0; i < chirp_ddt.size(); i++)
                chirp_ddt.at(i) = chirp_ddt.at(i)*1E12;
            // Calculate values for chirp:
            chirp_total = vec_max(chirp_y);

            // Initial and maximum system detuning (not taking into account custom chirps)
            init_detuning     = (omegaE - omegaC);
            max_detuning      = (init_detuning + chirp_total > init_detuning) ? init_detuning + chirp_total : init_detuning;

            // Adjust/calculate frequency range for spectrum
            akf_spec_max_w      = maxItMainRK/( (akf_everyXIt-1)*(1-akf_skip_omega) + 1 );
            if (akf_omega_center == -1)
                akf_omega_center  = omegaC;
            if (akf_omega_range == -1)
                akf_omega_range  = (std::abs(max_detuning)+g+k/2.)*3.;
            for(int w=0;w<akf_spec_max_w;w++) {
                akf_spec_possiblew.push_back(akf_omega_center-(akf_omega_range) + w/((double)maxItMainRK)*((akf_everyXIt-1)*(1-akf_skip_omega)+1)*( 2.*(akf_omega_range) ));
                akf_spec_out.push_back(0);
            }

            // Calculate total number of iterations necessary //FIXME: redundant
            maxItTotal        = maxItMainRK;
            if (doSpectrum) {
                int curIt = 1;
                for (double ii=t_start+t_step; ii<t_end; ii+=t_step) {
                    if (curIt%akf_everyXIt == 0) {
                        for (double iii=ii+t_step; iii<t_end; iii+=t_step) {
                            maxItTotal++;
                        }
                        curIt = 1;   
                    } else
                        curIt += 1;
                }
            }
            trace.reserve(maxItMainRK+5);

            // Done
            timer.end();
            logs.level2("Done! Elapsed time is {}ms\n\n",timer.getWallTime()*1E3);
        }
        void log() {
            logs.wrapInBar("Parameters"); logs("\n");
            logs.wrapInBar("Time borders and t_step",LOG_SIZE_HALF,LOG_LEVEL_1,LOG_BAR_1);
                logs("Timeborder left t_start = {} s\n",t_start);
                logs("Timeborder right t_end = {} s\n",t_end);
                logs("Timeborder t_step delta t = {} s\n",t_step);
                logs("Time iterations (main loop) = {}\n",maxItMainRK);
                logs("Total time iterations = {}\n\n",maxItTotal);//, (int)ceil(maxItMainRK/2.*maxItMainRK/((double)akf_everyXIt)) );
            logs.wrapInBar("System Parameters",LOG_SIZE_HALF,LOG_LEVEL_1,LOG_BAR_1);
                logs("Energy Level difference |g><g| - |e><e| = {} Hz -> {} eV -> {} nm\n",omegaE,Hz_to_eV(omegaE),Hz_to_wavelength(omegaE));
                logs("Cavity Frequency w_c = {} Hz -> {} eV -> {} nm\n",omegaC,Hz_to_eV(omegaC),Hz_to_wavelength(omegaC));
                logs("Coupling strengh g = {} Hz -> {} mueV\n",g,Hz_to_eV(g)*1E6);
                logs("Photon loss rate k = {} Hz -> {} mueV -> Q = {:.2f}\n",k,Hz_to_eV(k)*1E6,omegaC/k);
                logs("Atomic dephasing rate gamma_pure = {} Hz -> {} mueV\n",gamma_pure,Hz_to_eV(gamma_pure)*1E6);
                logs("RAD rate gamma = {} Hz -> {} mueV\n",gamma,Hz_to_eV(gamma)*1E6);
                logs("Initial state rho0 = |{},{}> with maximum number of {} photons\n\n", (rho_0%2 == 0 ? "g" : "e"), (int)std::floor(rho_0/2), maxPhotonNumber );
            logs.wrapInBar("Excitation Pulse",LOG_SIZE_HALF,LOG_LEVEL_1,LOG_BAR_1);
            if (pulse_center != -1 && pulse_amp != 0)
                logs("Exiting system at t_0 = {} with amplitude {} ({}meV), frequency {}eV ({}) and FWHM {}\n\n",pulse_center,pulse_amp,Hz_to_eV(pulse_amp)*1E3,pulse_freq,Hz_to_eV(pulse_freq),pulse_sigma*(2*std::sqrt(2*std::log(2))));
            else
                logs("Not using pulse to exite system\n\n");
            logs.wrapInBar("Energy Chirp",LOG_SIZE_HALF,LOG_LEVEL_1,LOG_BAR_1);
            if (chirp_total != 0) {
                for(int i = 0; i < chirp_t.size()-1; i++) {
                    logs("Chirp between t0 = {}ps and t1 = {}ps: {}mueV -> average rate {}mueV/ps\n",chirp_t.at(i),chirp_t.at(i+1), chirp_y.at(i+1)-chirp_y.at(i), (chirp_y.at(i+1)-chirp_y.at(i))*1e12/(chirp_t.at(i+1)-chirp_t.at(i)));
                }
                if (chirp_type.compare("none") != 0) 
                    logs("\nChirpfile of type '"+chirp_type+"' is used!");
            } else 
                logs("Not using chirp");
            logs("\n\n");
            
            logs.wrapInBar("Caluclated Frequencies");
                logs("\nInitial system detuning = {} Hz -> {} mueV\n",init_detuning,Hz_to_eV(init_detuning)*1E6);
                logs("Maximum system detuning = {} Hz -> {} mueV\n",max_detuning,Hz_to_eV(max_detuning)*1E6);
                logs("Initial Rabi Frequencies = {} Hz -> {} mueV\n",init_rabifrequenz,Hz_to_eV(init_rabifrequenz)*1E6);
                logs("Maximum Rabi Frequencies = {} Hz -> {} mueV\n\n",max_rabifrequenz,Hz_to_eV(max_rabifrequenz)*1E6);
            int works = 1;
            if ( (init_rabifrequenz != 0.0) && (3.*t_step > 2.*M_PI/init_rabifrequenz) )
                works = 0;
            else if (max_rabifrequenz != 0.0 && 3.*t_step > 2.*M_PI/max_rabifrequenz)
                works = 0;
            if (!works) {
                fmt::print("{} WARNING: Step may be too small to resolve predicted oscillation: dT needed vs dT: {:.10e} < {:.10e}\n",PREFIX_WARNING,2./3.*M_PI/std::max(init_rabifrequenz,max_rabifrequenz),t_step);   
                logs("WARNING: Step may be too small to resolve predicted oscillation: \n-> delta T needed: {:.10e} \n-> delta T used: {:.10e}\n\n",2./3.*M_PI/std::max(init_rabifrequenz,max_rabifrequenz),t_step);   
            }
            logs.wrapInBar("Spectrum");
            logs("\nCenter Frequency: {} Hz -> {} eV\n",akf_omega_center,Hz_to_eV(akf_omega_center));
            logs("Frequency Range: +/- {} Hz -> +/- {} mueV\n",akf_omega_range,Hz_to_eV(akf_omega_range)*1E6);
            logs("Skipping {} steps in grid\n\n",akf_everyXIt-1);
            
            logs.wrapInBar("Numerics");
            logs("\nOrder of Runge-Kutta used: Time: RK{}, Spectrum: RK{}\n",orderRKT,orderRKTau);
            logs("Use rotating wave approximation (RWA)? - {}\n",((doRWA == 1) ? "YES" : "NO"));
            logs("Use interaction picture for calculations? - {}\n",((doInteractionPicture == 1) ? "YES" : "NO"));
            logs("Threads used for AFK? - {}\n",akf_maxThreads);
            logs("Used pulsetype? - "+pulsetype+"\n");
            logs("\n");
            logs.wrapInBar("Program Log:",LOG_SIZE_FULL, LOG_LEVEL_2); logs("\n");
            logs.level2("OutputHandlerStrings: {}\n",outputHandlerStrings);
        }
};