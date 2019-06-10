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
        double chirp_total, chirp_delta, chirp_start, chirp_end; // FIXME: REDUNDANCY
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
            //TODO: inputs als --system -sys, --pulse -p, --chirp -c, --spectrum -s
            logs.wrapInBar("Conversion of input variables",LOG_SIZE_FULL,LOG_LEVEL_2,LOG_BAR_0); logs.level2("\n");
            logs.level2("Parsing input variables... ");
            int index = 1;
            t_start             = 0.0;
            t_end               = getNextInput<double>(arguments,"t_end",index);
            t_step              = getNextInput<double>(arguments,"t_step",index);
            omegaE              = getNextInput<double>(arguments,"omegaE",index);
            omegaC              = getNextInput<double>(arguments,"omegaC",index);
            g                   = getNextInput<double>(arguments,"g",index);
            k                   = getNextInput<double>(arguments,"k",index);
            gamma_pure          = getNextInput<double>(arguments,"gamma_pure",index);
            gamma               = getNextInput<double>(arguments,"gamma",index);
            chirp_t             = getNextInputVector<double>(arguments,"chirp_time",index); // TODO: ddt sollte optional sein
            chirp_y             = getNextInputVector<double>(arguments,"chirp_yval",index);
            chirp_ddt           = getNextInputVector<double>(arguments,"chirp_diff",index);
            chirp_type          = getNextInputString(arguments,"chirp_file_type",index);
            maxPhotonNumber     = getNextInput<int>(arguments,"maxPhotonNumber",index);
            rho_0               = getNextInput<int>(arguments,"rho_0",index);
            pulse_center        = getNextInput<double>(arguments,"pulse_center",index);
            pulse_amp           = getNextInput<double>(arguments,"pulse_amp",index);
            pulse_freq          = getNextInput<double>(arguments,"pulse_freq",index);
            pulse_sigma         = getNextInput<double>(arguments,"pulse_sigma",index);
            pulsetype           = getNextInputString(arguments,"pulsetype",index);
            double rkType       = getNextInput<double>(arguments,"rungekuttatype",index);
            doSpectrum          = getNextInput<double>(arguments,"doSpectrum",index);
            akf_everyXIt        = getNextInput<int>(arguments,"akf_everyXIt",index);
            akf_omega_center    = getNextInput<double>(arguments,"akf_omega_center",index);
            akf_omega_range     = getNextInput<double>(arguments,"akf_omega_range",index);
            akf_skip_omega      = getNextInput<int>(arguments,"akf_skip_omega",index);
            doInteractionPicture= getNextInput<int>(arguments,"doInteractionPicture",index);
            doRWA               = getNextInput<int>(arguments,"doRWA",index);
            akf_maxThreads      = getNextInput<int>(arguments,"akf_maxThreads",index);
            advancedLogging     = getNextInput<int>(arguments,"advancedLogging",index);
            subfolder           = getNextInputString(arguments,"subfolder",index);
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
                    t_step = 2./8.*M_PI/std::max(init_rabifrequenz,max_rabifrequenz);
                if (!doRWA)
                    t_step = 2./48.*M_PI/std::max(init_rabifrequenz,max_rabifrequenz);
            }

            // Calculate the maximum dimensions for operator matrices (max states)
            maxStates           = 2*(maxPhotonNumber+1);

            // Calculate stuff for RK
            maxItMainRK         = (int)std::ceil((t_end-t_start)/t_step);
            orderRKT            = (int)(rkType/10); // 4 or 5
            orderRKTau          = ((int)rkType)%10; // 4 or 5

            // Calculate values for chirp: (maybe) //TODO: Redundant, //REMOVE
            chirp_total = 1;
            chirp_start = 0;
            chirp_end = 1E-10;
            chirp_onofftime = 1E-10;

            // Mandatory: rescaling chirp ddt into chirp/ps
            for (long unsigned i = 0; i < chirp_ddt.size(); i++)
                chirp_ddt.at(i) = chirp_ddt.at(i)*1E12;

            // Initial and maximum system detuning (not taking into account custom chirps)
            init_detuning     = (omegaE - omegaC);
            max_detuning      = (init_detuning + chirp_total > init_detuning) ? init_detuning + chirp_total : init_detuning;
            
            // Chirp per iteration
            chirp_delta = chirp_ddt.at(std::distance(chirp_ddt.begin(), std::max_element(chirp_ddt.begin(), chirp_ddt.end(), abs_compare) ));
            //chirp_delta       = //(chirp_total != 0 && (chirp_start != chirp_end)) ? (double)(chirp_total)/(chirp_end-chirp_start)*t_step : 0;

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

            // Logging and handlerstrings:
            outputHandlerStrings=advancedLogging;
            advancedLogging = std::floor(advancedLogging/10.);
            outputHandlerStrings = outputHandlerStrings%10;
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
                logs("Total Chirp w_chirp = {} Hz -> {} eV\n",chirp_total,Hz_to_eV(chirp_total));
                logs("Chirp start t_chirp_start = {} s\n",chirp_start);
                logs("Chirp end t_chirp_end = {} s\n",chirp_end);
                logs("Chirp per Iteration w_chirp_delta = {} Hz/it -> {} mueV/it\n",chirp_delta,Hz_to_eV(chirp_delta)*1E6);
                logs("-> Chirp per ps = {} Hz/ps -> {} mueV/ps\n",chirp_delta*1E-12/t_step,Hz_to_eV(chirp_delta*1E-12/t_step)*1E6);
                logs("Chirp on/off time is {} s",chirp_onofftime);
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
            logs.level2("OutputHandlerStrings: {}, advancedlogging: {}\n",outputHandlerStrings,advancedLogging);
        }
};