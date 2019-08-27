#include "system.h"

void System_Parent::init() {
    terminate_message = global_message_normaltermination;
    // Adjusting inputs:
    Timer &timer_systeminit = createTimer("System Initialization");
    logs.level2("System initialization... ");
    timer_systeminit.start();
    if (!init_system()) {
        logs.level2("System initialization failed! Exitting program...\n");
        logs.close();
        exit(EXIT_FAILURE);
    }
    logs.level2("successful. Elapsed time is {}ms\n",timer_systeminit.getWallTime(TIMER_MILLISECONDS));
    timer_systeminit.end();
}

bool System_Parent::traceValid(MatrixXcd &rho, double t_hit, bool force = false) {
    double trace = std::real(rho.trace());
    parameters.trace.emplace_back(trace);
    if (trace < 0.99 || trace > 1.01 || force) {
        if (force)
            fmt::print("{} {} -> trace check failed at t = {} with trace(rho) = {}\n",PREFIX_ERROR,global_message_error_divergent,t_hit,trace);
        terminate_message = global_message_error_divergent;
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