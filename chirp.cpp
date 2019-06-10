#include "chirp.h"
//#include "misc/spline.h"

Chirp::Chirp(Parameters &p) {
    int n = (int)( (p.t_end - p.t_start)/p.t_step*2.0 + 5 );
    chirparray.reserve(n);
    timearray.reserve(n);
    step = p.t_step;
    //if (p.orderRKT == 5 || p.orderRKTau == 5)
    //    steps = { 0, 1./5.*p.t_step, 3./10.*p.t_step, 1./2.*p.t_step, 4./5.*p.t_step, 8./9.*p.t_step };
    //else 
        steps = {0,0.5*p.t_step};
    interpolant = Interpolant( p.chirp_t, p.chirp_y, p.chirp_ddt, p.chirp_type ); 
    generate(p);
}


/* Generate array of energy-values corresponding to the chirp */
void Chirp::generate(Parameters &p) {
    logs.level2("generating... ");
    
    // Array construction

    // Interpolant
    
    // Precalculate 
    double totalChirp = 0;
    double t;
    for (double t1 = p.t_start; t1 < p.t_end+10*p.t_step; t1 += p.t_step/2.0) {
        //for (int i = 0; i < steps.size(); i++ ) {
            /*t = t1;// + steps[i];
            long unsigned int n = 0;
            while (!(t >= inputs.chirp_start[n] && t < inputs.chirp_end[n]+inputs.chirp_onofftime[n]) && n < inputs.chirp_start.size())
                n++;
            double chirpDelta = (inputs.chirp_total[n] != 0 && (inputs.chirp_start[n] != inputs.chirp_end[n])) ? (double)(inputs.chirp_total[n])/(inputs.chirp_end[n]-inputs.chirp_start[n])*p.t_step/2.0 : 0;
            if (chirpDelta != 0) {
                if (t <= inputs.chirp_end[n] && t >= inputs.chirp_start[n]+inputs.chirp_onofftime[n] ){
                    totalChirp += chirpDelta;
                }
                // On:
                else if (t <= inputs.chirp_start[n]+inputs.chirp_onofftime[n] && t >= inputs.chirp_start[n] ){
                    totalChirp += chirpDelta/2. * (1.+std::sin(-M_PI/2. + M_PI* (t-inputs.chirp_start[n])/inputs.chirp_onofftime[n] ));
                }
                // Off:
                else if (t <= inputs.chirp_end[n]+inputs.chirp_onofftime[n] && t >= inputs.chirp_end[n] ){
                    totalChirp += chirpDelta/2. * (1.+std::cos(M_PI* (t-inputs.chirp_end[n])/inputs.chirp_onofftime[n] ));
                }
            }*/

            chirparray.push_back(interpolant.evaluate(t1));
            timearray.push_back(t1);
        //}
    }
    chirparray.shrink_to_fit();
    timearray.shrink_to_fit();
    size = chirparray.size();
    logs.level2("chirparray.size() = {}... ", size);
}

// Create chirparray from scaled parameter chirp
/*void Chirp::generateFromInputFile(Parameters &p) {
    logs.level2("Creating Chirp from scaled input chirp via chirpfile '" + p.subfolder + "'... ");
    FILE *chirpfile = std::fopen((p.subfolder+"chirpfile.txt").c_str(),"r");
    if (!chirpfile) {
        logs.level2("failed!\n");
        generateFromParam(p);
        return;
    }

    double deltaw,t0,t1,chirponoff;
    inputs = Inputs(p);
    logs.level2("Reading chirpfile... ");
    int n = 0;
    while (chirpfile) {
        std::fscanf(chirpfile, "%lf\t%lf\t%lf\t%lf\n", &deltaw,&t0,&t1,&chirponoff);
        if (deltaw + t0 + t1 + chirponoff != -4) {
            inputs.chirp_total.push_back( p.chirp_total*deltaw );    
            inputs.chirp_start.push_back( (p.t_end-p.t_start)*t0 );
            inputs.chirp_end.push_back( (p.t_end-p.t_start)*t1 );
            inputs.chirp_onofftime.push_back( p.chirp_onofftime*chirponoff );
            n++;
        } else {
            break;
        }
    }
    generate(p);
    logs.level2("done! Read %d scaled chirps.\n",n);
}*/

// Create single Chirp from parameters 
/*void Chirp::generateFromParam(Parameters &p) {
    logs.level2("Creating Chirp from input parameters only... ");
    inputs = Inputs(p);
    generate(p);
    logs.level2("done!\n");
}*/

/* 
    Create spline-interpolated chirp
    Chirp is created via total_chirp + spline_delta 
    Inputs: t:      
        time (absolut)
        deltaw: deltaw (relative)
*/
/*void Chirp::generateFromSpline(Parameters &p) {
    logs.level2("Creating Chirp via spline interpolation... ");
    FILE *chirpfile = std::fopen((p.subfolder+"chirpfile.txt").c_str(),"r");
    if (!chirpfile) {
        logs.level2("failed!\n");
        generateFromParam(p);
        return;
    } 
    inputs = Inputs(p);
    generate(p); 
    inputs.chirp_total.clear();
    inputs.chirp_start.clear();
    double t0 = p.chirp_end+p.chirp_onofftime;
    inputs.chirp_total.push_back(0); inputs.chirp_start.push_back(0);
    inputs.chirp_total.push_back(0); inputs.chirp_start.push_back(p.chirp_end + p.chirp_onofftime);
    logs.level2("Reading chirpfile... ");
    int n = 5;
    double t, deltaw;
    while (chirpfile) {
        std::fscanf(chirpfile, "%lf\t%lf\n", &t,&deltaw);
        if (deltaw + t != -2) {
            inputs.chirp_total.push_back( p.chirp_total*deltaw );    
            inputs.chirp_start.push_back( t );
            //fmt::print("t = {}, w = {}\n",t,p.chirp_total*deltaw);
            n++;
        } else {
            break;
        }
    }
    // generating spline interpolation
    logs.level2("Generating spline interpolation... ");
    for(long unsigned int i = 0; i < inputs.chirp_start.size()-2; i++) {
        if (inputs.chirp_start.at(i) >= inputs.chirp_start[i+1]) {
            logs.level2("Error: Input values for time overlap (x[t={}] = {:e} >= x[t={}] = {:e})! Program is probably going to crash... ", i, inputs.chirp_start.at(i), i+1, inputs.chirp_start.at(i));
        }
    }
    tk::spline s;
    s.set_points(inputs.chirp_start, inputs.chirp_total);
    logs.level2("adjusting chirp array... ");
    double last = 0;
    for (double t1 = t0; t1 < p.t_end+step*steps.size(); t1 += step) {
        for (int i = 0; i < steps.size(); i++ ) {
            t = t1 + steps[i];
            if (t <= inputs.chirp_start.at(inputs.chirp_start.size()-1)) {
                last = p.chirp_total + s(t);
            }
            int j = std::floor(t/step-1)*steps.size();
            while (j < chirparray.size() && timearray.at(j) < t) {
                j++;
            }
            chirparray.at(j) = last;
        }
    }
    chirparray.shrink_to_fit();
    timearray.shrink_to_fit();
    size = chirparray.size();
    logs.level2("chirparray.size() = {}... ", size);
    logs.level2("done! Read {} points for chirp interpolation.\n",n);
}*/

void Chirp::fileOutput(std::string filepath, Parameters &p) {
    FILE *chirpfile = std::fopen(filepath.c_str(),"w");
    if (!chirpfile) {
        logs.level2("Failed to open outputfile for chirp!\n");
        return;
    }
    for (long unsigned int i = 0; i < chirparray.size(); i++) {
        std::fprintf(chirpfile,"%.10e\t%.10e\n",timearray.at(i),chirparray.at(i));
    }
    std::fclose(chirpfile);
}

// Per index and rest is lerp delta
// TODO: option for new evaluation instead of precalculating
double Chirp::get(double t) const {
    int i = std::floor(t/step-1)*steps.size();
    while (timearray.at(i) < t) {
        i++;
    }
    if (i < 0 || i >= size) {
        logs.level2("!! Warning: requested chirpvalue at index {} is out of range! chirparray.size() = {}\n",i,chirparray.size());
        i = 0;
    }
    //logs.level2("Requested t = {:.10e}, returned t = {:.10e}\n",t,timearray.at(i));
    return chirparray.at(i);
    //return lerp(chirparray.at(j),chirparray.at(j+1),delta);
}

// Per index
double Chirp::get(int i) const {
    if (i < 0 || i >= size) {
        logs.level2("!! Warning: requested chirpvalue at index {} is out of range!\n",i);
        i = 0;
    }
    return chirparray.at(i);
}