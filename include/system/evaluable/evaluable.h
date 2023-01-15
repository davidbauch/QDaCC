#pragma once
#include "typedef.h"
#include "system/fileoutput.h"
#include "misc/log.h"

namespace QDLC {

class Evaluable {
   private:
    using storage_type = std::map<double, Scalar>;
    storage_type value_array;
    storage_type derivative_array;
    storage_type integral_array;
    storage_type fourier_value_array;
    Parameters::universal_config inputs;
    Scalar total_maximum = 0;
    Scalar total_minimum = 0;
    size_t counter_evaluated = 0;
    size_t counter_returned = 0;

   public:
    Evaluable() = default;
    Evaluable( Parameters::universal_config &inputs ) : inputs( inputs ){};
    // Dont generate Evaluate Objects directly using inputs, because those inputs cannot be used.
    Evaluable( Parameters::universal_config &inputs, Parameters &p ) = delete;

    virtual Scalar evaluate( double t ) = 0;
    virtual Scalar evaluate_derivative( double t, double dt = 0 ) = 0;
    virtual Scalar evaluate_integral( double t, double dt = 0 ) = 0;
    /**
     * @brief This function needs to precalculate f(t) for the desired time interval. generate() will then calculate
     * the integral, derivative and fourier transform accordingly.
     *
     */
    virtual ~Evaluable() = default;
    virtual void log() = 0;
    // Calculates the Fourier transformation of the input
    virtual void calculate_fourier( Parameters &p ) = 0;

    void log( const std::string &name ) const {
        Log::L2( "[System-{0}] {0} evaluations/returns: {1}/{2}\n", name, counter_evaluated, counter_returned );
    }

    double get_approximated_dt() {
        return std::real( std::get<1>( *std::next( value_array.begin() ) ) - std::get<1>( *value_array.begin() ) );
    }

    /**
     * @brief Output this evaluable function to a file
     *
     * @param path Output Path
     * @param name Output File Name
     * @param complex If true, both the real and the imaginary part will be output to file
     * @param spectrum If true, the spectrum will be evaluated and output
     */
    void to_file( const std::string &name = "Evaluable", bool complex = false, bool spectrum = false ) {
        // Output Temporal
        auto &file = FileOutput::add_file( name );
        if ( complex )
            file << "Time\tReal\tImag\tRealDerivative\tImagDerivative\tRealIntegral\tImagIntegral\n";
        else
            file << "Time\tValue\tDerivative\tIntegral\n";
        // This should work since the value array is an ordered map
        const auto dt_approx = get_approximated_dt();
        for ( const auto &[t, value] : value_array ) {
            const auto deriv = derivative( t, dt_approx );
            const auto integ = integral( t, dt_approx );
            if ( complex )
                file << fmt::format( "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t, std::real( value ), std::imag( value ), std::real( deriv ), std::imag( deriv ), std::real( integ ), std::imag( integ ) );
            else
                file << fmt::format( "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t, std::real( value ), std::real( deriv ), std::real( integ ) );
        }
        file.close();
        if ( not spectrum )
            return;
        // Output Spectral
        auto &sfile = FileOutput::add_file( name + "_fourier" );
        for ( const auto &[w, value] : fourier_value_array ) {
            sfile << fmt::format( "{:.8e}\t{:.8e}\t{:.8e}\n", w, std::real( value ), std::imag( value ) );
        }
    }

    /**
     * @brief Returns f(t), which is defined by the evaluate(t) function.
     *
     * @param t Current Time
     * @param force_evaluate If true, the function f(t) will be calculated instead of returned from cache.
     * @return Scalar
     */
    Scalar get( double t, bool force_evaluate = false ) {
        if ( force_evaluate or not value_array.contains( t ) ) {
            counter_evaluated++;
#pragma omp critical
            value_array[t] = evaluate( t );
        }
        counter_returned++;
        return value_array[t];
    }
    // Return Derivative
    Scalar derivative( const double t, const double dt = 0, const bool force_evaluate = false ) {
        if ( force_evaluate or not derivative_array.contains( t ) ) {
#pragma omp critical
            derivative_array[t] = evaluate_derivative( t, dt );
        }
        return derivative_array[t];
    }
    // Return Integral
    Scalar integral( const double t, const double dt = 0, const bool force_evaluate = false ) {
        if ( force_evaluate or not integral_array.contains( t ) ) {
#pragma omp critical
            integral_array[t] = evaluate_integral( t, dt );
        }
        return integral_array[t];
    };

    void generate( Parameters &p ) {
        // Precalculate time direction
        double t;
        const std::vector<double> steps = { 0, 0.5 * p.t_step };
        for ( double t1 = p.t_start; t1 <= p.t_end; t1 += p.t_step ) {
            for ( size_t i = 0; i < steps.size(); i++ ) {
                t = t1 + steps[i];
                Scalar val = get( t );
            }
        }
        // Calculate Maxima/Minima
        std::ranges::for_each( value_array, [&]( const std::pair<double, Scalar> &el ) {
            total_maximum = std::max<double>( std::abs( total_maximum ), std::abs( el.second ) );
            total_minimum = std::min<double>( std::abs( total_minimum ), std::abs( el.second ) );
        } );
        std::ranges::for_each( value_array, [&]( const std::pair<double, Scalar> &el ) {
            const auto t = std::get<0>( el );
            // Calculate Derivative at fixed timestep
            derivative( t, p.t_step );
            // Calculate Integral at fixed timestep
            integral( t, p.t_step );
        } );
        // Fourier Transform
        calculate_fourier( p );
    }

    // (Absolute) Minima and Maxima
    Scalar maximum() const {
        return total_maximum;
    }
    Scalar minimum() const {
        return total_minimum;
    }

    size_t size() const {
        return value_array.size();
    }

    /**
     * @brief Acces the fourier value array
     *
     * @param w
     * @param value
     */
    void set_fourier_value( const double w, const Scalar value ) {
        fourier_value_array[w] = value;
    }
    size_t get_fourier_size() const {
        return fourier_value_array.size();
    }

    const Parameters::universal_config &get_inputs() const {
        return inputs;
    }

    /**
     * @brief Returns the number of evaluated and returned values
     *
     * @return std::pair<size_t ,size_t> Evaluated, Returned
     */
    std::pair<size_t, size_t> stats() const {
        return { counter_evaluated, counter_returned };
    }
};
} // namespace QDLC