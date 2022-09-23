#pragma once
#include "global.h"

class Evaluable {
   private:
    std::map<double, Scalar> value_array;
    std::map<double, Scalar> derivative_array;
    std::map<double, Scalar> integral_array;
    std::map<double, Scalar> fourier_value_array;
    Scalar total_maximum = 0;
    Scalar total_minimum = 0;
    size_t size = 0;
    size_t counter_evaluated = 0;
    size_t counter_returned = 0;

   public:
    Evaluable() = default;
    Evaluable( Parameters::input_s &inputs, Parameters &p ) = default;

    virtual Scalar evaluate( double t ) = 0;
    virtual Scalar evaluate_derivative( double t, double dt ) = 0;
    virtual Scalar evaluate_integral( double t, double dt ) = 0;
    virtual Scalar get( double t, bool force_evaluate ) = 0;
    virtual void precalculate() = 0;

    /**
     * @brief Output this evaluable function to a file
     *
     * @param path Output Path
     * @param name Output File Name
     * @param complex If true, both the real and the imaginary part will be output to file
     * @param spectrum If true, the spectrum will be evaluated and output
     */
    void to_file( const std::string &path, const std::string &name, bool complex = false, bool spectrum = false ) {
        // Output Temporal
        auto &file = FileOutput::add_file( path + name );
        if ( complex )
            file << "Time\tReal\tImag\tRealDerivative\tImagDerivative\tRealIntegral\tImagIntegral\n";
        else
            file << "Time\tValue\tDerivative\tIntegral\n";
        for ( const auto &[t, value] : value_array ) {
            const auto deriv = derivative( t );
            const auto integ = integral( t );
            if ( complex )
                file << fmt::format( "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t, std::real( value ), std::imag( value ), std::real( deriv ), std::imag( deriv ), std::real( integ ), std::imag( integ ) );
            else
                file << fmt::format( "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t, std::real( value ), std::real( deriv ), std::imag( integ ) );
        }
        file.close();
        if ( not spectrum )
            return;
        // Output Spectral
        auto &sfile = FileOutput::add_file( path + "fourier_" + name );
        calculate_fourier();
        for ( const auto &[w, value] : fourier_value_array ) {
            file << fmt::format( "{:.8e}\t{:.8e}\t{:.8e}\n", w, std::real( value ), std::imag( value ) );
        }
    }

    // Calculate Derivative
    Scalar derivative( double t, bool force_evaluate = false ) {
        if ( force_evaluate or not derivative_array.contains( t ) )
            return evaluate_derivative( t );
        return derivative_array[t];
    }
    // Calculate Integral
    Scalar integral( double t, bool force_evaluate ) {
        if ( force_evaluate or not integral_array.contains( t ) )
            return evaluate_integral( t );
        return integral_array[t];
    };

    void generate() {
        // Precalculate Values. Can be fixed size or centered oround inputs.
        precalculate();
        // Calculate Maxima/Minima
        std::ranges::for_each( value_array, [&]( const std::pair<double, Scalar> &el ) {
            total_maximum = std::max<double>( total_maximum, std::abs( el.second ) );
            total_minimum = std::min<double>( total_minimum, std::abs( el.second ) );
        } );
        std::ranges::for_each( value_array, [&]( const std::pair<double, Scalar> &el ) {
            // Calculate Derivative at fixed timestep
            derivative( t );
            // Calculate Integral at fixed timestep
            integral( t );
        } );
    }

    // Calculates the Fourier transformation of the input
    void calculate_fourier() {
        // Calculate FourierTransform
    }

    // (Absolute) Minima and Maxima
    Scalar maximum() const {
        return total_maximum;
    }
    Scalar minimum() const {
        return total_minimum;
    }

    /**
     * @brief Returns the number of evaluated and returned values
     *
     * @return std::pair<size_t ,size_t> Evaluated, Returned
     */
    std::pair<size_t, size_t> stats() {
        return { counter_evaluated, counter_returned }
    }
}