#pragma once

#include <complex>
#include <cmath>
#include "typedef.h"

class Parameter {
   private:
    double value_si;
    double value_scaled;
    double scaling;
    int scale_type;
    void update() {
        if ( scaling != 0.0 ) {
            value_scaled = value_si * scaling;
        } else {
            value_scaled = value_si;
        }
    }

   public:
    const static int UNIT_SI = 1;
    const static int UNIT_ENERGY_EV = 2;
    const static int UNIT_ENERGY_MEV = 3;
    const static int UNIT_ENERGY_MUEV = 4;
    const static int UNIT_WAVELENGTH_NM = 5;
    const static int UNIT_TIME_NS = 6;
    const static int UNIT_TIME_PS = 7;
    const static int UNIT_TIME_FS = 8;
    const static int SCALE_ENERGY = 0;
    const static int SCALE_TIME = 1;

    Parameter() { scale_type = -1; };
    Parameter( double val, double scale ) : Parameter() { set( val, scale ); };
    Parameter( double val ) : Parameter( val, 0.0 ){};
    Parameter( Parameter &other ) : Parameter( other.value_si, other.scaling ){};
    Parameter( const Parameter &other ) : Parameter( other.value_si, other.scaling ){};
    template <typename T>
    Parameter( std::complex<T> val ) : Parameter( std::real( val ) ){};

    double get() const {
        return value_scaled;
    }
    double set( double newval, double scale = 0.0 ) {
        value_si = newval;
        scaling = scale;
        update();
        return value_scaled;
    }
    double setScale( double new_scaling, int scale = SCALE_ENERGY ) {
        if ( scale == SCALE_ENERGY )
            new_scaling = 1.0 / new_scaling;
        scale_type = scale;
        return set( value_si, new_scaling );
    }

    double getSI( const int index = UNIT_SI ) {
        switch ( index ) {
            case UNIT_ENERGY_EV:
                return value_si * 6.582119516885722624E-16;
                break;
            case UNIT_ENERGY_MEV:
                return value_si * 6.582119516885722624E-13;
                break;
            case UNIT_ENERGY_MUEV:
                return value_si * 6.582119516885722624E-10;
                break;
            case UNIT_WAVELENGTH_NM:
                return 299792458.0 * 2.0 * QDACC::Math::PI / value_si * 1E9;
                break;
            case UNIT_TIME_NS:
                return value_si * 1E9;
                break;
            case UNIT_TIME_PS:
                return value_si * 1E12;
                break;
            case UNIT_TIME_FS:
                return value_si * 1E15;
                break;
            default:
                return value_si;
                break;
        }
    }

    // Operator Overloads:
    operator double() const {
        return get();
    }
    template <typename T>
    operator std::complex<T>() const {
        return std::complex<T>( get(), 0 );
    }
    operator double() {
        return get();
    }
    template <typename T>
    operator std::complex<T>() {
        return std::complex<T>( get(), 0 );
    }
    double operator+( Parameter &other ) {
        return get() + other.get();
    }
    double operator-( Parameter &other ) {
        return get() - other.get();
    }
    double operator*( Parameter &other ) {
        return get() * other.get();
    }
    double operator/( Parameter &other ) {
        return get() / other.get();
    }
    bool operator<( Parameter &other ) {
        return get() < other.get();
    }
    bool operator>( Parameter &other ) {
        return get() > other.get();
    }
    double operator<=( Parameter &other ) {
        return get() <= other.get();
    }
    double operator>=( Parameter &other ) {
        return get() >= other.get();
    }
    friend std::ostream &operator<<( std::ostream &stream, Parameter &param ) {
        return stream << param.get();
    }
};

template <class T>
inline std::complex<T> operator+( const Parameter &p1, const std::complex<T> &p2 ) {
    return p2 + p1.get();
}
template <class T>
inline std::complex<T> operator+( const std::complex<T> &p1, const Parameter &p2 ) {
    return p1 + p2.get();
}
// Subtraction with arb. numerical class
template <class T>
inline std::complex<T> operator-( const Parameter &p1, const std::complex<T> &p2 ) {
    return p2 - p1.get();
}
template <class T>
inline std::complex<T> operator-( const std::complex<T> &p1, const Parameter &p2 ) {
    return p1 - p2.get();
}
// Multiplication with arb. numerical class
template <class T>
inline std::complex<T> operator*( const Parameter &p1, const std::complex<T> &p2 ) {
    return p2 * p1.get();
}
template <class T>
inline std::complex<T> operator*( const std::complex<T> &p1, const Parameter &p2 ) {
    return p1 * p2.get();
}
// Division with arb. numerical class
template <class T>
inline std::complex<T> operator/( const Parameter &p1, const std::complex<T> &p2 ) {
    return p2 / p1.get();
}
template <class T>
inline std::complex<T> operator/( const std::complex<T> &p1, const Parameter &p2 ) {
    return p1 / p2.get();
}

// template <class T>
// T std::max( const T &p1, const Parameter &p2 ) {
//     return std::max( p1, p2.get() );
// }
// template <class T>
// T std::max( const Parameter &p1, const T &p2 ) {
//     return std::max( p2, p1.get() );
// }
//
// template <class T>
// T std::min( const T &p1, const Parameter &p2 ) {
//     return std::min( p1, p2.get() );
// }
// template <class T>
// T std::min( const Parameter &p1, const T &p2 ) {
//     return std::min( p2, p1.get() );
// }