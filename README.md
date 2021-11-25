# Kwandebungdprogramm

# Features

- RK4/5/45  
- Polaron Frame Phonons
    - Double Markov
    - Single Markov, approximations
    - Backwards integral
- Path Integral Phonons  
- G1/2, single and double integrated
- Emission Spectra of any electronic or optical transition (mutliple transitions with "+")
- Indistinguishability of any electronic or optical transition  
- Concurrence of any two electronic or optical transition
- Wigner function of any single electronic or optical transition   

# Basic Syntax for defining a system
... 

# Usage

All numerical parameters support the following units:
Units of Energy:
- `Hz` - Hertz, can also be omitted, meaning any energy related value without a unit is treated as `Hz`
- `eV` - Electron Volt
- `meV` - Milli Electron Volt `= 1E-3eV`
- `mueV` - Mucro Electron Volt `= 1E-6eV`

Units of Time:
- `s` - Seconds, can also be omitted, meaning any time related value without a unit is treated as `s`
- `ns` - Nanoseconds `= 1E-9s`
- `ps` - Picoseconds `= 1E-12s`
- `fs` - Femtoseconds `= 1E-15s`

---

    --file [FILEPATH]
The program can read multiple parameter subsets from a given file. The first line of this file should start with a `#`, followed by the name of the project, which is also the subfolder that will be created for the project. Every other `#` indicates a comment, which can start at any point of the line.

Example:
    
    # NameOfTheProject
    [Parameterset2] # Comment
    # Random Comment
    [Parameterset2] # Comment

***

## System Parameters

    --S []


---
---
## Numerical Parameters

Different numerical Parameters

---
### Grid Resolution
Defines the t-tau grid on wich the G1/G2 correlation functions are evaluated.

    --gridres [NUM]
Defines the Grid timestep `(t1-t0)/NUM` such that exactly `NUM` timesteps are evaluated. The standard value is `-1`, which indicates that all timesteps should be used.

    --grid [...:t-dt:...]
Defines the Grid via several intervals. The syntax for the intervals is defined as `t-dt`, where `t` is the upper limit for where the timestep `dt` should be used. Indicated by the `...`, multiple intervals can be chained via `:`. Note, that this option overwrites `--gridres`.

Examples:

    --gridres 500
Uses 500 timesteps for the correlation grid.

    --grid '10ps-100fs:30ps-200fs:60ps-300fs'
Uses `10ps/100fs + 20ps/200fs + 30ps/300fs = 300` points for the grid, where for the first 10ps a timestep of 100fs is used, between 10ps and 20ps a timestep of 200fs is used, and between 20ps and 30ps a timestep of 300fs is used.

---
### Interpolation Parameters
The program uses interpolation for creating equidistant grids for the G1 and G2 calculations, as well as creating smooth outputs for the raw data that is calculated using varying timesteps. If not specified, the temporal outputs are not interpolated. The correlation functions are always at least linearly interpolated to guarantee all G1 and G2 functions are compatible with each other.

    -interpolate
Enables interpolation of the temporal output data. The standard interpolation method used is the [Monotone Hermite Spline Interpolation](https://jbrd.github.io/2020/12/27/monotone-cubic-interpolation.html) provided by the [ALGLIB](https://www.alglib.net/) package. If this flag ist not passed, the temporal outputs will **NOT** be interpolated.

    --interpolateOrder [orderTime,orderTau]
Defines the interpolation method for both the temporal output data as well as the data used by the G1 and G2 correlation functions. The standard method for interpolation for the correlation functions is a simple [Linear Interpolation](https://en.wikipedia.org/wiki/Linear_interpolation). Both parameters can either be `monotone` or `linear`. The correlation parameter can be omitted if only the temporal interpolation method should be changed.

Examples:

    --interpolateOrder 'monotone,monotone'
Changes the interpolation method for the correlation functions to use the Monotone Interpolation as well. Note that to change to interpolation method for the correlation function, the interpolation method for the temporal direction also has to be defined.

    --interpolateOrder linear
Changes the interpolation method for the temporal output data to use the Linear Interpolation.

---
### Logfile Output
The program generates `logfile.log` which contains information about the input parameters as well as the calculated outputs. While general information is always logged, higher tiers of logging can be enabled, mostly for debugging purposes because the creator of this program refuses to read into actual g++ debugging :)

    -L2
This flag enables subroutine logging, meaning more information about the parameters passed into a function are output to the console and to the logfile.

    -L3
This flag enables deep logging, where all subroutines log everything they do. Soley for debugging.

    -handler
This flag enables handler outputs for the [Threadhandler](https://github.com/davidbauch/Threadhandler), an ancient project.
