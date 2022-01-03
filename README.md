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

# Usage

All numerical parameters support the following units:
Units of Energy:
- `Hz` - Hertz, can also be omitted, meaning the default unit for energy is `Hz`
- `eV` - Electron Volt
- `meV` - Milli Electron Volt `= 1E-3eV`
- `mueV` - Micro Electron Volt `= 1E-6eV`

Units of Time:
- `s` - Seconds, can also be omitted, meaning the default value for time is `s`
- `ns` - Nanoseconds `= 1E-9s`
- `ps` - Picoseconds `= 1E-12s`
- `fs` - Femtoseconds `= 1E-15s`

---
## Inputting Parameters

    --file [FILEPATH]
The program can read multiple parameter subsets from a given file. The first line of this file should start with a `#`, followed by the name of the project, which is also the subfolder that will be created for the project. Every other `#` indicates a comment, which can start at any point of the line.


Example:
    
    # NameOfTheProject
    [Parameterset2] # Comment
    # Random Comment
    [Parameterset2] # Comment

---

## System Parameters
The program accepts two types of systems which will be transformed into a total system by using the tensor product of the single bases. A single electronic system can be tensored with multiple optical systems, where every electronic transition can be coupled to any optical system.

### Electronic Sytem - Fermionic States
The general syntax for a system with a single electronic state reads:
    
    --SE [...;Name:Energy:CoupledTo:DecayScaling:DephasingScaling:PhononCoupling;...]
Multiple electronic states can be chained by using the chain operator `;`. The parameters for a state read as follows:

- `Name` : **Single Uppercase** Character Identifier. **Must** be unique amoung all electronic states.
- `Energy` : Energy of the state. Supports the Units of Energy.
- `CoupledTo` : Comma seperated Single Character Indentifier of any other state that this state is coupled to. Using `-` indicates the state is not coupled to anything.
- `DecayScaling` : Scaling factor (float) for the radiative decay rate for this state
- `DephasingScaling` : Scaling factor (float) for the dephasing rate for this state
- `PhononCoupling` : Scaling factor (float) for the electron-phonon coupling

Example:

    --SE 'G:0:X:0:1:0;X:1.5eV:-:1:1:1'
- Ground State `Name = G` with `Energy = 0`. The groundstate cannot decay radiatively, such that `DecayScaling = 0`. Because the phase relation of this and other states can decay, `DephasingScaling = 1`. This state should not be influenced by phonons with `PhononCoupling = 0`.
- Excited State  `Name = X` with `Energy = 1.5eV`. The excited state can decay radiatively, with `DecayScaling = 1`. Because the phase relation of this and other states can decay, `DephasingScaling = 1`. This state is affected by electron-phonon coupling, hence `PhononCoupling = 1`.

### Photonic System - Bosonic Cavity Resonator
The general syntax for a system with a single optical resonator reads:

    --SO [...;Name:Energy:MaxPhotons:CoupledToTransition:CouplingScaling:DecayScaling;...]
Multiple photonic states can be chained by using the chain operator `;`. The parameters for a state read as follows:

- `Name` : **Single Lowercase** Character Identifier. **Must** be unique amoung all resonators.
- `Energy` : Energy of the resonator. Supports the Units of Energy.
- `MaxPhotons` : Maximum Photon Number for this resonator. 
- `CoupledToTransition` : Comma seperated Electronic Transitions this resonator is coupled to. Uses the *upwards* transition of two states. The transition has to exist in the electronic system.
- `CouplingScaling` : Scaling factor (float) for the transitions. When multiple transitions are provided, an qual number of scaling factors, seperated by commas, have to be defined.
- `DecayScaling` : Scaling factor (float) for the decay rate of photon population inside the resonator.

Example:

    --SE 'c:1.5eV:2:GX:1:1'
- Cavity with `Name = c` with `Energy = 1.5eV`. The maximum number of photons allowed is `MaxPhotons = 2`. The resonator is coupled to the electronic transition `CoupledToTransition = GX` with a coupling scaling of `CouplingScaling = 1`. The photons inside the resonator decay with a scaled rate of `DecayScaling = 1`.

### Further Examples
- Single Exciton with no Cavity. Hilbert Space Dimensions: `2 by 2`.
`--SE 'G:0:X:0:1:0;X:1.5eV:-:1:1:1'`
- Single Exciton coupled to a single mode Resonator with a maximum of two photons. Hilbert Space Dimensions: 6 by 6
`--SE 'G:0:X:0:1:0;X:1.5eV:-:1:1:1'`
`--SO 'c:1.5eV:2:GX:1:1'`     
- Three Level System with cavity at two-photon resonance. No direct Ground-to-Highest State transition is allowed in this system. Phonons scaled doubled for the upper state. Hilbert Space Dimensions: `9 by 9`.
`--SE 'G:0:L:0:1:0;L:1.5eV:U:1:1:1;U:2.8eV:-:1:1:2'`
`--SO 'c:1.6eV:2:GL,LU:1,1:1'`
- [Polarized Biexciton with polarization dependent cavity modes](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.104.085308). The excitons are seperated by a large finestructuresplitting energy of `100mueV` and the biexciton is reduced by a binding energy of `20meV`. The resonators lie at the two-photon resonance of the biexciton and only couple to their specific polarized exciton transitions. Hilbert Space Dimensions: `36 by 36`.
`--SE 'G:0:H,V:0:1:0;H:1.5001eV:Z:1:1:1;V:1.4999eV:Z:1:1:1;Z:2.98eV:-:1:1:1'`
`--SO 'h:1.49eV:2:GH,HZ:1,1:1;v:1.49eV:2:GV,VZ:1,1:1'`
- Polarized Biexciton with polarization dependend bi-modal cavity modes. The two energetically distinguishable resonator modes enhance the ground-exciton as well as the exciton-biexciton transitions. Hilbert Space Dimensions: `324 by 324`.
`--SE 'G:0:H,V:0:1:0;H:1.5001eV:Z:1:1:1;V:1.4999eV:Z:1:1:1;Z:2.98eV:-:1:1:1'`
`--SO 'h:1.5eV:2:GH,HZ:1,1:1;v:1.5eV:2:GV,VZ:1,1:1;k:1.48eV:2:GH,HZ:1,1:1;l:1.48eV:2:GV,VZ:1,1:1'`

### Optical Pulse - Optical Driving of the electronic or optical transitions
The electronic states as well as the resonator modes can be driven by an optical pulse. The general syntax reads:

    --SP '[...;Name:CoupledToTransition:Amp:Freq:Width:Center:Chirp:Type(:GaussAmp);...]'
Multiple transitions can be driven by multiple pulses, and multiple pulses can be chained by using the chain operator `;`. The parameters read as follows:

- `Name` : **Single Lowercase** Character Identifier. **Must** be unique amoung all pulses.
- `CoupledToTransition` : Comma seperated list of electronic transitions or optical resonators this set of pulses is coupled to.
- `Amp` : Amplitude of the pulse. When defined as a Gaussian pulse (see `Type`), this amplitude is assumed to be the pulse area in units of PI. In any other case, this value supports the Units of Energy.
- `Freq` : Frequency of the pulse using `e^(-i*omega)` where omega is the pulse frequency. This parameter supports the Units of Energy.
- `Width` : Temporal pulse width. This parameter supports the Units of Time.
- `Center` : Temporal center of the pulse. This parameter supports the Units of Time.
- `Chirp` : [Pulse Chirp](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.106.166801) rate. This parameter supports the Units of Energy.
- `Type` : Type of the pulse, can either be `gauss`, `gauss_pi` or `cw` (continuous wave)
- `GaussAmp` : **Optional Parameter**. Amplitude of the gaussian function `e^[(t-t0)^A]`, where `A = GaussAmp`. The default value is `GaussAmp = 2`.

Multiple pulses can be applied onto the same transitions by using comma sperated values for `Amp`,`Freq`,`Width`,`Center`,`Chirp`,`Type` (and optionally `GaussAmp`). In the case of `Type = cw`, the parameter `Width` does nothing.

Examples:

    --SP 'p:GX:1:1.5eV:2ps:20ps:0:gauss_pi'
Pulse for a single exciton with `Amp = PulseArea = 1*PI`, `Energy = 1.5eV`, centered at `Center = 20ps` with width `Width = 2ps`. The pulse is not chirped such that `Chirp = 0`. A Gaussian shape is used. Not that with `Type = gauss`, this pulse would have an energy equivalent of `Energy = 1Hz` due to the missing specification for the type. Since no `GaussAmp` was provided, the default value of `2` is used.

    --SP 'p:GH,HZ:5.4:1.48eV:5ps:30ps:0:gauss_pi:12'
Pulse for the horizontally polarized biexciton transition. A `GaussAmp = 12` is used to approximate a flat-top Gaussian shape.

    --SP 'p:GH,HZ:5.4,1,1:1.48eV,1.5eV,1.5eV:5ps,1ps,1ps:30ps,50ps,60ps:0,0,0:gauss_pi,gauss_pi,gauss_pi'
The same pulse as before, followed by two short PI-pulses at later times.

    --SP 'p:GL:1:1.5eV:2ps:20ps:0:gauss_pi;q:LU:1:1.4eV:2ps:24ps:0:gauss_pi'
Two pulses for different electronic transitions for the three level system.

    --SP 'p:h:0.15meV:1.5eV:0:69fs:0:cw'
A continuous wave pulse with an amplitude of `Amp = 0.15meV`, frequency `Freq = 1.5eV` driving the resonator `h`. Make sure `h` actually exists as a resonator, or the program will crash. The `Center = 69fs` equals a phase shift relative to zero.

### Electronic Chirp - Temporal Detuning of the electronic states
The energy of the electronic states can be shifted over time by the electronic chirp. The general syntax reads:

    --SC [...;Name:CoupledTo:AmpFactors:Amps:Times:DDTs:Type;...]
Multiple levels can be changed by the same chirp, and multiple chirps can be chained using the chain operator `;`. The parameters read as follows:

- `Name` : **Single Lowercase** Character indentifier. **Must** be unique amoung all chirps.
- `CoupledTo` : Comma seperated list of uppercase electronic state identifiers this chirp should influence.
- `AmpFactors` : Comma seperated list of scaling factors for the identifiers listet prior.
- `Amps` : Comma seperated list of amplitudes.
- `Times` : Comma seperated list of times.
- `DDTs` : Comma sperated list of derivatives. This list of parameters is only used if `hermite` is chosen as the chirp generation method.
- `Type` : Type of chirp, can be one of the following: `linear`, `spline`, `monotone` and `hermite`, which is the interpolation method the chirp class uses to generate curves from the initial points specified.


Examples:

    --SC 'c:X:1:0,1meV,0,0:50ps,100ps,150ps,200ps:0,0,0,0:monotone'
A chirp for a single state `CoupledTo = X` with an amplitude scaling of `AmpFactors = 1`. The curve will follow the list of points provided using a monotone interpolation method.

    --SC 'c:H,V,Z:1,1,2:0,1meV,0,0:50ps,100ps,150ps,200ps:0,0,0,0:monotone'
The same chirp but for the biexciton system. The two single excitons are shifted with factor `1`, while the biexciton shifts twice with factor `2` because of its energy definition.

---

## Temporal Settings
The program calculates the temporal dynamics of the provided system in a single given interval using either the fixed-stepsize RK4/RK5 or the variable order RK45 method. The timer interval is given by

    --time [START] [END] [STEP]
and supports the Units of Time. If `END` is set to `-1`, the temporal dynamics are calculated until the density matrix ground state is reached by at least `99.9%`. The density matrix diagonal index for the ground state can be provided by
    
    --groundstate [DM_INDEX]
where the default value for `DM_INDEX` is `0`.

If only a single parameter is to be changed, the following parameters can be used:

    --tstart [START]
    --tend [END]
    --tstep [STEP]

Examples:

    --time 0 2ns 100fs
Calculate from `0` to `2ns` using a timestep of `100fs`.

    --tstep 200fs
Use a timestep of `200fs`. If no `--tend [END]` is provided, the calculation will continue until the system is converged in its ground state. **Caution**: If the groundstate can never be reached, the program will fail to terminate.

---

## Electron-Phonon Coupling
The electron-phonon coupling can be calculated using two different approaches; The [Polaron approach](https://journals.aps.org/prb/supplemental/10.1103/PhysRevB.96.201201) utilizing at least one Markov approximation, where preferably two are used to significantly reduce calculation times; And the more suffisticated [Path-Integral Formalism](https://journals.aps.org/prb/supplemental/10.1103/PhysRevB.96.201201), where no such approximations need to be done.

### General Syntax for Phonons
In general, phonon contributions are not included unless a temperature is provided:

    --temperature [TEMP]
Where `TEMP` is any value greater than zero. For values smaller than zero, the electron-phonon coupling will be disabled.

The spectral properties of the phonons are defined by

![](https://latex.codecogs.com/svg.image?%5Cleft.J_%7Bij%7D%5Comega=%5Csum_%7Bk%7D%5Cgamma_k%5Ei%5Cgamma_k%5Ej%5Cdelta(%5Comega-%5Comega_j)=%5Cfrac%7B%5Comega%5E3%7D%7B4%5Cpi%5E2%5Crho%5Chbar%20c_s%5E5%7D%5Cleft(D_ee%5E%7B-%5Comega%5E2a_e%5E2/(4c_s%5E2)%7D-D_he%5E%7B-%5Comega%5E2a_h%5E2/(4c_s%5E2)%7D%5Cright)%5E2=%5Calpha_p%5Comega%5E3e%5E%7B-%5Cfrac%7B%5Comega%5E2%7D%7B2%5Comega_b%7D%7D%5Cright.).

While the specific spectral shape varies slightly for electrons (`e`) and holes (`h`), the algorythm used in this program assumes a unified distribution combinding both electron and hole parameters into two parameters:

    --phononalpha [COUPLING]
The effective phononcoupling. The standard value for GaAs QDs lies around `COUPLING = 0.03ps^-1`.

    --phononwcutoff [CUTOFF]
The energy for which maximum coupling occurs. The standard value for GaAs QDs lies around `CUTOFF = 1meV`

## Correlation functions
Two different types of correlation functions can be evaluated: `G1(t,tau)` and `G2(t,tau)`. From these two, quantum properties such as the `indistinguishability` or the `two-photon concurrence` can be evaluated. If the evaluation of any property that needs either G1 or G2 is provided, the correlation functions are calculated automatically. In this case, explicitely specifying to calculate the corresponding G1 or G2 functions is not neccessary. If the output of the corresponding function is needed, specification becomes neccessary again.

The G1 function is calculated via  

![](https://latex.codecogs.com/svg.latex?G_i%5E%7B%281%29%7D%28t%2Ct%27%29%20%3D%5Ctext%7BTr%7D%5B%7B%5Crho%27%28t%27%29%5Chat%7Bb%7D_i%5E%5Cdagger%280%29%7D%5D%5Ctext%7B%20with%20%7D%5Crho%27%280%29%20%3D%20%5Chat%7Bb%7D_i%280%29%5Crho%28t%29).

The G2 function is calculated via  

![](https://latex.codecogs.com/svg.latex?G_%7Bi%2Cj%7D%5E%7B%282%29%7D%28t%2Ct%27%29%20%3D%20%5Ctext%7BTr%7D%5B%7B%5Crho%27%28t%27%29%5Chat%7Bb%7D_i%5E%5Cdagger%280%29%5Chat%7Bb%7D_j%280%29%7D%5D%5Ctext%7B%20with%20%7D%5Crho%27%280%29%20%3D%20%5Chat%7Bb%7D_j%280%29%5Crho%28t%29%5Chat%7Bb%7D_i%5E%5Cdagger%280%29).

The general Syntax for both functions reads:

    --GF [Modes:Orders:Integration]
The parameters read as follows:

- `Modes` : Comma sperated list of modes to calculate the function for. Miltiple transitions can be combined via the `+` operator.
- `Orders` : Comma sperated list of the order corresponding to the specified mode.
- `Integration` : `0` for the `G(t,tau)` function, `1` for both the `t` and `tau` integrated functions and `2` for both. Note that the `G(t,tau)` output will be a possibly giant matrix resulting in filesizes of several hundreds of megabytes.

To specify the resolution for the `G1` and `G2` grid, see `Numerical Parameters` -> `Grid Resolution`

Examples:

    --GF 'h,GX:1,1:2,2'
Calculating the G1 functions for the resonator mode `h` and the radiative electronic mode `GX`. Both of order `1` and both output as a matrix and integrated.

    --GF 'h,h:1,2:2:2'
Calculating the G1 and G2 functions for the resonator mode `h`.

---

## Emission Spectrum
The emission characteristics of any transition can be calculated by using the [Eberly-Wodkiewicz](https://www.osapublishing.org/josa/abstract.cfm?uri=josa-67-9-1252) Spectrum and is calculated via  

![](https://latex.codecogs.com/svg.latex?%5Cmathcal%7BS%7D_i%28t_%5Ctext%7Bmax%7D%2C%5Comega%29%20%3D%20%5CRe%20%5Cint_%7B0%7D%5E%7Bt_%5Ctext%7Bmax%7D%7D%5Cint_0%5E%7Bt_%5Ctext%7Bmax%7D-t%7DG_i%5E%7B%281%29%7D%28t%2Ct%27%29e%5E%7B-i%5Comega%20t%27%7D%5Ctext%7Bd%7Dt%27%5Ctext%7Bd%7Dt).

If the neccessary `G1` function is not yet available, it will be calculated automatically. The fourier transformation is brute-forced for a set number of points. The general syntax for spectrum calculation reads:

    --GS [Mode:Center:Resolution:Points]

The parameters read as follows:
- `Mode` : Comma sperated list of resonator modes or electronic transitions or any superposition of any of the available transitions. Multiple transitions can be added using the `+` operator.
- `Center` : Comma sperated list of spectral centers for the Fourier transformation of the corresponding mode. This parameter supports the Units of Energy.
- `Range` : Comma sperated list of spectral ranges for the Fourier transformation of the corresponding mode. This paramter suppports the Units of Energy.
- `Points` : Number of points the Fourier transformation should be avaluated on. The `dw` step of the integral is then determined by `dw=2w_Range/N_Points`.

Examples:

    --GS 'h:1.5eV:500mueV:1000'
Calculating the emission spectrum for the resonator mode `h` centered around `1.5meV` expanding for `500mueV` in both directions. The Fourier transformation is evaluated for `1000` points, resulting in $d\omega=1\mu eV$.

    --GS 'h,GX,h+GX:1.5eV,1.5eV,1.5eV:1meV,1meV,1meV:2000,2000,2000'
Calculating the emission spectra for the resonator mode `h`, the radiative transition `GX` and the superposition of both `h+GX` using the same parameters as before.

---

## Single Photon Indistinguishability and Visibility
The single photon [HOM indistinguishability](https://en.wikipedia.org/wiki/Hong%E2%80%93Ou%E2%80%93Mandel_effect) is a figure of merit for the consistentcy and quality of a single photon source. The indistinguishability is calculated [via](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.045309)  

![](https://latex.codecogs.com/svg.latex?%5Cmathcal%7BI%7D_%7Bi%7D%20%3D%201-p_%7Bc%2Ci%7D%20%3D%201-%5Cfrac%7B%5Cint_0%5E%7Bt_%5Ctext%7Bmax%7D%7D%5Cint_0%5E%7B%7Bt_%5Ctext%7Bmax%7D%7D-t%7D%202G_%7B%5Ctext%7BHOM%7D%2Ci%7D%5E%7B%282%29%7D%28t%2Ct%27%29dt%27dt%7D%7B%5Cint_0%5E%7Bt_%5Ctext%7Bmax%7D%7D%5Cint_0%5E%7B%7Bt_%5Ctext%7Bmax%7D%7D-t%7D%5Cleft%28%202G_%7B%5Ctext%7Bpop%7D%2Ci%7D%5E%7B%282%29%7D%28t%2Ct%27%29-%20%7C%5Clangle%5Chat%7Ba%7D_i%28t&plus;t%27%29%5Crangle%20%5Clangle%5Chat%7Ba%7D_i%5E%5Cdagger%28t%29%5Crangle%7C%5E2%20%5Cright%29dt%27dt%7D)  
with  

![](https://latex.codecogs.com/svg.latex?G_%7B%5Ctext%7BHOM%7D%2Ci%7D%5E%7B%282%29%7D%28t%2Ct%27%29%20%3D%20%5Cfrac%7B1%7D%7B2%7D%5Cleft%28%20G_%7B%5Ctext%7Bpop%7D%2Ci%7D%5E%7B%282%29%7D%28t%2Ct%27%29&plus;G%5E%7B%282%29%7D_%7Bi%2Ci%7D%28t%2Ct%27%29%20-%20%7CG_i%5E%7B%281%29%7D%28t%2Ct%27%29%7C%5E2%20%5Cright%29)  
and  

![](https://latex.codecogs.com/svg.latex?G_%7B%5Ctext%7Bpop%7D%2Ci%7D%5E%7B%282%29%7D%20%3D%20%5Clangle%5Chat%7Bb%7D_i%5E%5Cdagger%5Chat%7Bb%7D_i%5Crangle%28t%29%5Clangle%5Chat%7Bb%7D_i%5E%5Cdagger%5Chat%7Bb%7D_i%5Crangle%28t&plus;t%27%29).

The single photon visibility is a more simple figure of merit for the photon quality and is calculated [via](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.233605)    

![](https://latex.codecogs.com/svg.latex?%5Cmathcal%7BV%7D_%7Bi%7D%20%3D%20%5Cfrac%7B%5Cint_0%5E%7Bt_%5Ctext%7Bmax%7D%7D%5Cint_0%5E%7B%7Bt_%5Ctext%7Bmax%7D%7D-t%7D%202%7CG_%7Bi%7D%5E%7B%281%29%7D%28t%2Ct%27%29%7C%5E2%5Ctext%7Bd%7Dt%27%5Ctext%7Bd%7Dt%7D%7B%5Cleft%28%5Cint_0%5E%7Bt_%5Ctext%7Bmax%7D%7D%5Clangle%5Chat%7Ba%7D%5E%5Cdagger_i%28t%29%5Chat%7Ba%7D_i%28t%29%5Crangle%5Ctext%7Bd%7Dt%5Cright%29%5E2%7D).  

If any of the neccessary `G1` or `G2` correlation functions is not yet available, it will be automatically calculated. The general syntax reads:

    --GI [Modes]

The parameters read as follows:
- `Modes` : Comma seperated list of modes to calculate both the indistinguishability as well as the visibility for. The superposition of multiple modes can be calculated using the `+` operator.

Since both properties require at least `G1(t,tau)`, they can always both be calculated without significantly increasing the required calculation times.

Examples:

    --GI 'h'
Calculate both properties for a single resonator mode `h`.

    --GI 'h,GH,h+GH'
Calculate both properties for the resonator mode `h`, the radiative transition `GH` and the superposition of both.

---

## Two-Photon Concurrence - Measurement for the entanglement of two transitions

The [Two-Photon Concurrence](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.80.2245) is calculated by reducing the densitymatrix from its total Hilbert Space onto a two photon densitymatrix. The elements of this densitymatrix are calculated using the `G2` correlation function. A single element is calculated by integrating `G2_ijkl`, where `i`,`j`,`k` and `l` are one of the two modes that the concurrence is evaluated for.

The Two-Photon Matrix results in  

![](https://latex.codecogs.com/svg.latex?%5Crho_%5Ctext%7B2phot%7D%3D%5Cbegin%7Bbmatrix%7D%5Crho_%7Biiii%7D%20%26%200%20%26%200%20%26%20%5Crho_%7Biijj%7D%20%5Ccr%200%20%26%20%5Crho_%7Bijij%7D%20%26%20%5Crho_%7Bijji%7D%20%26%200%20%5Ccr%200%20%26%20%5Crho_%7Bjiij%7D%20%26%20%5Crho%5E*_%7Bijij%7D%20%26%200%20%5Ccr%20%5Crho%5E*_%7Biijj%7D%20%26%200%20%26%200%20%26%20%5Crho_%7Bjjjj%7D%5Cend%7Bbmatrix%7D%20%5Ctext%7B%20with%20%7D%20%5Crho_%7Bijkl%7D%20%3D%20%5Cint_0%5E%7Bt_%5Ctext%7Bmax%7D%7D%5Cint_0%5E%7Bt_%5Ctext%7Bmax%7D-t%7DG_%7Bij%2Ckl%7D%5E%7B%282%29%7D%28t%2Ct%27%29%5Ctext%7Bd%7Dt%27%5Ctext%7Bd%7Dt)  

Because of symmetry of the operators, only 6 of the 8 `G2` functions have to be calculated. The concurrence is then evaluated by calculating the Eigenvalues of  

![](https://latex.codecogs.com/svg.latex?R%3D%5Csqrt%7B%5Csqrt%7B%5Crho_%5Ctext%7B2ph%7D%7D%28%5Csigma_y%20%5Cotimes%20%5Csigma_y%29%5Crho_%5Ctext%7B2ph%7D%5E*%28%5Csigma_y%20%5Cotimes%20%5Csigma_y%29%5Csqrt%7B%5Crho_%5Ctext%7B2ph%7D%7D%7D)  

where

![](https://latex.codecogs.com/svg.latex?%5Cmathcal%7BC%7D%20%3D%20%5Ctext%7Bmax%7D%5Cleft%5C%7B%7B0%2C%5Clambda_4-%5Clambda_3-%5Clambda_2-%5Clambda_1%7D%5Cright%5C%7D%20%5Ctext%7B%20with%20%7D%20%5Clambda_i%20%5Ctext%7B%20%5Ctextit%7Bi%7D%27th%20Eigenvalue%20of%20%7D%20R)  

If any of the neccessary `G2` correlation functions is not yet available, it will be automatically calculated. Note that if the Hilbert space of the system does not allow `G2` to be non-zero (e.g only a single photon allowed in a resonator mode), the concurrence will be undefined. The general syntax reads:

    --GC [Mode1-Mode2]

The parameters read as follows:
- `Mode1-Mode2` : Comma seperated list of modes. Two modes are connected using the `-` operator, indicating the Concurrence between those two modes is calculated. The superposition of multiple modes can be calculated using the `+` operator.

Examples:

    --GC 'h-v'
Calculating the concurrence for the two resonator modes of the polarization dependent biexciton.

    --GC 'h+GH+HZ-v+GV+VZ'
Calculating the concurrence for every cavity or radiative mode of the polarization dependent biexciton.

---

## Wigner Function
The Wigner function is calculated using the weighted superposition of Laguerre polynomials. A reduced denistymatrix for a given mode is calculated by using a partial trace over the total Hilbert space. This densitymatrix is then assumed to be in a Fock-base containing the different weights for the Laguerre polynomials.

The general syntax to calculate the Wigner function reads:

    --GW [Modes,X(,Y,Resolution,Skip)]

The parameters read as follows:
- `Modes` : Comma seperated list of (resonator) modes to calculate the Wigner function for
- `X` : Comma sperated list of X-Resolutions, where the total interval will be `-X:X`
- `Y` : Optional Parameter, comma seperated list of Y-Resolutions. As default it is assumed that X=Y
- `Resolution` : Optional Parameter, comma sperated list of resolutions, where `dx,dy = 2X/Resolution, 2Y/Resolution`. The default value is `100`.
- `Skip` : Number of timesteps to skip. Higher values will result in faster calculations. The default value is `1`.

Examples:  

    --GW 'h:5'
Calculate the Wigner function using the parameterset `h:5:5:100:1`.

    --GW 'h:4:6:200:2'
Calculating the Wigner function using `-X:X = -4:4`, `-Y:Y = -6:6` with `200` Points for each axis, resulting in `dx = 0.04` and `dy = 0.06`. The Wigner function is only evaluated for every second timestep, even if those timesteps are not equidistant (e.g. for RK45 or a modified grid without interpolation).

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
### Interpolation
The program uses interpolation for creating equidistant grids for the G1 and G2 calculations, as well as creating smooth outputs for the raw data that is calculated using varying timesteps. If not specified, the temporal outputs are not interpolated. The correlation functions are always at least linearly interpolated onto equidistant grids to guarantee all G1 and G2 functions are compatible with each other.

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
