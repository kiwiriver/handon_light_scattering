# Basic scattering methods for radiative transfer simulations

This package demonstrate basic single and multiple light scattering method.
The simulations require the knowledge of the single scattering properties from individual particles, and can provide the diffuse reflected and transmitted radiance from the multiple scattering among a collection of the particles.
Details description of this package is in the document (https://github.com/kiwiriver/handon_light_scattering/blob/main/Handon_LightScattering_Practice.pdf).


## Numerical Methods
Here we demonstrate the basic concepts and numerical simulation methods for both single and multiple scattering methods including the 2D geometrical ray tracing method, and the Monte Carlo radiative transfer method. The Monte Carlo code is also extended to discuss the method of sucessive orders of scattering (SOS).

### 2D Geometrical Ray Tracing
A ray is incident on a sphere, after each reflection and transmission its intensity decreases

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/Ray_Tracing/example/ray_path.png" alt="drawing" width="400"/>

The simulated phase function:

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/Ray_Tracing/example/sphere_nr_1.33/phase_function.png" alt="drawing" width="400"/>

### Efficient Monte Carlo Radiative Transfer Simulation
The diffuse radiance simulated by Monte Carlo method in a Rayleigh scattering plane parallel system with a normal optical depth of 1.0.

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/Monte_Carlo/example/rayleigh_tau_1.000_a_1.0/radiance.png" alt="drawing" width="400"/>


### Successive Order of Scattering (SOS)
The diffuse radiance contributed from the first 1, 2, 4, 8 orders of scattering (corresponding to the dot plot upward from bottom) for a turbid medium with an optical depth 1

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/SOS/example/sos.png" alt="drawing" width="400"/>


## Compiling
Each code is provided with makefile and examples for plot and analysis.

## Authors
This package was developed by Meng Gao (https://github.com/kiwiriver) with support by Prof. Ping Yang at Texas A&M University in 2014. 
