# Basic scattering methods for radiative transfer simulations

This package demonstrate basic single and multiple light scattering method.
The simulations require the knowledge of the single scattering properties from individual particles, and can provide the diffuse reflected and transmitted radiance from the multiple scattering among a collection of the particles.
Details description of this package is in the [document (https://github.com/kiwiriver/handon_light_scattering/blob/main/Handon_LightScattering_Practice.pdf)].


## Numerical Methods
Here we demonstrate the basic concepts and numerical simulation methods for both single and multiple scattering methods including the 2D geometrical ray tracing method, and the Monte Carlo radiative transfer method. The Monte Carlo code is also extended to discuss the method of sucessive orders of scattering (SOS).

### 2D Geometrical Ray Tracing
A ray is incident on a sphere, after each reflection and transmission its intensity decreases

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/Ray_Tracing/example/sphere_nr_1.33/ray_path.png" alt="drawing" width="500"/>

The simulated phase function:

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/Ray_Tracing/example/sphere_nr_1.33/phase_function.png" alt="drawing" width="500"/>

### Efficient Monte Carlo Radiative Transfer Simulation
### Successive Order of Scattering (SOS)

## Compiling
Each code is provided with makefile and examples for plot and analysis.

## Authors
This package was developed by Meng Gao (https://github.com/kiwiriver) with support by Prof. Ping Yang at Texas A&M University in 2014. 
