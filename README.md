# Basic scattering methods for radiative transfer simulations

This package demonstrates basic single and multiple light scattering methods.
Details description of this package is in the document (https://github.com/kiwiriver/handon_light_scattering/blob/main/Handon_LightScattering_Practice.pdf).


## Numerical Methods
The following methods are included:

### 2D Geometrical Ray Tracing
A ray is incident on a sphere, after each reflection and transmission its intensity decreases

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/Ray_Tracing/example/ray_path.png" alt="drawing" width="400"/>

The simulated phase function:

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/Ray_Tracing/example/sphere_nr_1.33/phase_function.png" alt="drawing" width="400"/>

### Efficient Monte Carlo Radiative Transfer Simulation
The diffuse radiance simulated by Monte Carlo method in a Rayleigh scattering plane parallel system with various optical depths.

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/Monte_Carlo/example/mc.png" alt="drawing" width="400"/>


### Successive Order of Scattering (SOS)
The diffuse radiance contributed from the first 1, 2, 4, 8 orders of scattering (corresponding to the dot plot upward from bottom) for a turbid medium with an optical depth 1

<img src="https://github.com/kiwiriver/handon_light_scattering/blob/main/SOS/example/sos.png" alt="drawing" width="400"/>


## Compiling
Each code is provided with makefile and examples for plot and analysis.

## Authors
This package was developed by Meng Gao (https://github.com/kiwiriver) with support by Prof. Ping Yang at Texas A&M University in 2014. 
