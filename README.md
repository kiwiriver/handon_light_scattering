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
This package was developed in 2014 by Meng Gao (https://github.com/kiwiriver) with supports from Prof. Ping Yang and George Kattawar at Texas A&M University. 

## References
- Born, M. and E. Wolf (1999). Principles of Optics: Electromagnetic Theory of Propagation, Interference and Diffraction of Light.

- Chandrasekhar, S. (2011). Radiative Transfer, Dover Publications, Inc.

- Gao, M., X. Huang, P. Yang and G. W. Kattawar (2013). "Angular distribution of diffuse reflectance from incoherent multiple scattering in turbid media." Applied Optics 52(24): 5869-5879.

- Gao, M. (2012), "Physics of the structural color on the skin of cephalopods", PhD Thesis, Texas A&M University

- Wendisch, M. and P. Yang (2012). Theory of Atmospheric Radiative Transfer.

- Yang, P. and K. N. Liou (1995). "Light scattering by hexagonal ice crystals: comparison of finite-difference time domain and geometric optics models." Journal of the Optical Society of America A 12(1): 162-176.

- Zhai, P.-W., G. W. Kattawar and P. Yang (2008). "Impulse response solution to the three-dimensional vector radiative transfer equation in atmosphere-ocean systems. I. Monte Carlo method." Applied Optics 47(8): 1037-1047.

