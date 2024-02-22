# WaveCascade
Exact calculations of cylindrical waveguides and coaxial waveguides connected in series.

Electromagnetic waves are assumed to travel from left to right inside a waveguide. The waveguide may consist of several sections, and each section is either a circular (cylindrical) waveguide or a circular (cylindrical) coaxial waveguide. 

The main functions used to compile such a waveguide are _scattering_matrix_coaxials()_ and _scattering_matrix_mixed()_. These functions are described below, while details on the physics can be found in the attached PDF report, specifically chapter 2. See the simple example scripts _example_coaxial_fsweep.m_ and _example_hw3_waveguide.m_ for basic usage of the code. More details on code usage can be found in Chapter 5 of the PDF report, including explanations of the two example scripts.

The function _scattering_matrix_coaxials()_ calculates the scattering matrix between a junction of two coaxials with shared inner conductor but different outer wall radii, see Chapter 2.1 in PDF. Input parameters are
- Operating frequency of the electromagnetic wave.
- Radius of both coaxials' inner conductor. NOTE: Two connecting coaxials are assumed to share the same inner conductor.
- Radius of the outer walls of the first coaxial.
- Radius of the outer walls of the second coaxial.
- Length of the first coaxial. Only relevant if the user wishes to include the propagation from coaxial start to end (junction).
- Length of the second coaxial. Only relevant if the user wishes to include the propagation from coaxial start (junction) to its end.
- Number of modes to include in the calculation.
- Optional input parameters include (i) printing the cut-off frequencies for each waveguide segment (ii) whether to include propagation in first and/or second segment.

The function _scattering_matrix_mixed()_ calculates the scattering matrix between a junction of a coaxial and a circular waveguide with shared outer radii, see Chapter 2.3 in PDF. Input parameters are
- Operating frequency of the electromagnetic wave.
- Radius of coaxial's inner conductor. 
- Radius of the outer walls of both waveguides. NOTE: Two connecting waveguides of mixed type are assumed to share the same outer wall radius.
- Length of the first waveguide. Only relevant if the user wishes to include the propagation from coaxial start to end (junction).
- Length of the second waveguide. Only relevant if the user wishes to include the propagation from coaxial start (junction) to its end.
- Number of modes to include in the calculation.
- Junction case of type 1 or 2. A type 1 junction describes going from a coaxial waveguide to a circular waveguide. A type 2 junction describes going from a circular waveguide to a coaxial.
- Optional input parameters include (i) printing the cut-off frequencies for each waveguide segment (ii) whether to include propagation in first and/or second segment.

A 'coaxial' is meant as a circular waveguide with an inner conductor. A circular waveguide is meant as a hollow metallic waveguide with circular cross-section geometry. The outer walls of both waveguide types are perfect electric conductors.
