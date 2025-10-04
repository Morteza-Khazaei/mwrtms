(update on Sep. 10, 2014; Sep. 18, 2014)

DMRT-QMS (QCA Mie scattering of Sticky spheres) v0.1

The code supports both active and passive remote sensing of layered snowpack.
It is an implementation of the dense media radiative transfer theory, applying 
the scattering model of QCA Mie of densely packed Sticky spheres. 
The radiative transfer equation is solved using the discrete ordinate method by 
eigen-quadrature analysis.

Ground roughness is partially taken into account for the emissivity in passive remote 
sensing and backscatter in active remote sensing.

In the snowpack description file, each row describes one snow layer,  
starting from top, going downwards. 
Each row contains five columns, specifying the
layer_thickness (cm), density (gm/cc), temperature (K), grain_diameter (cm), and stickiness,
respectively.


Please adapt passive/test_DMRT_QMS_passive.m to deploy the passive code,
and modify active/test_DMRT_QMS_active.m to run the active code.

The following two files are maintained seperated. Check http://www.ee.washington.edu/research/laceo/NMM3D-LUT regularly for updates.
* Common/NMM3D_LUT_NRCS_40degree.dat 	(bare soil backscatter NMM3D look up table)
* Common/NMM3D_LUT_NRCS_40degree_interp.m	(correponding interpolation subroutine)

Copyright: University of Washington, 
Laboratory of Applications and Computations in Electromagnetics and Optics

All these materials could be found at 
http://www.ee.washington.edu/research/laceo/DMRT-QMS.html