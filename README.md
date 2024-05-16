# fullkerr
The fullkerr XSPEC module: continuum emission from within the ISCO 

## More detailed README
See the pdf file README_fullkerr for more detailed instructions. (in prep.)

## Paper
This repository contains the XSPEC model that was developed in [Mummery et al. 2024](https://academic.oup.com/mnras/article/531/1/366/7671518). 

Please cite this paper if you use this model. 

## Code 
* The main program is a fortran90 file ``fullkerr.f90``
* Raytracing modules are included in the fortran90 file ``amodules.f90``, and associated ``X.mod`` files. 
* Files to load this into XSPEC ``load.xcm`` and ``lmodel.dat`` is also included. 


## Full model input parameter description

All of the model parameters are listed below, along with their units and limits of validity. 

| Parameter | Units | range | Description |
 --- | --- | --- | --- |
| $`M `$ |  Solar masses  | $`0 < M `$ | Black hole mass |
| $`a_\star `$ | Dimensionless | $`-1 < a_\star < 1`$ | Black hole spin |
| $`i `$  | Degrees | $`0 < i < 90`$ | Disc-observer inclination angle |
| $`\dot M`$  | $`L_{\rm edd}/c^2 `$  | $`0 < \dot M`$ | Accretion rate |
| $`\delta_{\cal J}`$  | Dimensionless  |  $`0 < \delta_{\cal J} < 1`$ | ISCO stress parameter |
| $`f_1`$  | Dimensionless  |  $`1 < f_1`$ | Colour-correction factor in main body of the disc |
| $`f_2`$  | Dimensionless  |  $`1 < f_2`$ | Colour-correction factor at the ISCO |
| $`\xi`$  | Dimensionless  |  $`0 < \xi`$ | Colour-correction index within the ISCO $`f_c = f_2 (r/r_I)^{-\xi}`$ |
| $`{\tt norm}`$  |  1/kpc$`^2`$   | $`0 < {\tt norm}`$ | Disc-observer distance |

## Loading into XSPEC 
* For use in the XSPEC software (https://heasarc.gsfc.nasa.gov/xanadu/xspec/)
* XSPEC requires the files: ``lmodel.dat``, ``load.xcm``, ``fullkerr.f90``, ``amodules.f90`` and all ``.mod`` files
* ``fullkerr.f90`` contains the actual fortran implementation of the model 
* See https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html for more information on how to load the model into XSPEC
* For any questions about the model, please email andrew.mummery@physics.ox.ac.uk 
