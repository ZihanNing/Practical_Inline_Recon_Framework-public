**Package for MRI reconstruction. Supports:**
  - non-Cartesian regridding
  - iterative SENSE (Cartesian or non-Cartesian)
  - iterative SENSE for higher order models that may include: 
    - a B0 map
    - time varying spherical harmonics of phase accrual
    
**Before Usage**
  - run setPath.m to add all the directories to the Matlab path
    
**Includes a basic image viewer and ROI drawing tool called "bview"**
  - enter "bview" in the Matlab command prompt to use
    
**Tips**
  - step through the demos in the demos folder
    - recommended order is demo_regridding, demo_regSENSE, demo_highOrder

**Notes**
  - Matlab gpuArray functionality is used whenever possible, so a compatible GPU is strongly recommended

**Acknowledgement**
  - Cite as: Baron, Corey A. (2021, February 2). MatMRI: A GPU enabled package for model based MRI image reconstruction (Version 0.1.00). Zenodo. http://doi.org/10.5281/zenodo.4495477
  - Additionally, please reference the following works for usage of the below functions:
    - findDelAuto: Dubovan PI, Baron CA. Model-based determination of the synchronization delay between MRI and trajectory data. Magn Reson Med. 2022 Sep 26. doi:10.1002/mrm.29460.  PMID: 36161333.
    - nufftOp: Baron CA, Dwork N, Pauly JM, Nishimura DG. Rapid compressed sensing reconstruction of 3D non-Cartesian MRI. Magn. Reson. Med. 2018;79:2685–2692.
    - sampHighOrder: Wilm BJ, Barmet C, Pruessmann KP. Fast higher-order MR image reconstruction using singular-vector separation. IEEE Trans Med Imaging. 2012 Jul;31(7):1396-403. doi: 10.1109/TMI.2012.2190991. Epub 2012 Mar 14. Erratum in: IEEE Trans Med Imaging. 2012 Sep;31(9):1833.

(c) 2020, Corey Baron
