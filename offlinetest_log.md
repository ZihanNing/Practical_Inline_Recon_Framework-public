# offline test log
## 1. offline reconstruction (without inline implementation)
### 1.1 SENSE
1. SENSE recon with SWI (using ACS)
   demo data set: RAW_SWI_forSENSE_integratedACS.dat 
   - status: code passed successfully 
   - results: checked! (x and xPh saved)

### 1.2 AlignSENSE
1. AlignSENSE recon with T1 MPRAGE (using ACS)
   demo data set: RAW_T1_MPRAGE_forAlignSENSE_seperateACS.dat
   - status: code passed successfully 
   - results: checked! (An-Ve generated)
2. AlignSENSE recon with SWI small FOV (using ACS)
   demo data set: RAW_SWI_forAlignSENSE_smallFOV_seperateACS.dat
   - status: code passed successfully
   - results: checked!
3. AlignSENSE: recon with SWI small FOV (using exterREF)
   demo data set: RAW_SWI_forAlignSENSE_smallFOV_seperateACS.dat
   demo data set: RAW_external_REF_forAlignSENSE_alignedSWIsmallFOV.dat
   - status: code passed successfully
   - results: checked!

## 2. offline test of the inline implementation version
### 2.1 SENSE
1. SENSE recon with SWI (using ACS)
   - status: code passed successfully
   - results: checked!

### 2.2 AlignSENSE
1. AlignSENSE recon with T1 MPRAGE (using ACS)
   - status: code passed successfully 
   - results: checked! 
2. AlignSENSE recon with SWI small FOV (using ACS)
   - status: code passed successfully
   - results: checked!
3. AlignSENSE: recon with SWI small FOV (using exterREF)
   - status: code passed successfully
   - results: checked!
