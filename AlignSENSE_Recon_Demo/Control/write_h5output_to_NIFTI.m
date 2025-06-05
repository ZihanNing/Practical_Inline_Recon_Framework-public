

function [] = write_h5output_to_NIFTI (connection)

    disph5Info=0;
    writeNII =1;
    fileName = 'DISORDEROut';
    ISMRMRD2NII (fileName, writeNII, disph5Info);
    
end