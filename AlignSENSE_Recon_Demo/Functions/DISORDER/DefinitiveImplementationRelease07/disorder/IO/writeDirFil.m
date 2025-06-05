function writeDirFil(filen,a)

%WRITEDIRFIL   Writes a DWI b-val/PEs file
%   WRITEDIRFIL(FILEN,A)
%   * FILEN is a filename
%   * A is the variable to be written
%

dlmwrite(filen,a,'delimiter',' ','precision',6);