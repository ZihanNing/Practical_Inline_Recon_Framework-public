function tdfwrite(filename,st)

%TDFWRITE   Writes a structure into a .tab file. The code extends the
%function in https://uk.mathworks.com/matlabcentral/fileexchange/25622-tdfwrite--export-data-in-tab-delimited-format-?requestedDomain=www.mathworks.com
%from Rafael Palacios
%   TDFWRITE(FILENAME,ST)
%   * FILENAME is the file name
%   * ST is the structure to write. It is a structure created with
%   ST=TDFREAD('FILE.TAB'); It is a structure with several fields. Each 
%   field is a vector of numbers or a matrix of char. Field names are used 
%   as headers of each column
%

%ERROR CONTROL
assert(ischar(filename),'First argument must be the name of the file');%LCG modification
assert(isstruct(st),'Second argument must be a strcuture');%LCG modification

%FIELD NAMES
names=fieldnames(st);
rows=size(st.(names{1}),1);
for j=2:length(names);assert(rows==size(st.(names{j}),1),'Field $s has a different length than first field (%s)',names{j},names{1});end 

%STRUCTURE CREATION AND WRITING
[fp,message]=fopen(filename,'w');
assert(fp>=0,'Error opening file: %s',message);
fprintf(fp,'%s',names{1});fprintf(fp,'\t%s',names{2:end});fprintf(fp,'\n');
for i=1:rows
    for j=1:length(names)
        if j~=1;fprintf(fp,'\t');end
        v=st.(names{j});
        if iscell(v)
            fprintf(fp,'%s',v{i});
        else%LCG modification
            if (ischar(v(1,1)));fprintf(fp,'%s',v(i,:));
            %The instruction below will only work from MatlabR2016b
            %elseif (isstring(v(1,1)));fprintf(fp,'%s',v(i));
            elseif (isinteger(v(1,1)));fprintf(fp,'%d',v(i));
            else fprintf(fp,'%.20g',v(i));%General format
            end
        end
    end
    fprintf(fp,'\n');
end
fclose(fp);
