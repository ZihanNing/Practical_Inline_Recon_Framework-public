
function [pathIn, refIn, B0In, B1In, fileIn] = extract_studies(pathIn, refIn,B0In, B1In, fileIn, id_path, id_ref,id_B0,id_B1, id_file )

assert(length(id_file)==length(id_ref),'id_ref and id_file should have the same size.')

if ~isempty(id_path)
    id_path(id_path>length(pathIn))=[];
    pathIn=pathIn(id_path);
    fileIn=fileIn(id_path);
    refIn=refIn(id_path);
    B0In=B0In(id_path);
    B1In=B1In(id_path);
end
if ~iscell(id_file) && ~isempty(id_file)
    for p=1:length(fileIn)
        id_file(id_file>length(fileIn{p}))=[];        
        fileIn{p}=fileIn{p}(id_file);
    end
end
if ~iscell(id_ref) && ~isempty(id_ref)
    for p=1:length(refIn)
        id_ref(id_ref>length(refIn{p}))=[];        
        refIn{p}=refIn{p}(id_ref);
    end
end
if ~iscell(id_B0) && ~isempty(id_B0)
    for p=1:length(B0In)
        id_B0(id_B0>length(B0In{p}))=[];        
        B0In{p}=B0In{p}(id_B0);
    end
end
if ~iscell(id_B1) && ~isempty(id_B1)
    for p=1:length(B1In)
        id_B1(id_B1>length(B1In{p}))=[];        
        B1In{p}=B1In{p}(id_B1);
    end
end
