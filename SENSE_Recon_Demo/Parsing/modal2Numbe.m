function numbe=modal2Numbe(modal)

%MODAL2NUMBE   Converts from a modality as an string to an index for the
%modality. Modalities are defined in the file modalList.m
%   NUMBE=MODAL2NUMBE(MODAL)
%   * MODAL is a string / strings containing the modality / modalities
%   names
%   * NUMBE is a vector containing the indexes of the modalities
%

modals=modalList;
M=length(modals);
if isempty(modal);numbe=1:M;return;end
if ~iscell(modal);modalC{1}=modal;else modalC=modal;end

numbe=zeros(1,length(modalC));
for m=1:M
    for n=1:length(modals)
        if strcmp(modalC{m},modals{n})
            numbe(m)=n;
            break;
        end
    end
    if numbe(m)==0;fprintf('The requested modality %s could not be found so the resulting vector will be squeezed\n',modalC{m});end
end
numbe(numbe==0)=[];

if isempty(numbe)
    fprintf('None of the requested modalities (');for m=1:length(modalC);fprintf('%s ',modalC{m});end;fprintf('could be found\n');
    modalList(1);error('See the modality list and correct the name');
end