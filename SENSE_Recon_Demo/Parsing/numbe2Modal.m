function modal=numbe2Modal(numbe,ce)

%NUMBE2MODAL   Converts from an index for the modality to the modality as a
%string. Modalities are defined in the file modalList.m
%   MODAL=NUMBE2MODAL(NUMBE,CE)
%   * NUMBE is a vector containing the indexes of the modalities
%   * MODAL is a string / strings containing the modality / modalities
%   names
%   * {CE} indicates whether the return a single modality as a cell or not
%   (default)
%

if nargin<2 || isempty(ce);ce=0;end

modals=modalList;
N=length(modals);
numbesq=numbe;numbesq(numbesq<0)=[];numbesq(numbesq>N)=[];
if isempty(numbesq)
    fprintf('None of the requested modality numbers (');for m=1:length(numbe);fprintf('%s ',numbe(m));end;fprintf('could be found\n');    
    modalList(1);error('See the modality list and correct the numbers');
end
    
if any(numbe)>N || any(numbe)<0;fprintf('Non positive or too big numbers for modalities found. They have been squeezed\n');end
numbe=numbesq;

N=length(numbe);
if ce || N>1;modal=cell(1,N);end
for n=1:N
    if ce || n>1;modal{n}=modals{numbe(n)};else modal=modals{numbe(n)};end
end
