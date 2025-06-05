function modals=modalList(printL)

%MODALLIST   Returns the list of the modalities contemplated in the
%reconstruction.
%   MODALS=MODALLIST({PRINTL})
%   * {PRINTL} indicates whether to print the list with an explanation or
%   not (default)
%   * MODAL is the list of returned modalities
%

if nargin<1 || isempty(printL);printL=0;end

modals{1}='Re-Su';%Survey
modals{2}='Re-Se';%Coils
modals{3}='Re-F0';%BO
modals{4}='Re-F1';%B1
modals{5}='An-S2';%Multislice T2
modals{6}='An-S1';%Multislice T1
modals{7}='An-Ve';%Volumetric encoding (3D)
modals{8}='';%Not used
modals{9}='Dy-Fu';%Functional
modals{10}='Dy-Di';%Diffusion
modals{11}='Dy-Ca';%Cine cardiac
modals{12}='XX-GT';%This is not a modality but indicates ReconFrame reconstructions for testing purposes

if printL
    expln{1}='Survey data';
    expln{2}='Coil data';
    expln{3}='B0 data';
    expln{4}='B1 data';
    expln{5}='Multi-slice T2';
    expln{6}='Multi-slice T1';
    expln{7}='Volumetric encoding 3D';
    expln{8}='Undefined';%Not defined
    expln{9}='Functional';
    expln{10}='Diffusion';
    expln{11}='Cine cardiac';
    expln{12}='ReconFrame reconstruction';
            
    fprintf('List of modalities: \n');
    fprintf('Number Name Description\n');
    for n=1:length(modals);fprintf('%d %s %s\n',n,modals{n},expln{n});end;fprintf('\n');
end
