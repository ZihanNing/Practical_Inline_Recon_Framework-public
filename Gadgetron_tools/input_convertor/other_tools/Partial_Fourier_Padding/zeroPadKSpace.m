function [data,kspace_Pad, Lin, Par] = zeroPadKSpace(data,matrix_size,Acq_matrix_size,Lin_center,Par_center,Lin, Par, display)
%%% This function handling the zero-padding in kspace in terms of partial
%%% Fourier or (recon_resol < acq_resol)
%%%
%%% INPUT:
%%% data: kspace data [Col, Cha, Lin, Par, Ave, Con]
%%% matrix_size: expected matrix size after padding
%%% Acq_matrix_size: original matrix size before padding
%%% Lin_center: center point of Lin
%%% Par_center: center point of Par
%%% Lin: sampled idx 
%%% Par: sampled idx
%%% Display: to log
%%% 
%%% OUTPUT:
%%% data - padded kspace data
%%% kspace_Pad: status of zero-padding
%%% Lin & Par: accordingly modifed samplind idx
%%%
%%% by Zihan @ king's
%%% March-2025
    
    if nargin < 5 || isempty(display); display=1;end 

    % Calculate padding
    kspace_Pad = [0 0; 0 0]; % ZN: [Pre-pad Lin, Pre-pad Par; Post-pad Lin, Post-pad Par]
    kspace_Pad(1,1) = centerIdx(matrix_size(2)) - Lin_center; % Pre-padding in Lin 
    kspace_Pad(1,2) = centerIdx(matrix_size(3)) - Par_center; % Pre-padding in Par
    kspace_Pad(2,1) = matrix_size(2) - Acq_matrix_size(2) - kspace_Pad(1,1);
    kspace_Pad(2,2) = matrix_size(3) - Acq_matrix_size(3) - kspace_Pad(1,2);

    % Pad
    % only for Lin & Par, suppose that asymmetric echo is handled by AsymmetricEcho gadget
    fprintf('Zero-padding the kspace from [%s %s %s] to [%s %s %s]. \n',...
        num2str(Acq_matrix_size(1)),num2str(Acq_matrix_size(2)),num2str(Acq_matrix_size(3)),...
        num2str(Acq_matrix_size(1)),num2str(matrix_size(2)),num2str(matrix_size(3))); % ZN: here we do not pad RO direction
    if multDimSum(kspace_Pad(1,:))>0; data = padArrayND(data, [0 0 kspace_Pad(1,1) kspace_Pad(1,2)],[],0,'pre'); end % Pre-padding
    if multDimSum(kspace_Pad(2,:))>0; data = padArrayND(data, [0 0 kspace_Pad(2,1) kspace_Pad(2,2)],[],0,'post'); end % Post-padding
    
    % handling Lin & Par
    Lin = Lin + kspace_Pad(1,1); 
    Par = Par + kspace_Pad(1,2);

end