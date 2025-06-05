function [fileToSend_all,status] = save_for_retrieval(x,twix,config,flag_savewometa)
% This is a function to save the file/image to be retrieved by a retrieval
% dummy scan or a retro-recon
% Input:
%   x - image to be send [RO Lin Par Con]
%   twix - twix_like raw, to generate the essential header for retrieval  
%   image
%   config - configurations for saving
%       send_mag
%       send_phs
%       save_path
%       fileName
%       seq_type
%
% Output:
%   fileToSend_all: file to be retrieved (image & headers)
%   status: 0-fail; 1-success
%
% by Zihan Ning @ King's college london
% 2 June 2025

% flag_savewometa = 1, just save the reconstructed image without its
% ISMRMRD metadata generated from the raw
% This will be essential for the sequence that the scanner cannot provide a
% conventional reconstruction (e.g., non-cartesian sequences)
% the metadata of the image will be generated from the raw of the retrieval
% dummy scan, which should be aligned for the matrix size (and FOV) with
% the reconstructed image
if nargin < 4 || isempty(flag_savewometa); flag_savewometa=0;end 

status = 0;
%% POST-PROCESSING OF THE IMAGE TO BE SEND OUT
if length(size(x))>3; NCon = size(x,4); else; NCon = 1; end % if multi-inv or multi-echo
for con = 1:NCon
    x(:,:,:,con) = flipPermute(x(:,:,:,con), [1 1 1]);
    x(:,:,:,con) = circshift(x(:,:,:,con), [1 1 1]);
end
% ZN: shift the image due to slice shift in geom (to be fixed)
x = circshift(x,-1,3);
% ZN: make the matrix size of the reconstructed image align with
% the sequence to receive the image from gadgetron 
if ~isequal(config.seq_type,'radial'); [x,~,~] = alignMatrixSize(twix.hdr,x,1); end % might need to be modified for non-Cartesian traj


%% CREATE FILE TO SEND
ref_func = @referenceFromReconData_general;

% GENERATE THE FILE TO BE SENT
fileToSend_all = {};
if config.send_mag; fileToSend_all.Mag = {}; end
if config.send_phs; fileToSend_all.Phs = {}; end
% if more types of outputs, add them here...

for con = 1:NCon
    if config.send_mag % Magnitude
        if flag_savewometa % save the image without meta
            fileToSend.data = permute(abs(x(:,:,:,con)), [4 1:3]); %Coil profiles first by convention of this (gadgetron) function
        else % generate headers and meta with the image
            fileToSend = ...
                gadgetron.types.Image.from_data(...
                                                permute(abs(x(:,:,:,con)), [4 1:3]), ...%Coil profiles first by convention of this (gadgetron) function
                                                ref_func(twix) ... % generate the header of the reconstructed image based on twix-like raw
                                                );
        end
        fileToSend.header.image_type = gadgetron.types.Image.MAGNITUDE;
        fileToSend.header.data_type = gadgetron.types.Image.FLOAT;
        fileToSend.header.image_series_index = 1; % ZN: to be within a image tag (MAG)
        fileToSend.header.image_index = con; % ZN: to have recognizable sub-tag (contrasts/invs/echoes...)
        fileToSend_all.Mag{end+1} = fileToSend; % add to the structure
    end
    if config.send_phs % Phase
        if flag_savewometa % save the image without meta
            fileToSend.data = permute(angle(x(:,:,:,con)), [4 1:3]); %Coil profiles first by convention of this (gadgetron) function
        else % generate headers and meta with the image 
        fileToSend = ...
            gadgetron.types.Image.from_data(...
                                            permute(angle(x(:,:,:,con)), [4 1:3]), ...%Coil profiles first by convention of this (gadgetron) function
                                            ref_func(twix) ... % generate the header of the reconstructed image based on twix-like raw
                                            );
        end
        fileToSend.header.image_type = gadgetron.types.Image.PHASE;
        fileToSend.header.data_type = gadgetron.types.Image.FLOAT;
        fileToSend.header.image_series_index = 2; % ZN: to be within a image tag (PHS)
        fileToSend.header.image_index = con; % ZN: to have recognizable sub-tag (contrasts/invs/echoes...)
        fileToSend_all.Phs{end+1} = fileToSend; % add to the structure
    end
    % ZN: if more types of ouputs, add them here...
    % recommend you to add them with different image_series_index
end

% SAVE THE FILE (WITH CORRECT FORMAT)
% save the fileToSend for now and close the current program;
% waiting for another dummy scan to trigger the retrieving
save_path = config.save_path;
% ZN: this version is saved for backup, but not read by auto dummy scan (can be used for callback dummy by manually selecting this file)
image_name = ['image_to_be_send_',config.fileName,'.mat'];
nameRef = fullfile(save_path, image_name); 
save(nameRef, 'fileToSend_all','-v7.3');
% ZN: this version is saved for auto dummy scan, but will be overwritten and save the last one with the same seq type
% to note, the dummy scan should be as the same seq_type for a right pick
nameRef = fullfile(save_path, ['image_to_be_send_',config.seq_type,'.mat']); 
save(nameRef, 'fileToSend_all','-v7.3');
status = 1; 