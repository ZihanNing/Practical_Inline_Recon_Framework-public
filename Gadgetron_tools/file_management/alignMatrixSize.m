function [image,diff,status] = alignMatrixSize(hdr,image,log)
%%% Align the matrix size of the reconstructed image with the expected
%%% recon matrix size of the sequence (otherwise, the reconstructed image
%%% will be rejected when it sent back to the scanner)
%%%
%%% INPUT:
%%% hdr: in twix-like structure, use twix.hdr.encoding; or directly use
%%% connection.header.encoding
%%% image: reconstructed image with the dimension order [Col,Lin,Par,...]
%%% log: 1-give the log
%%% 
%%% OUTPUT:
%%% image: image size which align with the expected recon matrix size of
%%% the sequence (in twix.hdr.encoding)
%%% diff: differences between the matrix size of reconstructed image & the
%%% expected matrix size of the sequence
%%% status: No_need_for_processing or cropped or interpolated
%%%
%%% by Zihan @ king's
%%% March-2025

 
Seq_matrixSize = [hdr.encoding.reconSpace.matrixSize.x,...
    hdr.encoding.reconSpace.matrixSize.y,...
    hdr.encoding.reconSpace.matrixSize.z];% ZN: the expected recon matrix size of the sequence (to received the reconstructed image by gadgetron)
Acq_matrixSize = size(image,[1:3]); % ZN: the matrix size of the image to be sent [Col,Lin,Par,...]
status = 'No_need_for_processing';
diff = zeros(size(Seq_matrixSize));
if ~isequal(Seq_matrixSize,Acq_matrixSize)% ZN: check if the matrix size of the seq & acq is identical (may cause send failure if not)
    diff = Seq_matrixSize - Acq_matrixSize;
    if log; fprintf('The matrix size of the reconstructed image is not identical to the it of the sequence (seq - acq = [%s]).\n',...
            num2str(diff)); end
    
    % CROP IF SEQ IMAGE IS SMALLER
    if Seq_matrixSize(1)<= Acq_matrixSize(1) && Seq_matrixSize(2)<= Acq_matrixSize(2) && Seq_matrixSize(3)<= Acq_matrixSize(3)
        % crop the region in the middle
        start_point = [floor((- Seq_matrixSize(1) + Acq_matrixSize(1))/2)+1,...
            floor((- Seq_matrixSize(2) + Acq_matrixSize(2))/2)+1,...
            floor((- Seq_matrixSize(3) + Acq_matrixSize(3))/2)+1];
        end_point = [start_point(1)+Seq_matrixSize(1)-1,start_point(2)+Seq_matrixSize(2)-1,start_point(3)+Seq_matrixSize(3)-1];
        % crop
        dims = ndims(image); % [Col,Lin,Par,...]
        for d = 1:dims;if (d<=3);idx{d} = start_point(d):end_point(d); else; idx{d} = 1:size(image,d);end; end
        image = image(idx{:});
        if log; fprintf('The image has been cropped in the middle as the size as the sequence matrix size [%s]!\n',...
            num2str(size(image,[1:3]))); end
        status = 'Cropped';
    else % INTERPORATE IF SEQ IMAGE IS LARGER
        fprintf('!!!!!!!!! WARNING !!!!!!!!!!!\n');
        fprintf('The reconstructed image is smaller than the expected matrix size of the sequence at least in one dimension!\n');
        fprintf('The reconstructed image is forced to be the same size as the sequence!\n');
        % interpolate slab by slab
        ori_size = size(image); % [Col,Lin,Par,...]
        squeezed_dim = prod(ori_size(4:end));
        image_squeezed = reshape(image,[ori_size(1:3),squeezed_dim]);
        for i = 1:squeezed_dim
            image_squeezed(:,:,:,i) = imresize(image_squeezed(:,:,:,i),Seq_matrixSize,'bilinear'); % linear interpolation
        end
        image = reshape(image_squeezed,ori_size);
        fprintf('The image has been interpolated (bilinear) to the size of the sequence matrix size [%s]!\n',...
            num2str(size(image,[1:3])));
        status = 'Interpolated';
    end
else
    if log; fprintf('[Success] The matrix size of the reconstructed image is identical to the expected matrix size of the sequence now.\n'); end
end