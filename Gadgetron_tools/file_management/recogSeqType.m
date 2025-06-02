function [seq_type,status] = recogSeqType(sequence_type,tScanningSequence,log)
%%% This function is to regonize the sequence type based on
%%% connection.header provided by gadgetron connection
%%%
%%% INPUT:
%%% sequence_type: usually at connection.header.sequenceParameters.sequence_type
%%% tScanningSequence: this is a customer addded parameter, you need to
%%% modify the stylesheet to achieve this: connection.header.userParameters.userParameterString.name
%%% 
%%% OUTPUT:
%%% seq_type: currently work for - MEGE, MPRAGE (MP2RAGE) sequence
%%% more sequence type to be added
%%%
%%% by Zihan @ king's
%%% March-2025

if nargin < 3 || isempty(log); log=0;end 

seq_type = [];
switch sequence_type
    case 'Flash'
        if isequal(tScanningSequence,'GR\IR')
            seq_type = 'mprage'; % for the cases of mprage or mp2rage, we consider they are the same
            status = 'Found';
        else
            seq_type = 'mege'; % multi echo gre
            status = 'Found';
        end
    case 'TurboSpinEcho' % TSE sequence
        if isequal(tScanningSequence,'SE\IR')
            seq_type = 'flair';
            status = 'Found';
        else
            seq_type = 'tse';
            status = 'Found';
        end
    otherwise % Other cases to be added...
        status = 'Failed';
end

if log
    if isequal(status,'Failed')
        fprintf('Cannot recognize such sequence type yet - might influence the following file management and reconstruction!\n');
    else
        fprintf('The sequence type for current data is recognized as %s.\n',seq_type);
    end
end
