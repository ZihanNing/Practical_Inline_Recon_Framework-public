function [value,status] = searchUPfield(hdr_UP,type,name,debug)
%%% This is the function used to find the user defined parameters in the
%%% connection_hdr
%%% To note: you need to define those parameters in the ParameterMap when
%%% convert twix to ismrmrd (for detail: see siemens_to_ismrmrd git repo)
%%%
%%% INPUT:
%%% hdr_UP - the field contains all the user defined headers in current
%%% connection header
%%% type - the type of the user-defined header (double, long, string)
%%% name - the name of the user-defined header (should be pre-defined in
%%% parameterMap)
%%% debug - display debug log
%%%
%%% OUTPUT:
%%% value: the value of the user-defined header
%%% status: found the parameter or the parameter is not exist
%%%
%%% by Zihan Ning @ King's
%%% last updated time: 4-Mar-2025

if nargin < 4 || isempty(debug); debug=0;end

% hdr_UP.Double = connection_hdr.userParameters.userParameterDouble;
% hdrUP.Long = connection_hdr.userParameters.userParameterLong;
% hdrUP.String = connection_hdr.userParameters.userParameterString;

switch type
    case 'Double'
        if isfield(hdr_UP,'Double');hdr = hdr_UP.Double;end
        if isfield(hdr_UP,'userParameterDouble');hdr = hdr_UP.userParameterDouble;end
    case 'Long'
        if isfield(hdr_UP,'Long');hdr = hdr_UP.Long;end
        if isfield(hdr_UP,'userParameterLong');hdr = hdr_UP.userParameterLong;end
    case 'String'
        if isfield(hdr_UP,'String');hdr = hdr_UP.String;end
        if isfield(hdr_UP,'userParameterString');hdr = hdr_UP.userParameterString;end
    otherwise
        disp('Unrecognizable type!!Abort!');
        status = 'aborted';
        return
end

status = [];
if ~isempty(hdr)
    for i = 1:length(hdr)
        if strcmp(hdr(i).name,name)
            value = hdr(i).value;
            status = 'found';
            if debug; fprintf('Found %s in the user-defined %s field.\n',name, type);end
        end
    end
else
    disp('Empty user-defined parameter field!! Abort!');
    status = 'aborted';
    return
end
if isempty(status) % no matched user-defined header
    value = []; status = 'not-found';
    if debug; fprintf('Cannot found %s in the user-defined %s field.\n',name, type);end
end

