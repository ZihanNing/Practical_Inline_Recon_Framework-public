% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 22-May-2025
%
% getReconAlg  Return the recon_alg for a given Recon_ID in Framework_config.xml
%   recon_alg = getReconAlg(Recon_ID) looks for Framework_config.xml on the
%   MATLAB path, parses it, finds the <Component> whose name attribute matches
%   Recon_ID, and returns the contents of its <recon_alg> element.
%
%   Errors if the XML file is not found, if no matching component exists,
%   or if the recon_alg element is missing or empty.
function recon_alg = getReconAlg(Recon_ID)

    xmlFile = 'Framework_config.xml';

    % 1. Check existence on path or in current folder
    if exist(xmlFile, 'file') ~= 2
        error('getReconAlg:NoConfigFile', ...
              'Cannot find %s on MATLAB path or current folder.', xmlFile);
    end

    % 2. Parse XML
    try
        xDoc = xmlread(xmlFile);
    catch ME
        error('getReconAlg:ParseError', ...
              'Failed to parse %s: %s', xmlFile, ME.message);
    end

    % 3. Find all <Component> nodes
    comps = xDoc.getElementsByTagName('Component');
    recon_alg = '';
    for k = 0:comps.getLength-1
        comp = comps.item(k);
        nameAttr = char(comp.getAttribute('name'));
        if strcmp(nameAttr, Recon_ID)
            % Found the component; look for recon_alg child
            algNodes = comp.getElementsByTagName('recon_alg');
            if algNodes.getLength < 1
                error('getReconAlg:NoReconAlg', ...
                      'Component "%s" has no <recon_alg> defined.', Recon_ID);
            end
            textNode = algNodes.item(0).getFirstChild;
            if isempty(textNode) || all(isspace(char(textNode.getData)))
                error('getReconAlg:EmptyReconAlg', ...
                      '<recon_alg> for component "%s" is empty.', Recon_ID);
            end
            recon_alg = strtrim(char(textNode.getData));
            return
        end
    end

    % If we reach here, no matching component was found
    error('getReconAlg:UnknownReconID', ...
          'No <Component name="%s"> found in %s.', Recon_ID, xmlFile);
end
