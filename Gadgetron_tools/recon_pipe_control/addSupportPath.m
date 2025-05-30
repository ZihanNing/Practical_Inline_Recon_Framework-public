function supportPath = addSupportPath(Recon_ID)
% addSupportPath  Add support_func_path for a given Recon_ID to MATLAB path
%
%   supportPath = addSupportPath(Recon_ID) parses Framework_config.xml in the
%   current folder, locates <Component name="Recon_ID">, extracts its
%   <support_func_path> value, and adds that folder (and its subfolders) to
%   the MATLAB search path. Returns the full path string added.
% 
%   by Zihan Ning

    xmlFile = 'Framework_config.xml';
    if exist(xmlFile,'file')~=2
        error('addSupportPath:NoFile', 'Cannot find %s in current folder.', xmlFile);
    end

    % Parse XML
    try
        xDoc = xmlread(xmlFile);
    catch ME
        error('addSupportPath:ParseError', 'Failed to parse %s: %s', xmlFile, ME.message);
    end

    % Locate the <Component> node
    comps = xDoc.getElementsByTagName('Component');
    compNode = [];
    for k = 0:comps.getLength-1
        node = comps.item(k);
        if strcmp(char(node.getAttribute('name')), Recon_ID)
            compNode = node;
            break;
        end
    end
    if isempty(compNode)
        error('addSupportPath:NoComponent', ...
              'No <Component name="%s"> found in %s.', Recon_ID, xmlFile);
    end

    % Extract <support_func_path>
    nodes = compNode.getElementsByTagName('support_func_path');
    if nodes.getLength < 1
        error('addSupportPath:MissingTag', ...
              '<support_func_path> not found under Component "%s".', Recon_ID);
    end
    textNode = nodes.item(0).getFirstChild;
    if isempty(textNode) || all(isspace(char(textNode.getData)))
        error('addSupportPath:EmptyPath', ...
              '<support_func_path> for Component "%s" is empty.', Recon_ID);
    end
    supportPath = strtrim(char(textNode.getData));

    % Verify directory exists
    if ~isfolder(supportPath)
        error('addSupportPath:InvalidPath', ...
              'support_func_path "%s" does not exist or is not a folder.', supportPath);
    end

    % Add to MATLAB path
    addpath(genpath(supportPath));
    fprintf('Added support path: %s (and subfolders)\n', supportPath);
end
