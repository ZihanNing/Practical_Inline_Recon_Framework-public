function supportPaths = addSupportPath(Recon_ID)
% addSupportPath  Add one or more support_func_path entries for a given Recon_ID to MATLAB path
%
%   supportPaths = addSupportPath(Recon_ID) parses Framework_config.xml in the
%   current folder, locates <Component name="Recon_ID">, extracts all
%   <support_func_path> values, verifies they exist, and adds each (and its
%   subfolders) to the MATLAB search path. Returns a cell array of the full
%   paths added.
%
%   Example XML fragment:
%     <Component name="NUFFT_radial">
%         ...
%         <support_func_path>/home/zn23/matlab/merina_recon</support_func_path>
%         <support_func_path>/home/zn23/matlab/another_toolbox</support_func_path>
%     </Component>
%
%   Then calling addSupportPath('NUFFT_radial') will add both folders (and
%   subfolders) to the MATLAB path and return:
%     {'/home/zn23/matlab/merina_recon', '/home/zn23/matlab/another_toolbox'}
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

    % Extract all <support_func_path> nodes
    nodes = compNode.getElementsByTagName('support_func_path');
    nNodes = nodes.getLength;
    if nNodes < 1
        error('addSupportPath:MissingTag', ...
              'No <support_func_path> tags found under Component "%s".', Recon_ID);
    end

    % Preallocate cell array for paths
    supportPaths = cell(1, nNodes);

    for idx = 0:(nNodes-1)
        textNode = nodes.item(idx).getFirstChild;
        if isempty(textNode) || all(isspace(char(textNode.getData)))
            error('addSupportPath:EmptyPath', ...
                  '<support_func_path> entry #%d under Component "%s" is empty.', idx+1, Recon_ID);
        end
        pathStr = strtrim(char(textNode.getData));
        if ~isfolder(pathStr)
            error('addSupportPath:InvalidPath', ...
                  'support_func_path "%s" does not exist or is not a directory.', pathStr);
        end
        % Add this folder (and subfolders) to MATLAB path
        addpath(genpath(pathStr));
        supportPaths{idx+1} = pathStr;
        fprintf('Added support path: %s (and subfolders)\n', pathStr);
    end
end
