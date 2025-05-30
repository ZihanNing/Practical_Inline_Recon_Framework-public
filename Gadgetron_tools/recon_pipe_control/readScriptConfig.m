function cfg = readScriptConfig(Recon_ID)
% readScriptConfig  Read GPU‐wrapper params for a given Recon_ID
%   cfg = readScriptConfig(Recon_ID) parses Framework_config.xml in the
%   current folder, locates <Component name="Recon_ID">, and returns a struct
%   with fields: logFile, maxWaitTime, interval, matlabFunction.

    xmlFile = 'Framework_config.xml';
    if exist(xmlFile,'file')~=2
        error('readScriptConfig:NoFile', 'Cannot find %s in current folder.', xmlFile);
    end

    try
        xDoc = xmlread(xmlFile);
    catch ME
        error('readScriptConfig:ParseError', 'Failed to parse %s: %s', xmlFile, ME.message);
    end

    % Find the matching <Component> node
    comps = xDoc.getElementsByTagName('Component');
    compNode = [];
    for k = 0:comps.getLength-1
        thisNode = comps.item(k);
        if strcmp(char(thisNode.getAttribute('name')), Recon_ID)
            compNode = thisNode;
            break;
        end
    end
    if isempty(compNode)
        error('readScriptConfig:NoComponent', ...
              'No <Component name="%s"> found in %s.', Recon_ID, xmlFile);
    end

    % Helper to get text from a child tag of compNode
    function txt = getTag(tagName)
        nodes = compNode.getElementsByTagName(tagName);
        if nodes.getLength < 1
            error('readScriptConfig:MissingTag', ...
                  '<%s> not found under Component "%s".', tagName, Recon_ID);
        end
        txt = strtrim(char(nodes.item(0).getFirstChild.getData));
    end

    % Populate output struct
    cfg.logFile        = getTag('logFile');
    cfg.maxWaitTime    = str2double(getTag('maxWaitTime'));
    cfg.interval       = str2double(getTag('interval'));
    cfg.matlabFunction = getTag('matlabFunction');

    % Validate numeric fields
    if isnan(cfg.maxWaitTime) || isnan(cfg.interval)
        error('readScriptConfig:BadNumber', ...
              'maxWaitTime or interval for Component "%s" is not numeric.', Recon_ID);
    end
end
