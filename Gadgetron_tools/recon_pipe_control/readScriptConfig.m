function cfg = readScriptConfig(Recon_ID)
% readScriptConfig  Read GPU‐wrapper params (and optional input‐check flags) from Framework_config.xml
%   cfg = readScriptConfig(Recon_ID) parses Framework_config.xml in the
%   current folder, locates <Component name="Recon_ID">, and returns a struct
%   with fields:
%     .logFile          (string)
%     .maxWaitTime      (double)
%     .interval         (double)
%     .matlabFunction   (string)
%     .flagCheckInputs  (0 or 1; default = 0 if missing)
%     .inputList        (cell array of string patterns; empty if flagCheckInputs==0)
%
% Example XML structure expected under Component:
%   <Component name="SENSE_exterREF">
%     <logFile>...</logFile>
%     <maxWaitTime>900</maxWaitTime>
%     <interval>10</interval>
%     <matlabFunction>handle_connection_background_sense_acs</matlabFunction>
%     <flag_checkInputs>1</flag_checkInputs>
%     <InputList>
%       <Input>ExterREF_CMS_*</Input>
%       <Input>AnotherPattern_*.mat</Input>
%     </InputList>
%   </Component>

    xmlFile = 'Framework_config.xml';
    if exist(xmlFile,'file')~=2
        error('readScriptConfig:NoFile', 'Cannot find %s in current folder.', xmlFile);
    end

    try
        xDoc = xmlread(xmlFile);
    catch ME
        error('readScriptConfig:ParseError', 'Failed to parse %s: %s', xmlFile, ME.message);
    end

    % Find the matching <Component> node by name
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

    % Helper: get text content of a single-child tag; throws if missing
    function txt = getTag(tagName)
        nodes = compNode.getElementsByTagName(tagName);
        if nodes.getLength < 1
            error('readScriptConfig:MissingTag', ...
                  '<%s> not found under Component "%s".', tagName, Recon_ID);
        end
        txt = strtrim(char(nodes.item(0).getFirstChild.getData));
    end

    % 1) Required fields
    cfg.logFile        = getTag('logFile');
    cfg.maxWaitTime    = str2double(getTag('maxWaitTime'));
    cfg.interval       = str2double(getTag('interval'));
    cfg.matlabFunction = getTag('matlabFunction');

    if isnan(cfg.maxWaitTime) || isnan(cfg.interval)
        error('readScriptConfig:BadNumber', ...
              'maxWaitTime or interval for Component "%s" is not numeric.', Recon_ID);
    end

    % 2) Optional: flag_checkInputs (default = 0 if missing)
    try
        rawFlag = getTag('flag_checkInputs');
        cfg.flagCheckInputs = double(str2double(rawFlag) == 1);
    catch
        % tag not present → assume 0
        cfg.flagCheckInputs = 0;
    end

    % 3) If flagCheckInputs == 1, collect all <Input> under <InputList>
    cfg.inputList = {};  % default: empty cell
    if cfg.flagCheckInputs == 1
        % Find the <InputList> node
        inputLists = compNode.getElementsByTagName('InputList');
        if inputLists.getLength < 1
            error('readScriptConfig:MissingTag', ...
                  '<InputList> not found under Component "%s", but flag_checkInputs==1.', Recon_ID);
        end
        inputListNode = inputLists.item(0);
        % Now find every <Input> child
        inputs = inputListNode.getElementsByTagName('Input');
        nInputs = inputs.getLength;
        if nInputs < 1
            error('readScriptConfig:NoInputs', ...
                 'No <Input> entries under <InputList> for Component "%s".', Recon_ID);
        end
        temp = cell(1, nInputs);
        for ii = 0:(nInputs-1)
            node = inputs.item(ii);
            % text content of <Input>
            txt = strtrim(char(node.getFirstChild.getData));
            temp{ii+1} = txt;
        end
        cfg.inputList = temp;
    end
end
