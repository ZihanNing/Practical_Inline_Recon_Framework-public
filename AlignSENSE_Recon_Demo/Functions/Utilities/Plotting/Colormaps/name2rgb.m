% Function that converts an HEX string (format: #FFFFFF) to the
% corresponding RGB color
function rgb = name2rgb(name)

    if strcmp(name, 'yellow'); rgb = [1 1 0];
    elseif strcmp(name, 'magenta'); rgb = [1 0 1];
    elseif strcmp(name, 'cyan'); rgb = [0 1 1];
    elseif strcmp(name, 'red'); rgb = [1 0 0];
    elseif strcmp(name, 'green'); rgb = [0 1 0];
    elseif strcmp(name, 'blue'); rgb = [0 0 1];
    elseif strcmp(name, 'white'); rgb = [1 1 1];
    elseif strcmp(name, 'black'); rgb = [0 0 0];
    else; error('name2rgb:: %s not a color or not supported. %s', name);
    end
    
end
