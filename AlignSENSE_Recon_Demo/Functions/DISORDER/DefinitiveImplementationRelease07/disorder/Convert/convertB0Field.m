function x=convertB0Field(x,TE,ES,inp,out,N)

%CONVERTFIELD   Converts units of B0 field
%   X=CONVERTFIELD(X,TE,ES,INP,OUT,{N})
%   * X is a B0 field in units given by INP
%   * TE is the echo time at which a given field-induced dephasing has been observed %YB: In ms
%   * ES is the inter-echo spacing which scales the distortion produced by the field %YB: In ms
%   * INP are the input units of the field
%   * OUT are the output units of the field
%   * {N} is the grid size in the phase encoding direction
%   ** X is a B0 field in units given by OUT
%

if nargin<6;N=size(x,2);end

if strcmp(inp,'Hz')
    if strcmp(out,'pix');x=x*N*ES/1000;elseif strcmp(out,'rad');x=bsxfun(@times,x,2*pi*TE/1000);elseif strcmp(out,'res');x=x*ES/1000;end
elseif strcmp(inp,'pix')
    if strcmp(out,'Hz');x=x*1000/(ES*N);elseif strcmp(out,'rad');x=bsxfun(@times,x,2*pi*TE/(ES*N));elseif strcmp(out,'res');x=x/N;end
elseif strcmp(inp,'rad')
    if strcmp(out,'Hz');x=bsxfun(@times,x,1000./(2*pi*TE));elseif strcmp(out,'pix');x=bsxfun(@times,x,N*ES./(2*pi*TE));elseif strcmp(out,'res');x=bsxfun(@times,x,ES./(2*pi*TE));end
elseif strcmp(inp,'res')
    if strcmp(out,'Hz');x=x*1000/ES;elseif strcmp(out,'pix');x=x*N;elseif strcmp(out,'rad');x=bsxfun(@times,x,(2*pi*TE)./ES);end
end