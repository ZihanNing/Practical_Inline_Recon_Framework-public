function [kloc,dcf] = projection(Nproj,NperProj,ndim)
% Generate trajectory and density compensation function for projection trajectory.
%   
% (c) Corey Baron 2015
%

switch ndim
    case 2
        phi0 = linspace(0,pi,Nproj+1);
        phi = repmat(phi0(1:end-1),[NperProj, 1]);
        r0 = linspace(-0.5,0.5,NperProj+1);
        r0 = r0(1:end-1) + (r0(2)-r0(1))/2;
        r = repmat(r0', [1 Nproj]);
        kloc = cat(2,r(:).*cos(phi(:)), r(:).*sin(phi(:))); 
        dcf = abs(r(:)*pi/Nproj);
        dcf = dcf/max(dcf(:));
    otherwise 
        error('requested dim not coded!')
end

end

