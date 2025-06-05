function e=clique_energy_ho(d,p,pot)
%clique_energy_ho  Computes clique energy: e=th^(p-2)*d.^2.*mask + d.^p.*(1-mask)
%        e=clique_energy_ho(d,p,th,quant)
%
%  Input arguments --------------------
%  d             -> clique difference
%  p             -> power law exponent
%  pot.thres            -> it defines a region over which the potential grows quadratically
%  pot.quant         -> it defines whether or not the potential is quantized

switch pot.quant
    case 'no';d=abs(d);% non quantized potential        
    case 'yes';d=abs(round(d/2/pi)*2*pi);% quantized potential (2pi Quantization of phase difference)        
end

if pot.thres~=0
    mask = (d<=pot.thres);
    e = pot.thres^(p-2)*d.^2.*mask + d.^p.*(1-mask);
else
    e = d.^p;
end


