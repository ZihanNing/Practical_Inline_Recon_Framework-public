function erg = energy_hoND(k,ext_psi,ROI,p,c,ext_w,pot)
%energy_ho   Energy from k labeling and psi phase measurements.
%   erg = energy_ho(k,psi,base,p,c,disc_bar,p,th,quant) returns the energy of kappa labeling given the 
%   psi measurements image, the base ROI image (having ones in the region of interest (psi) and a passe-partout
%   made of zeros), the exponent p, the cliques matrix (each row indicating a displacement vector corresponding
%   to each clique), the disc_bar (complement to one of the quality maps), a threshold pot.thres defining a region for
%   which the potential (before a possible quantization) is quadratic, and pot.quant which is a flag defining whether
%   the potential is or is not quantized.
%   (see J. Bioucas-Dias and G. Valadao, "Phase Unwrapping via Graph Cuts"
%   submitted to IEEE Transactions Image Processing, October, 2005).
%   SITE: www.lx.it.pt/~bioucas/ 

ND=numDims(ext_psi);
NC = size(c);
maxCl = max(abs(c(:)));

ext_k=padarray(k,(maxCl+1)*ones(1,ND));
NN=size(ext_w);NN(end+1:3)=1;

erg=0;
for t=1:NC(1)
    ROIPos=circshift(ROI,c(t,:)).*ROI;
    dk= ext_k-circshift(ext_k,c(t,:));          
    dpsi = circshift(ext_psi,c(t,:)) - ext_psi;   
    a= (2*pi*dk-dpsi).*ROIPos;%Bug in PUMA code: .*reshape(ext_w(t,:,:,:,:,:),NN(2:end))
    b=clique_energy_ho(a,p,pot).*reshape(dynInd(ext_w,t,1),NN(2:end));
    erg=erg+sum(b(:));
end


