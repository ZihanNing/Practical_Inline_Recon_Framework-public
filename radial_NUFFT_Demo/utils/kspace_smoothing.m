function smoothed_volume = kspace_smoothing(volume)

kspace = fftn(volume);
s = size(kspace);
for i=1:s(1)
    for j=1:s(2)
        for k=1:s(3)
            if ((abs(32-i) > 6) || (abs(32-j) > 6) || (abs(32-k) > 6))
                kspace(i,j,k) = 0 + 0j;
            end
        end
    end
end

smoothed_volume = ifftn(kspace); 

end

