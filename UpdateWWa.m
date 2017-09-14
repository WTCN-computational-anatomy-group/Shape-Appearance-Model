function WWa = UpdateWWa(Wa,s)

Ka  = size(Wa,5);
WWa = zeros(Ka);
for k=1:Ka,
    for l=1:size(Wa,4)
        a0 = spm_field('vel2mom', Wa(:,:,:,l,k), [s.vx s.a_settings]);
        for k1=k:Ka,
            WWa(k,k1) = WWa(k,k1) + sum(sum(sum(a0.*Wa(:,:,:,l,k1))));
            WWa(k1,k) = WWa(k,k1);
        end
    end
end

