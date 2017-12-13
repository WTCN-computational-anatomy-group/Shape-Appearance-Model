function WWv = UpdateWWv(Wv,s)

Kv  = size(Wv,5);
WWv = zeros(Kv);
for k=1:Kv,
    spm_diffeo('boundary',0);
    u0 = spm_diffeo('vel2mom', Wv(:,:,:,:,k), [s.vx s.v_settings]);
    for k1=k:Kv,
        WWv(k,k1) = sum(sum(sum(sum(u0.*Wv(:,:,:,:,k1)))));
        WWv(k1,k) = WWv(k,k1);
    end
end

