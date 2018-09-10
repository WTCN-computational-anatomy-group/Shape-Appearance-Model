function WWa = UpdateWWa(Wa,s)
% Compute WWa from appearance basis functions
% FORMAT WWa = UpdateWWa(Wa,s)
%
% Wa  - Shape bases
% s   - Settings. Uses s.vx and s.a_settings.
%
% WWa - Wa'*La*Wa.
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

Ka  = size(Wa,5);
WWa = zeros(Ka);
for k=1:Ka
    for l=1:size(Wa,4)
        spm_field('boundary',0);
        a0 = spm_field('vel2mom', Wa(:,:,:,l,k), [s.vx s.a_settings]);
        for k1=k:Ka
            WWa(k,k1) = WWa(k,k1) + sum(sum(sum(a0.*Wa(:,:,:,l,k1))));
            WWa(k1,k) = WWa(k,k1);
        end
    end
end

