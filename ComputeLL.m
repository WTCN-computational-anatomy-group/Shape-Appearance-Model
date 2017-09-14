function ll = ComputeLL(f,iphi,a0,s,noise)
% Compute log-likelihood
% FORMAT ll = ComputeLL(f,iphi,a0,s,noise)
% f     - This observation
% iphi  - Deformation field
% a0    - Appearance part of the model
% s     - Settings.
% noise - Noise information
%
% ll    - log-likelihood
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if true
    ll = ComputeLL_native(f,iphi,a0,s,noise);
else
    ll = ComputeLL_warped(f,iphi,a0,s,noise);
end


%==========================================================================
%
%==========================================================================
function ll = ComputeLL_warped(f,iphi,a0,s,noise)
% Compute ll of warped image data - less accurate
[f1,rho] = Push(f,iphi);
msk      = isfinite(f1) & isfinite(a0);
switch lower(s.likelihood)
case {'normal','gaussian'}
    D  = sum(rho(:));
    ll = 0;
    for l=1:d(4)
        a0l = a0(:,:,:,l);
        fl  = f1(:,:,:,l);
        ll  =  ll + noise.nu_factor*0.5*(D*(log(noise.lam(l)) - log(2*pi)) - noise.lam(l)*sum(sum(sum((fl(msk)-a0l(msk).*rho(msk)).^2./(rho+eps)))));
    end

case {'laplace'}
    b  = reshape(sqrt(1./(2*noise.lam)),[1,1,1,d(4)]);
    D  = sum(rho(:));
    r  = bsxfun(@times,a0-bsxfun(@rdivide,f1,rho+eps),1./b);
    ll = noise.nu_factor*(-sum(sum(sum(rho.*sum(abs(r),4),1),2),3) - D*sum(log(2*b)));

case {'binomial','binary'}
    ll =  noise.nu_factor*(sum(sum(sum(a0(msk).*f1(msk) - rho(msk).*log(1+exp(a0(msk))),3),2),1));

case {'multinomial','categorical'}
    [~,ll] = SoftMax(a0,f1);
    ll     = noise.nu_factor*ll;

otherwise
    error('Unknown likelihood function.');
end


%==========================================================================
%
%==========================================================================
function ll = ComputeLL_native(f,iphi,a0,s,noise)
% Compute ll of image in native space
d   = [size(a0) 1 1]; d = d(1:4);
a1  = Resamp(a0,iphi,s.bs_args);

msk = all(isfinite(f),4) & all(isfinite(a1),4);
switch lower(s.likelihood)
case {'normal','gaussian'}
    D   = sum(msk(:));
    ll  = 0;
    for l=1:d(4)
        a1l = a1(:,:,:,l);
        fl  = f(:,:,:,l);
        ll  = ll + noise.nu_factor*0.5*(D*(log(noise.lam(l)) - log(2*pi)) - noise.lam(l)*sum(sum(sum((fl(msk)-a1l(msk)).^2))));
    end

case {'laplace'}
    b   = reshape(sqrt(1./(2*noise.lam)),[1,1,1,d(4)]);
    D   = sum(msk(:));
    r   = bsxfun(@times,a1-f,1./b);
    ad  = sum(abs(r),4);
    ll  = noise.nu_factor*(-sum(ad(msk)) - D*sum(log(2*b))); % Check

case {'binomial','binary'}
    a1  = min(a1,40);
    ll  = noise.nu_factor*(sum(sum(sum(a1(msk).*f(msk) - log(1+exp(a1(msk))),3),2),1));

case {'multinomial','categorical'}
    [~,ll] = SoftMax(a1,f);
    ll     = noise.nu_factor*ll;

otherwise
    error('Unknown likelihood function.');
end


