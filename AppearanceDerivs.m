function [ll,g,H] = AppearanceDerivs(f,a0,iphi,noise,s)
% Compute gradients and Hessian w.r.t. voxelwise changes in appearance
% FORMAT [ll,g,H] = AppearanceDerivs(f,a0,iphi,noise,s)
%
% f     - image data (3/4D single precision)
% a0    - atlas data (3/4D single)
% iphi  - Deformation
% noise - noise precision (for Gaussian noise only)
% s     - settings
% s.likelihood - either 'normal', 'binomial' or 'multinomial'.
%
% ll  - log-likelihood
% g   - Voxel-wise 1st derivatives 
% H   - Voxel-wise 2nd derivatives
%
% Note that for multi-channel images, H encodes a diagonal hessian matrix
% except for 'multinomial', where the diagonal followed by the off-diagonal
% elements are returned.
%
% Best to read the code to understand what it does.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

d   = [size(a0) 1 1 1];
d   = d(1:4);
a1  = Pull(a0,iphi);
msk = isfinite(f) & isfinite(a1);

switch lower(s.likelihood)
case {'normal','gaussian'}
    % The L2 norm works reasonably well, but future versions should probably make use
    % of non-stationary variances.
    lam = reshape(noise.lam,[1,1,1,d(4)]);
    ll  = 0;
    for l=1:d(4)
        a1l  = a1(:,:,:,l);
        fl   = f(:,:,:,l);
        mskl = msk(:,:,:,l);
        D    = sum(mskl(:));
        ll   = ll + noise.nu_factor*0.5*(D*(log(noise.lam(l)) - log(2*pi)) - noise.lam(l)*sum(sum(sum((fl(mskl)-a1l(mskl)).^2))));
    end
    if nargin<=1, return; end

    r   = a1-f;
    g   = bsxfun(@times,noise.nu_factor*lam,r);
    H   = bsxfun(@times,noise.nu_factor*lam,single(msk));

case {'binomial','binary'}
    % For simple binary (or almost binary) images
    a1  = min(a1,40);
    ea  = exp(a1);
    ll  = noise.nu_factor*(sum(sum(sum(a1(msk).*f(msk) - log(1+ea(msk)),3),2),1));
    if nargin<=1, return; end

    sig = ea./(1+ea);
    sig(~msk) = NaN;
    g   = noise.nu_factor*(sig-f);
    H   = noise.nu_factor*((sig.*(1-sig) + 1e-3));

case {'multinomial','categorical'}
    % Could be made more efficient in terms of disk I/O etc.  This representation
    % has some redundancy, which could be eliminated using an approach like the
    % one in spm12/toolbox/Shoot/spm_shoot_blur.m , which works with the subspace
    % null(ones(1,d(4))).  This approach may also increase stability as the Hessian
    % of the likelihood would not be singular.
    [sig,ll] = SoftMax(a1,f);
    ll       = noise.nu_factor*ll;
    if nargin<=1, return; end

    sig(~msk) = NaN;
    g   = sig-f;
    H   = zeros([d(1:3) d(4)*(d(4)+1)/2],'single');
    for l1 = 1:d(4)
        H(:,:,:,l1) = sig(:,:,:,l1) - sig(:,:,:,l1).^2 + 1e-4*d(4)^2;
    end
    l = d(4);
    for l1 = 1:d(4)
        for l2 = (l1+1):d(4)
            l          = l + 1;
            H(:,:,:,l) = -sig(:,:,:,l1).*sig(:,:,:,l2);
        end
    end
    g   = noise.nu_factor*g;
    H   = noise.nu_factor*H;
otherwise
    error('Unknown likelihood function.');
end

g = Push(g,iphi);
H = Push(H,iphi);

