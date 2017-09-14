function [g,H] = AppearanceDerivs(f1,rho,a0,noise,s)
% Compute gradients and Hessian w.r.t. voxelwise changes in appearance
% FORMAT [g,H] = AppearanceDerivs(f1,rho,a0,noise,s)
%
% f1    - warped (pushed) image data (3/4D single precision)
% rho   - counts of original voxels (3D single)
% a0    - atlas data (3/4D single)
% noise - noise precision (for Gaussian noise only)
% s     - settings
% s.likelihood - either 'normal', 'laplace', 'binomial' or 'multinomial'.
%
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

d = [size(a0) 1 1 1];
d = d(1:4);

switch lower(s.likelihood)
case {'normal','gaussian'}
    % The L2 norm works reasonably well, but future versions should probably make use
    % of non-stationary variances.
    lam = reshape(noise.lam,[1,1,1,d(4)]);
    r   = bsxfun(@times,rho,a0)-f1;
    g   = bsxfun(@times,noise.nu_factor*lam,r);
    H   = bsxfun(@times,noise.nu_factor*lam,rho);

case {'laplace'}
    % The L1 norm does not work so well in practice.  Maybe its a bug, or maybe
    % its simply not a great approach.
    lam = reshape(noise.lam,[1,1,1,d(4)]);
    r   = a0 - bsxfun(@rdivide,f1,rho+1e-6);
    q   = 2./sqrt(bsxfun(@times,lam/2,r.^2)+0.01);
    H   = bsxfun(@times,bsxfun(@times,noise.nu_factor*lam,q),rho);
    g   = H.*r;

case {'binomial','binary'}
    % For simple binary (or almost binary) images
    ea  = exp(a0);
    sig = ea./(1+ea);
    g   = noise.nu_factor*(rho.*sig-f1);
    H   = noise.nu_factor*(rho.*(sig.*(1-sig) + 1e-3));

case {'multinomial','categorical'}
    % Could be made more efficient in terms of disk I/O etc.  This representation
    % has some redundancy, which could be eliminated using an approach like the
    % one in spm12/toolbox/Shoot/spm_shoot_blur.m , which works with the subspace
    % null(ones(1,d(4))).  This approach may also increase stability as the Hessian
    % of the likelihood would not be singular.
    sig = SoftMax(a0);
    g   = bsxfun(@times,rho,sig)-f1;
    H   = zeros([d(1:3) d(4)*(d(4)+1)/2],'single');
    for l1 = 1:d(4)
        H(:,:,:,l1) = rho.*(sig(:,:,:,l1) - sig(:,:,:,l1).^2 + 1e-3*d(4)^2);
    end
    l = d(4);
    for l1 = 1:d(4)
        for l2 = (l1+1):d(4)
            l          = l + 1;
            H(:,:,:,l) = -rho.*sig(:,:,:,l1).*sig(:,:,:,l2);
        end
    end
    g   = noise.nu_factor*g;
    H   = noise.nu_factor*H;
otherwise
    error('Unknown likelihood function.');
end

g(~isfinite(g)) = 0;
H(~isfinite(H)) = 0;

