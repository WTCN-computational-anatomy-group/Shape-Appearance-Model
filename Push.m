function [f1,rho] = Push(f,psi)
% Push an image (or set of images) accorging to a spatial transform
% FORMAT [f1,rho] = Push(f,psi)
%
% f    - Image (3D or 4D)
% psi - Spatial transform
%
% f1   - "Pushed" image
% rho  - Count of voxels pushed
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if ~isempty(psi)
    [f1,rho] = spm_diffeo('pushc',single(f),psi);
else
    f1  = f;
    msk = all(isfinite(f1),4);
    rho = single(msk);
end
f1(~isfinite(f1))   = 0;
rho(~isfinite(rho)) = 0;

