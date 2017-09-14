function [f1,rho] = Push(f,iphi)
% Push an image accorging to a spatial transform
% FORMAT [f1,rho] = Push(f,iphi)
%
% f    - Image (3D or 4D)
% iphi - Spatial transform
%
% f1   - "Pushed" image
% rho  - Count of voxels pushed
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if ~isempty(iphi)
    [f1,rho] = spm_diffeo('pushc',single(f),iphi);
else
    f1  = f;
    msk = all(isfinite(f1),4);
    rho = single(msk);
end
f1(~isfinite(f1))   = 0;
rho(~isfinite(rho)) = 0;

