function a1 = Resamp(a0,iphi,bs_args)
% resample an image
% FORMAT a1 = Resamp(a0,iphi,bs_args)
%
% a0      - Input image(s)
% iphi    - Deformation
% bs_args - B-spline settings
%
% a1      - Output image(s)
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if isempty(a0)
    a1 = a0;
elseif isempty(iphi)
    a1 = a0;
else
    d  = [size(a0) 1 1];
    a1 = zeros(d,'single');
    for l=1:d(4)
        c           = spm_diffeo('bsplinc',a0(:,:,:,l),bs_args);
        a1(:,:,:,l) = spm_diffeo('bsplins',c,iphi,bs_args);
    end
end

