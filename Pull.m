function a1 = Pull(a0,psi)
% Resample an image or set of images
% FORMAT a1 = Pull(a0,psi)
%
% a0      - Input image(s)
% psi    - Deformation
%
% a1      - Output image(s)
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if isempty(a0)
    a1 = a0;
elseif isempty(psi)
    a1 = a0;
else
    d  = [size(a0) 1 1];
    a1 = zeros(d,'single');
    for l=1:d(4)
        a1(:,:,:,l) = spm_diffeo('bsplins',a0(:,:,:,l),psi,[1 1 1  1 1 1]);
    end
end

