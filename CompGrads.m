function G = CompGrads(a)
% Compute gradients of an image
% FORMAT G = CompGrads(a)
%
% a - The input image or images.
% G - A D x 3 cell array of gradients
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

d  = [size(a) 1 1];
id = Identity(d);
G  = cell(d(4),3);

for l=1:d(4)
    [~,G{l,1},~,~] = spm_diffeo('bsplins',a(:,:,:,l),id, [2 0 0 1 1 1]);
    [~,~,G{l,2},~] = spm_diffeo('bsplins',a(:,:,:,l),id, [0 2 0 1 1 1]);
    [~,~,~,G{l,3}] = spm_diffeo('bsplins',a(:,:,:,l),id, [0 0 2 1 1 1]);
end

