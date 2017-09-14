function G = CompGrads(a0,bs_args)
% Compute gradients of an image
% FORMAT G = CompGrads(a0,bs_args)
% a0 - The input image or images.
% bs_args - Interpolation model (see spm_bsplinc/spm_bsplins)
% G       - A D x 3 cell array of gradients
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

d  = [size(a0) 1 1];
id = Identity(d);
G  = cell(d(4),3);

if nargin<2
    bs_args = [1 1 1  1 1 1];
end

for l=1:d(4)
    ac             = spm_diffeo('bsplinc',a0(:,:,:,l),[bs_args(1)   0 0 bs_args(4:6)]);
    [~,G{l,1},~,~] = spm_diffeo('bsplins',ac,id,      [bs_args(1)+1 0 0 bs_args(4:6)]);

    ac             = spm_diffeo('bsplinc',a0(:,:,:,l),[0 bs_args(2)   0 bs_args(4:6)]);
    [~,~,G{l,2},~] = spm_diffeo('bsplins',ac,id,      [0 bs_args(2)+1 0 bs_args(4:6)]);

    ac             = spm_diffeo('bsplinc',a0(:,:,:,l),[0 0 bs_args(3)   bs_args(4:6)]);
    [~,~,~,G{l,3}] = spm_diffeo('bsplins',ac,id,      [0 0 bs_args(3)+1 bs_args(4:6)]);
end

