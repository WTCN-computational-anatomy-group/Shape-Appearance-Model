function id = Identity(d)
% Identity transform for deformations
% FORMAT id = Identity(d)
%
% d  - Dimensions
%
% id - Identity transform (d(1) x d(2) x d(3) x 3 single)
%
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2013-2017)

% John Ashburner
% $Id$

d  = [d(:)' 1 1];
d  = d(1:3);
id = zeros([d,3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));

