function I = Horder(d)
% Identify where to find relevant elements of Hessians.
% FORMAT I = Horder(d)
%
% d - dimension of Hessian matrix.
%
% I - Lookup table for Hessian elements.
%
% Only the upper triangle of the hessian matrix fields are saved. This
% function identifies where to find each element.
%     e.g. H(i,j) = stored(I(i,j))
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

I = diag(1:d);
l = d;
for i1=2:d
    for i2=1:(i1-1)
        l = l + 1;
        I(i1,i2) = l;
        I(i2,i1) = l;
    end
end

