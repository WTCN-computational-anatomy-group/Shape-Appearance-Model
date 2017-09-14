function dat = AddToZ(dat,za)
% Addition of vector to latent variables
% FORMAT dat = AddToZ(dat,za)
%
% dat - data (with fields dat(n).z & dat(n).S)
% za  - vector to add
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

N = numel(dat);
for n=1:N
    dat(n).z = dat(n).z + za;
end

