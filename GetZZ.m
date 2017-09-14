function varargout = GetZZ(dat)
% Compute Z Z^T
% FORMAT [N,Z,ZZ,S] = GetZZ(dat)
%
% dat - 
% N   - 
% ZZ  - 
% S   -
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

Z  = single(0);
ZZ = single(0);
S  = single(0);
N  = numel(dat);
for n=1:N
    z  = dat(n).z;
    Z  = Z + z;
    ZZ = ZZ + z*z';
    S  = S + dat(n).S;
end
varargout{1} = N;
varargout{2} = Z;
varargout{3} = ZZ;
varargout{4} = S;

