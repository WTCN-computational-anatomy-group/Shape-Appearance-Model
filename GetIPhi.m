function iphi = GetIPhi(varargin)
% Compute a diffeomorphism from shape basis functions
% FORMAT iphi = GetIPhi(z,Wv,s)
%
% z    - Latent variables
% Wv   - Shape basis functions
% s    - Settings. Uses s.vx, s.v_settings & s.int_args.
% iphi - The resulting diffeomorphism
%
% FORMAT iphi = GetIPhi(v0,s)
%
% v0   - Initial velocity
% s    - Settings. Uses s.vx, s.v_settings & s.int_args.
% iphi - The resulting diffeomorphism
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if nargin==3 % && isstruct(varargin{3})
    z   = varargin{1};
    Wv  = varargin{2};
    s   = varargin{3};
    v0  = GetV0(z,Wv);
elseif nargin==2
    v0  = varargin{1};
    s   = varargin{2};
else
    error('Unknown option.');
end
if ~iscell(v0)
    iphi = Shoot(v0,[s.vx s.v_settings],s.int_args);
else
    for n=1:numel(v0)
        v0{n} = Shoot(v0{n},[s.vx s.v_settings],s.int_args);
    end
    iphi = v0;
end

