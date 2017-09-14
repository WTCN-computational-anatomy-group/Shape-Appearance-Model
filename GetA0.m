function a0 = GetA0(z,Wa,mu)
% Compute mean plus linear combination of appearance bases
% FORMAT a0 = GetA0(z,Wa,mu)
% z  - Latent variables (weights for linear combination)
% Wa - Appearance basis functions
% mu - Appearance mean
% a0 - mean plus linear combination of appearance basis functions
%
% If z is a cell array of vectors, a0 will be a cell array of
% the corresponding images.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

d = [size(Wa) 1 1 1];
d = d(1:5);

if ~iscell(z) % z is a vector
    if nargin<3
        a0 = zeros(d(1:4),'single');
    else
        a0 = mu;
    end
    for k=1:d(5),
        a0 = a0 + Wa(:,:,:,:,k)*z(k);
    end
else         % z is a cell array
    a0 = cell(size(z));
    if nargin<3
        % Create one volume and copy it. May save memory.
        a0{1} = zeros(d(1:4),'single');
        for n=2:numel(a0)
            a0{n} = a0{1};
        end
    else
        % Create one volume and copy it. May save memory.
        a0{1} = mu(:,:,:,:);
        for n=2:numel(a0)
            a0{n} = a0{1};
        end
    end
    for k=1:d(5),
        Wk = Wa(:,:,:,:,k);
        for n=1:numel(a0)
            a0{n} = a0{n} + Wk*z{n}(k);
        end
    end
end

