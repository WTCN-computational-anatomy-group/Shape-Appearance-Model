function v0 = GetV0(z,Wv)
% Compute velocity field(s) from shape basis functions
% FORMAT v0 = GetV0(z,Wv)
%
% z  - Latent variables
% Wv - Shape basis functions
%
% v0 - The resulting velocity field
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

d   = [size(Wv) 1];
d   = d(1:5);
if ~iscell(z)
    if prod(d)==0
        v0 = [];
    else
        Koff = numel(z)-d(5);
        v0   = zeros([d(1:3),3],'single');
        for k=1:d(5)
            v0 = v0 + z(k+Koff)*Wv(:,:,:,:,k);
        end
    end
else
    v0 = cell(size(z));
    if prod(d)==0
        for n=1:numel(v0)
            v0{n} = [];
        end
    else
        v0{1} = zeros([d(1:3) 3],'single');
        for n=2:numel(v0)
            v0{n} = v0{1};
        end
        Koff = numel(z{1})-d(5);
        for k=1:d(5)
            Wk = Wv(:,:,:,:,k);
            for n=1:numel(v0)
                v0{n} = v0{n} + Wk*z{n}(k+Koff);
            end
        end
    end
end

