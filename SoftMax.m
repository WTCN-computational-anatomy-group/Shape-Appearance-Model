function [sig,ll] = SoftMax(a1,f)
% Soft-Max likelihoods
% FORMAT [sig,ll] = SoftMax(a1,f)
%
% a1  - Model parameters
% f   - Data
%
% sig - Softmax of a1
% ll  - Log-likelihood
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

ll  = 0;
mxa = max(a1,[],4);
a1  = bsxfun(@minus,a1,mxa);
sig = exp(a1);
s   = sum(sig,4);

if nargin>=2 && nargout>=2
    f(~isfinite(f)) = 0;
    ll = sum(sum(sum(sum(f.*a1,4)-log(s).*sum(f,4), 3),2),1);
end
sig = bsxfun(@rdivide,sig,s);
