function [sig,ll] = SoftMax(a0,f1)
% Soft-Max likelihoods
% FORMAT [sig,ll] = SoftMax(a0,f1)
%
% a0  - Model parameters
% f1  - Data
%
% sig - Softmax of a0
% ll  - Log-likelihood
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

ll  = 0;
mxa = max(a0,[],4);
a0  = bsxfun(@minus,a0,mxa);
sig = exp(a0);
s   = sum(sig,4);
if nargin>=2 && nargout>=2
    f1(~isfinite(f1)) = 0;
    ll = sum(sum(sum(sum(f1.*a0,4)-log(s).*sum(f1,4), 3),2),1);
end
sig = bsxfun(@rdivide,sig,s);

