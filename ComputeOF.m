function nll = ComputeOF(dat,mu,Wa,Wv,noise,s)
% Compute part of objective function
% FORMAT nll = ComputeOF(dat,mu,Wa,Wv,noise,s)
% dat   - 
% mu    - 
% Wa    - 
% Wv    - 
% noise - 
% s     - 
% nll   - 
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

batchsize = 1;
if isfield(s,'batchsize'), batchsize = s.batchsize; end

nll = 0;
for n1=1:batchsize:numel(dat)
    nn    = n1:min(n1+(batchsize-1),numel(dat));
    z     = {dat(nn).z};
    cell1 = GetV0(z,Wv);
    cell2 = GetA0(z,Wa,mu);
    dat1  = dat(nn);

    parfor n=1:numel(nn)
        iphi = GetIPhi(cell1{n},s);
        nll  = nll - AppearanceDerivs(GetDat(dat1(n),s),cell2{n},iphi,noise,s);
    end
end

