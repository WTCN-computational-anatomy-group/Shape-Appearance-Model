function [gmu,Hmu,nll] = muGradHess(dat,mu,Wa,Wv,noise,s)
% Return derivatives w.r.t. appearance basis functions
% FORMAT [nll,gmu,Hmu] = WaGradHess(dat,mu,Wa,Wv,noise,s)
%
% dat   - Structure containing various information about each image.
%         Fields for each image n are:
%         dat(n).f - Image data.
%         dat(n).z - Expectations of latent variables.
%         dat(n).S - Covariances of latent variables.
% mu    - Mean
% Wa    - Appearance basis functions
% Wv    - Shape basis functions
% noise - Noise information
% s     - Settings. Uses s.likelihood and s.batchsize.
%
% gmu   - Gradients
% Hmu   - Hessians
% nll   - Negative log-likelihood
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if isempty(dat), gmu = []; Hmu = []; nll = 0; return; end

d  = [size(mu), 1,1];
d  = d(1:4);

switch lower(s.likelihood)
case {'normal','gaussian','binomial','binary'}
    d4 = d(4);
case {'multinomial','categorical'}
    d4 = d(4)*(d(4)+1)/2;
end

batchsize = 1;
if isfield(s,'batchsize'), batchsize = s.batchsize; end
gmu = zeros( d(1:4)    ,'single');
Hmu = zeros([d(1:3) d4],'single');

nll = 0; 

for n1=1:batchsize:numel(dat)
    nn    = n1:min(n1+(batchsize-1),numel(dat));
    z     = {dat(nn).z};
    cell1 = GetV0(z,Wv);
    cell2 = GetA0(z,Wa,mu);
    dat1  = dat(nn);

    parfor n=1:numel(nn)
        iphi     = GetIPhi(cell1{n},s);
        a0       = cell2{n};
        f        = GetDat(dat1(n),s);
        [ll,g,H] = AppearanceDerivs(f,a0,iphi,noise,s);
        nll      = nll - ll;
        gmu      = gmu + g;
        Hmu      = Hmu + H;
    end
end

