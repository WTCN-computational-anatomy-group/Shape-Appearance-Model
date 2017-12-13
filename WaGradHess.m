function [ga,Ha,nll] = WaGradHess(dat,mu,Wa,Wv,noise,s)
% Compute derivatives w.r.t. appearance basis functions
% FORMAT [ga,Ha,nll] = WaGradHess(dat,mu,Wa,Wv,noise,s)
% dat   - Array of data
% mu    - Mean
% Wa    - Appearance basis functions
% Wv    - Shape basis functions
% noise - Noise information
% s     - Settings. Uses s.likelihood, s.ondisk, s.result_dir & s.result_name
%
% ga    - Gradients
% Ha    - Hessians
% nll   - Negative log-likelihood
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if isempty(dat), ga = []; Ha = []; nll = 0; return; end

Ka  = size(Wa,5);
if Ka==0, ga = []; Ha = []; nll = []; return; end

d  = [size(mu), 1,1];
d  = d(1:4);

switch lower(s.likelihood)
case {'normal','laplace'}
    d4 = d(4);
case {'binomial','binary'}
    d4 = 1;
case {'multinomial','categorical'}
    d4 = d(4)*(d(4)+1)/2;
end

batchsize = 1;
if isfield(s,'batchsize'), batchsize = s.batchsize; end
if isfield(s,'ondisk') && s.ondisk
    ga = file_array(fullfile(s.result_dir,[s.result_name '_ga.dat']),[d(1:4)    Ka],'float32',352);
    Ha = file_array(fullfile(s.result_dir,[s.result_name '_Ha.dat']),[d(1:3) d4 Ka],'float32',352);
else
    ga = zeros([d(1:4)    Ka],'single');
    Ha = zeros([d(1:3) d4 Ka],'single');
end

nll = 0;

% Compute 1st and 2nd derivatives w.r.t. appearance model
for k=1:Ka,
    ga(:,:,:,:,k) = 0;
    Ha(:,:,:,:,k) = 0;
end

for n1=1:batchsize:numel(dat)
    nn    = n1:min(n1+(batchsize-1),numel(dat));
    z     = {dat(nn).z};
   %S     = {dat(nn).S};
    cell1 = GetV0(z,Wv);
    cell2 = GetA0(z,Wa,mu);
    dat1  = dat(nn);

    parfor n=1:numel(nn)
        iphi      = GetIPhi(cell1{n},s);
        a0        = cell2{n};
        f         = GetDat(dat1(n),s);
        [tmp,g,H] = AppearanceDerivs(f,a0,iphi,noise,s);
        nll       = nll - tmp;
        cell1{n}  = g;
        cell2{n}  = H;
    end

    for k=1:Ka,
        g1 = single(0);
        H1 = single(0);
        for n=1:numel(cell1)
            g1  = g1 + cell1{n}*z{n}(k);
            ezz = z{n}(k)^2;% + S{n}(k,k);
            H1  = H1 + cell2{n}*ezz;
        end
        ga(:,:,:,:,k) = ga(:,:,:,:,k) + g1;
        Ha(:,:,:,:,k) = Ha(:,:,:,:,k) + H1;
    end
end

