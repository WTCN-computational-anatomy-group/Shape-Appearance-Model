function [dat,stats] = UpdateLatentVariables(dat,mu,Wa,Wv,noise,A,s)
% Re-estimate all latent variables and return useful suffient statistics
% FORMAT [dat,stats] = UpdateLatentVariables(dat,mu,Wa,Wv,noise,A,s)
%
% dat   - Structure containing various information about each image.
%         Fields for each image n are:
%         dat(n).f - Image data.
%         dat(n).z - Expectations of latent variables.
%         dat(n).S - Covariances of latent variables.
% mu    - Mean image
% Wa    - Appearance basis functions
% Wv    - Shape basis functions
% noise - Noise information
% A     - Precision of Z
% s     - Settings
%
% stats - An assortment of statistics from Z
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

%CompSmo = false;
%CompMu  = false;

if isempty(dat), stats=[]; return; end

batchsize = 1;
if isfield(s,'batchsize'), batchsize = s.batchsize; end

stats = struct('N',numel(dat), 'Z',0, 'ZZ',0, 'S',0, 's0',0, 's1',0, 'L', 0);

%if CompSmo
%    stats.SmoSuf = 0;
%end
%if CompMu
%    d            = size(mu);
%    stats.gmu    = zeros(d,'single');
%    stats.Hmu    = single(0);
%end

for n1=1:batchsize:numel(dat)
    nn    = n1:min(n1+(batchsize-1),numel(dat));
    z     = {dat(nn).z};
    f     = cell(size(z));
    for n=1:numel(nn)
        f{n} = GetDat(dat(nn(n)),s);
    end

    [z,S,l,omisc] = UpdateZpar(z,f,mu,Wa,Wv,A,s,noise);
    [dat(nn).z] = z{:};
    [dat(nn).S] = S{:};

    z      = double(cat(2,z{:}));
    S      = sum(double(cat(3,S{:})),3);

    stats.L      = stats.L      + l;
    stats.s1     = stats.s1     + omisc.s1;
    stats.s0     = stats.s0     + omisc.s0;
    stats.Z      = stats.Z      + sum(z,2);
    stats.ZZ     = stats.ZZ     + z*z';
    stats.S      = stats.S      + S;
    %if CompSmo
    %    stats.SmoSuf = stats.SmoSuf + omisc.SmoSuf;
    %end
    %if CompMu
    %    stats.gmu    = stats.gmu    + omisc.gmu;
    %    stats.Hmu    = stats.Hmu    + omisc.Hmu;
    %end
end

