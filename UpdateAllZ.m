function [dat,stats] = UpdateAllZ(dat,mu,Wa,Wv,noise,A,s)
% Re-estimate all latent variables and compute useful suffient statistics
% FORMAT [dat,stats] = UpdateAllZ(dat,mu,Wa,Wv,noise,A,s)
%
% dat   - Data
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

if isempty(dat), stats=[]; return; end
batchsize = 1;
if isfield(s,'batchsize'), batchsize = s.batchsize; end

d      = size(mu);
gmu    = zeros(d,'single');
Hmu    = single(0);
L      = 0;
s1     = 0;
s0     = 0;
SmoSuf = 0;
sS     = 0;
ZZ     = 0;
Z      = 0;

for n1=1:batchsize:numel(dat)
    nn    = n1:min(n1+(batchsize-1),numel(dat));
    z     = {dat(nn).z};
    f     = {};
    for n=nn
        f = [f, {GetDat(dat(n),s)}];
    end

    [z,S,l,omisc] = UpdateZpar(z,f,mu,Wa,Wv,A,s,noise);
    [dat(nn).z] = z{:};
    [dat(nn).S] = S{:};

    z      = double(cat(2,z{:}));
    S      = sum(double(cat(3,S{:})),3);

    L      = L      + l;
    gmu    = gmu    + omisc.gmu;
    Hmu    = Hmu    + omisc.Hmu;
    s1     = s1     + omisc.s1;
    s0     = s0     + omisc.s0;
    SmoSuf = SmoSuf + omisc.SmoSuf;
    Z      = Z      + sum(z,2);
    ZZ     = ZZ     + z*z';
    sS     = sS     + S;
end

stats = struct('N',numel(dat), 'Z',Z, 'ZZ',ZZ,...
               'sS',sS, 'L',L, 'gmu',gmu, 'Hmu',Hmu,...
               's0',s0, 's1',s1, 'SmoSuf',SmoSuf);




function [dat,stats] = UpdateAllZnonpar(dat,mu,Wa,Wv,noise,A,s)
% Re-estimate all latent variables and compute useful suffient statistics
% FORMAT [dat,stats] = UpdateAllZ(dat,mu,Wa,Wv,noise,A,s)
%
% dat   - Data
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

if isempty(dat), stats=[]; return; end

d      = size(mu);
gmu    = zeros(d,'single');
Hmu    = single(0); %zeros(d,'single');
L      = 0;
s1     = 0;
s0     = 0;
SmoSuf = 0;
sS     = 0;
ZZ     = 0;
Z      = 0;
for n=1:numel(dat),
    f      = GetDat(dat(n),s);
    z      = dat(n).z;
    [z,S,l,omisc] = UpdateZ(z,f,mu,Wa,Wv,A,s,noise);
    dat(n).z  = z;
    dat(n).S  = S;

    L      = L      + l;
    gmu    = gmu    + omisc.gmu;
    Hmu    = Hmu    + omisc.Hmu;
    s1     = s1     + omisc.s1;
    s0     = s0     + omisc.s0;
    SmoSuf = SmoSuf + omisc.SmoSuf;
    Z      = Z      + double(z);
    ZZ     = ZZ     + double(z*z');
    sS     = sS     + double(S);
end

stats = struct('N',numel(dat), 'Z',Z, 'ZZ',ZZ,...
               'sS',sS, 'L',L, 'gmu',gmu, 'Hmu',Hmu,...
               's0',s0, 's1',s1, 'SmoSuf',SmoSuf);

