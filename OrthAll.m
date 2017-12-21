function [Wa,Wv,WWa,WWv,ss,WW] = OrthAll(Wa,Wv,WWa,WWv,ss,s)
% Orthogonalise everything
% FORMAT [Wa,Wv,WWa,WWv,ss,WW] = OrthAll(Wa,Wv,WWa,WWv,ss,s)
%
% Wa    - Appearance basis functions.
% Wv    - Shape basis functions.
% WWa   - Wa'*La*Wa.
% WWv   - Wv'*Lv*Wv.
% ss    - Sufficient statistics, which includes ss.ZZ, which is computed
%         from Z*Z'.
% s     - Settings.
%
% Wa    - Appearance basis functions.
% Wv    - Shape basis functions.
% WWa   - Wa'*La*Wa.
% WWv   - Wv'*Lv*Wv.
% ss    - Sufficient statistics, which includes ss.ZZ, which is computed
%         from Z*Z'.
% WW    - Combined WWa and WWv.
%
% Computes a suitable transform, T, which is used to transform W and Z so
% that WW and ZZ are diagonal:
%     W        <- W*inv(T)
%     WW       <- inv(T)'*WW*inv(T)
%     dat(n).z <- T*dat(n).z
%     dat(n).S <- T*dat(n).S*T'
%     ss.Z     <- T*ss.Z 
%     ss.ZZ    <- T*ss.ZZ*T'
%     ss.sS    <- T*ss.sS*T'
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


d             = [size(Wa) 1 1 1];
dv            = [size(Wv) 1 1 1];
Ka            = d(5);
Kv            = dv(5);
K             = size(ss.ZZ,1);
if Ka<K || Kv<K,
    Koff   = Ka;
    linked = false;
else
    Koff   = 0;
    linked = true;
end
inda  = 1:Ka;
indv  = Koff+(1:Kv);

WW            = zeros(K);
WW(inda,inda) = WW(inda,inda) + WWa;
WW(indv,indv) = WW(indv,indv) + WWv;

[T,iT] = OrthogonalisationMat(ss.ZZ,ss.sS,WW,ss.N,s);

if linked
    iTa = iT;
    iTv = iT;
else
    iTa = iT(inda,inda);
    iTv = iT(indv,indv);
end

for i=1:size(Wa,3)
    Wa(:,:,i,:,:) = reshape(reshape(Wa(:,:,i,:,:),[prod(d([1 2 4])) numel(inda)])*iTa,[d(1:2) 1 d(4) numel(inda)]);
end
WWa = iTa'*WWa*iTa;
for i=1:size(Wv,3)
    Wv(:,:,i,:,:) = reshape(reshape(Wv(:,:,i,:,:),[prod(d(1:2))*3   numel(indv)])*iTv,[d(1:2) 1   3  numel(indv)]);
end
WWv           = iTv'*WWv*iTv;
WW            = zeros(K);
WW(inda,inda) = WW(inda,inda) + WWa;
WW(indv,indv) = WW(indv,indv) + WWv;

ss.Z     = T*ss.Z;
ss.ZZ    = T*ss.ZZ*T';
ss.sS    = T*ss.sS*T';
PGdistribute('TransfZ',T);

