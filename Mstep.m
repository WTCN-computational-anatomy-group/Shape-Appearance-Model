function [Wa,Wv,WW,omega] = Mstep(mu,Wa,Wv,noise,B,ZZ,A,WWa,WWv,s)
% Update the basis functions
% FORMAT  [Wa,Wv,WW,omega] = Mstep(mu,Wa,Wv,noise,B,ZZ,A,WWa,WWv,s)
%
% mu    - Mean.
% Wa    - Appearance basis functions.
% Wv    - Shape basis functions.
% noise - Noise information.
% B     - Precision matrix for W (sort of).
% ZZ    - Z*Z' sufficient statistic.
% A     - Expectation of precision matrix for Z.
% WWa   - Wa'*La*Wa.
% WWv   - Wv'*Lv*Wv.
% s     - Settings. Uses s.likelihood, s.ondisk, s.result_dir and
%         s.result_name.
%
% Wa    - Appearance basis functions.
% Wv    - Shape basis functions.
% WW    - Combined WWa and WWv.
% omega - step size for Gauss-Newton update of Wv
%
% Updates Wa and Wv via Gauss-Newton updates.  This normally involves two
% passes through the data, although more may be required if the Gauss-
% Newton updates of the shape bases functions (Wv) overshoot.  A new WW is
% also computed.  
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

Ka   = size(Wa,5);
inda = 1:Ka;
Kv   = size(Wv,5);
indv = 1:Kv;
K    = size(ZZ,1);
if Ka<K || Kv<K
    indv = indv + Ka;
end

RegW = double(s.lambda(1)*B  + s.lambda(2)*ZZ);

WW            = zeros(K);
WW(inda,inda) = WW(inda,inda) + WWa;
WW(indv,indv) = WW(indv,indv) + WWv;

if numel(indv)>0

    [gv,Hv,nll] = PGdistribute('ShapeDerivatives',mu,Wa,Wv,noise,s);

    if ~isempty(nll)

        if isfield(s,'ondisk') && s.ondisk
            prev.Wv       = Wv;
            prev.Wv.fname = fullfile(s.result_dir,[s.result_name '_prevWv.dat']);
            for k=1:numel(indv)
                prev.Wv(:,:,:,:,k) = Wv(:,:,:,:,k);
            end
        else
            prev.Wv = Wv;
        end
        prev.WWv    = WWv;

        nll0        = nll+0.5*s.lambda(1)*(trace(WW*B) + trace(A*ZZ)) + 0.5*s.lambda(2)*trace(WW*ZZ);
        fprintf('%9.6g ', -nll0);

        for subit=1:8
            Wv            = UpdateWv(Wv,gv,Hv,RegW(indv,indv),s);
            WWv           = UpdateWWv(Wv,s);
            WW            = zeros(K);
            WW(inda,inda) = WW(inda,inda) + WWa;
            WW(indv,indv) = WW(indv,indv) + WWv;

            if numel(inda)>0
                [ga,Ha,nll1] = PGdistribute('AppearanceDerivatives',mu,Wa,Wv,noise,s);
            else
                nll1         = PGdistribute('ComputeOF',mu,Wa,Wv,noise,s);
            end


            nll              = nll1 + 0.5*s.lambda(1)*(trace(WW*B) + trace(A*ZZ)) + 0.5*s.lambda(2)*trace(WW*ZZ);

           %if ~isfinite(nll), error('Something went wrong.'); end
            if nll<max(nll0*0.999999,nll0*1.000001)
                s.omega = min(s.omega*1.1,1);
                break;
            else
                s.omega = max(s.omega*0.5,0.001);
                for k=1:numel(indv)
                    Wv(:,:,:,:,k) = prev.Wv(:,:,:,:,k);
                end
                WWv   = prev.WWv;
            end
        end
        fprintf('%9.6g ', -nll);
    end
else
    if numel(inda)>0
        [ga,Ha]   = PGdistribute('AppearanceDerivatives',mu,Wa,Wv,noise,s);
    end
end
if numel(inda) > 0
    Wa            = UpdateWa(Wa,ga,Ha,RegW(inda,inda),s);
    WWa           = UpdateWWa(Wa,s);
end

WW            = zeros(K);
WW(inda,inda) = WW(inda,inda) + WWa;
WW(indv,indv) = WW(indv,indv) + WWv;

omega         = s.omega;

if isfield(s,'ondisk') && s.ondisk && exist('prev','var')
    DeleteFA(prev.Wv,gv,Hv);
    if numel(inda)>0
        DeleteFA(ga,Ha);
    end
end
 
