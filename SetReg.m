function [A,B,lb_qA,lb_pA] = SetReg(EZZ,N,~)
% Set the regularisation of Z
% FORMAT [A,B,lb_qA,lb_pA] = SetReg(EZZ,N,s)
%
% EZZ   - Expectation of Z*Z' (where Z encodes the modes of the estimates)
% N     - Number of observations
% s     - Settings. Uses s.nu0 & s.Lambda0
%
% A     - Expectation of A
% B     - B regularisation
% lb_qA - Lower bound stuff (needs more work)
% lb_pA - Lower bound stuff (needs more work)
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

K       = size(EZZ,1);
%if nargin<3, s = struct; end
B       = eye(K);
A       = inv((EZZ+max(diag(EZZ))*1e-5*eye(size(EZZ)))/N);

if false
if isfield(s,'nu0')
    % Get the Wishart Priors
    nu0 = s.nu0;
    if isfield(s,'Lambda0')
        Lambda0  = s.Lambda0;
    else
        Lambda0  = eye(K)/nu0;
    end

    nu     = N + nu0;
    Laminv = inv(double(Lambda0))+EZZ;
else
    % Maximum likelihood
    nu     = N;
    Laminv = EZZ;
end

% Stable inverse, as no ARD-style pruning is used
[V,D]    = eig(Laminv);
D        = diag(D);
D        = diag(max(D,max(D)*K*1e-6));
Lambda   = real(V*inv(D)*V');
%Lambda   = diag(diag(Lambda));
%reg      = mean(diag(Laminv))*eye(size(Laminv));
%Lambda   = inv((1-0.01)*Laminv + 0.01*reg);

A        = Lambda*nu;
end
if false
    % Unused stuff - possibly for fixing later
    ElndetA  = Elogdet(Lambda,nu);
    lb_qA    = 0.5*(nu -K-1)*ElndetA - 0.5*trace(Lambda \A)...
              -0.5*nu *K*log(2) - 0.5*nu *LogDet(Lambda ) - MultiGammaLn(nu /2,K);
    lb_pA    = 0.5*(nu0-K-1)*ElndetA - 0.5*trace(Lambda0\A)...
              -0.5*nu0*K*log(2) - 0.5*nu0*LogDet(Lambda0) - MultiGammaLn(nu0/2,K);
else
    lb_qA = 0;
    lb_pA = 0;
end

