function [A,B,lb_qA,lb_pA,ElndetA] = SetReg(EZZ,N,s)
% Set the regularisation of Z
% FORMAT [A,B,lb_qA,lb_pA,ElndetA] = SetReg(EZZ,N,s)
%
% EZZ     - Expectation of Z*Z' (where Z encodes the modes of the estimates)
% N       - Number of observations
% s       - Settings. May use s.nu0 & s.Lambda0 if available.
%
% A       - Expectation of A
% B       - B regularisation
% lb_qA   - Lower bound stuff (needs more work)
% lb_pA   - Lower bound stuff (needs more work)
% ElndetA - E[log(det(A))]
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

K       = size(EZZ,1);
B       = eye(K)*N;

if nargin<3, s = struct; end

if isfield(s,'nu0') && s.nu0 > 0
    % Get the Wishart Priors
    nu0 = s.nu0;
    if isfield(s,'Lambda0')
        Lambda0  = s.Lambda0;
    else
        Lambda0  = eye(K)/nu0;
    end

    nu     = N + nu0;
    Laminv = inv(double(Lambda0))+EZZ;
    Lambda = inv(Laminv);
    A      = Lambda*nu;

    if s.nu0 > K
        ElndetA  = Elogdet(Lambda,nu);
        lb_qA    = 0.5*(nu -K-1)*ElndetA - 0.5*trace(Lambda \A) -0.5*nu *K*log(2) - 0.5*nu *LogDet(Lambda ) - MultiGammaLn(nu /2,K);
        lb_pA    = 0.5*(nu0-K-1)*ElndetA - 0.5*trace(Lambda0\A) -0.5*nu0*K*log(2) - 0.5*nu0*LogDet(Lambda0) - MultiGammaLn(nu0/2,K);
    else
        ElndetA = LogDet(A);
        lb_qA   = 0;
        lb_pA   = 0;
    end

else
    A       = inv((EZZ+max(diag(EZZ))*1e-5*eye(size(EZZ)))/N);
    ElndetA = LogDet(A);
    lb_qA   = 0;
    lb_pA   = 0;
end

