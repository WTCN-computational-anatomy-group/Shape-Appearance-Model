function noise = NoiseModel(ss,s,d)
% Compute noise model from sufficient statistics
% FORMAT noise = NoiseModel(ss,s,d)
%
% ss    - Data structure returned by UpdateZ. Uses the fields s0 & s1.
% s     - Settings. Uses the fields s.nu_factor, s.likelihood, s.alpha0
%         and s.beta0.
% d     - Image dimensions (not currently used).
%
% noise - A data structure with the following fields:
%         lam       - Reciprocal of the variance. Used only for the
%                     Gaussian noise model.
%         lb_lam    - Lower bound relating to the computation of lambda
%                     (Gaussian noise model only).
%         nu_factor - Fudge factor weighting the matching term likelihood.
%                     Can be used for accounting for image smoothness.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

CompSmo         = false; % Smoothness estimation is not currently used
noise.lam       = 1;
noise.lb_lam    = 0;
noise.nu_factor = 1;
if isfield(s,'nu_factor'), noise.nu_factor = s.nu_factor; end

if CompSmo
    % See the appendix of:
    % Groves AR, Beckmann CF, Smith SM, Woolrich MW. Linked independent component analysis
    % for multimodal data fusion. Neuroimage. 2011 Feb 1;54(3):2198-217.
    SmoSuf          = ss.SmoSuf;
    noise.fwhm      = sqrt(4*log(2)*(SmoSuf(2)/SmoSuf(1))/(sum(SmoSuf([4 6 8]))/sum(SmoSuf([3 5 7]))));
    noise.fwhm      = sqrt(max(noise.fwhm^2 - 0.5, 4*log(2)/pi)); % An adjustment that helps empirically
end

switch lower(s.likelihood)
case {'normal','gaussian'}

    if CompSmo
        % Note that a slightly different adjustment is used for the Gaussian noise model
        % This one is a heavier adjustment, which is based on taking iid data and smoothing
        % in order to estimate the original variance of the iid noise.
        noise.nu_factor = (sqrt(2*log(2)/pi)/noise.fwhm)^sum(d(1:3)>1);
    end

    alpha0 = 0.0001*ones(size(ss.s1)); if isfield(s,'alpha0'), alpha0 = s.alpha0; end
    beta0  = 0.0001*ones(size(ss.s1)); if isfield(s,'beta0' ), beta0  = s.beta0;  end

    noise.alpha = alpha0 + ss.s0*noise.nu_factor;
    noise.beta  = beta0  + ss.s1*noise.nu_factor;
    noise.lam   = noise.alpha./noise.beta;

    lb_lam  = 0;
    lb_qlam = 0;
    lb_plam = 0;
    for l=1:d(4)
        lb_lam  = lb_lam + 0.5*ss.s0(1)*(psi(noise.alpha(l)) - log(noise.beta(l))) ...
                         - 0.5*ss.s0(1)*log(2*pi) - 0.5*noise.lam(l)*ss.s1(1);

        lb_qlam = lb_qlam + noise.alpha(l) *log(noise.beta(l)) ...
                          + (noise.alpha(l)-1)*(psi(noise.alpha(l)) - log(noise.beta(l))) ...
                          -(noise.alpha(l)/noise.beta(l))*noise.beta(l) ...
                          - gammaln(noise.alpha(l));
        lb_plam = lb_plam + alpha0(l)*log(beta0(l)) ...
                          + (alpha0(l)-1)*(psi(noise.alpha(l)) - log(noise.beta(l))) ...
                          -(noise.alpha(l)/noise.beta(l))*beta0(l) ...
                          - gammaln(alpha0(l));
    end
    noise.lb_lam  = lb_lam + lb_plam - lb_qlam;

case {'binomial','binary','multinomial','categorical'}
    if CompSmo
        %The more conventional adjustment for the neumber of "degrees of freedom".
        noise.nu_factor = (sqrt(4*log(2)/pi)/noise.fwhm)^sum(d(1:3)>1);
    end

otherwise
    error('Unknown likelihood function.');
end


