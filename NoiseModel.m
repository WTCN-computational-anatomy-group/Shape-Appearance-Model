function [noise,lb_lam] = NoiseModel(ss,s,d)

% See the appendix of:
% Groves AR, Beckmann CF, Smith SM, Woolrich MW. Linked independent component analysis
% for multimodal data fusion. Neuroimage. 2011 Feb 1;54(3):2198-217.

SmoSuf          = ss.SmoSuf;
noise.fwhm      = sqrt(4*log(2)*(SmoSuf(2)/SmoSuf(1))/(sum(SmoSuf([4 6 8]))/sum(SmoSuf([3 5 7]))));
noise.fwhm      = sqrt(max(noise.fwhm^2 - 0.5, 4*log(2)/pi)); % An adjustment that helps empirically

switch lower(s.likelihood)
case {'normal','gaussian','laplace'}

    % Note that a slightly different adjustment is used for the Gaussian noise model
    % This one is a heavier adjustment, which is based on taking iid data and smoothing
    % in order to estimate the original variance of the iid noise.
    noise.nu_factor = (sqrt(2*log(2)/pi)/noise.fwhm)^sum(d(1:3)>1);
    if isfield(s,'nu_factor'), noise.nu_factor = s.nu_factor; end

    %fprintf('  %g,%g',noise.fwhm, noise.nu_factor);

    noise.alpha     = s.alpha0 + ss.s0*noise.nu_factor;
    noise.beta      = s.beta0  + ss.s1*noise.nu_factor;
    noise.lam       = noise.alpha./noise.beta;

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
        lb_plam = lb_plam + s.alpha0(l)*log(s.beta0(l)) ...
                          + (s.alpha0(l)-1)*(psi(noise.alpha(l)) - log(noise.beta(l))) ...
                          -(noise.alpha(l)/noise.beta(l))*s.beta0(l) ...
                          - gammaln(s.alpha0(l));
    end
    lb_lam  = lb_lam + lb_plam - lb_qlam;

case {'binomial','binary','multinomial','categorical'}
    %The more conventional adjustment for the neumber of "degrees of freedom".
    noise.nu_factor = (sqrt(4*log(2)/pi)/noise.fwhm)^sum(d(1:3)>1);
    if isfield(s,'nu_factor'), noise.nu_factor = s.nu_factor; end

    fprintf('  %g,%g',noise.fwhm, noise.nu_factor);
    noise.lam = 1;
    lb_lam    = 0;
otherwise
    error('Unknown likelihood function.');
end


