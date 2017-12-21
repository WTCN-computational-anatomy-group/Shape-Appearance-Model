function [mu,noise] = ComputeMean(s0,s1,s2,s)
% Compute initial mean (and variance) from sufficient statistics
% FORMAT [mu,noise] = PGinit(s0,s1,s2,s)
%
% s0    - Zeroeth moment (number of observations)
% s1    - First moment   (sum over observations)
% s2    - Second moment  (sum of squares of observations)
% s     - Settings. Uses s.likelihood, s.vx, s.mu_settings, s.mg_its (and
%         maybe s.alpha0, s.beta0)
%
% mu    - mean
% noise - an assortment of stuff used when s.likelihood = 'normal'.
%
% This function is used as part of initialising the model fitting
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


noise = struct;
d      = [size(s1) 1 1];
d      = d(1:4);
mu     = zeros(d,'single');
N      = double(nanmax(s0(:)));

noise.nu_factor = 1;
if isfield(s,'nu_factor'), noise.nu_factor = s.nu_factor; end

spm_field('boundary',0);

for it=1:8

    switch lower(s.likelihood)
    case {'normal','gaussian'}
        alpha0 = 0.0001; if isfield(s,'alpha0'), alpha0 = s.alpha0; end
        beta0  = 0.0001; if isfield(s,'beta0' ), beta0  = s.beta0;  end
        gmu   = bsxfun(@times,s0,mu) - s1;
        hmu   = s0;
        noise.beta  = beta0  + reshape(sum(sum(sum(max(s2-bsxfun(@rdivide,s1.^2,s0),0),1),2),3),[1 size(s1,4)]);
        noise.alpha = alpha0 + sum(sum(sum(s0,1),2),3);
        noise.lam   = double(noise.alpha./noise.beta);
        for l=1:d(4)
            gmul          = gmu(:,:,:,l) + spm_field('vel2mom', mu(:,:,:,l), [s.vx N*s.mu_settings/noise.lam(l)/noise.nu_factor]);
            mu(:,:,:,l)   = mu(:,:,:,l)  - spm_field(hmu, gmul,              [s.vx N*s.mu_settings/noise.lam(l)/noise.nu_factor s.mg_its]);
        end

    case {'binomial','binary'}
        sig   = 1./(1+exp(-mu));
        gmu   = s0.*sig - s1;
        hmu   = s0.*(sig.*(1-sig)+1e-3);
        noise.lam = 1;
        gmu  = gmu + spm_field('vel2mom', mu, [s.vx N*s.mu_settings/noise.nu_factor]);
        mu   = mu  - spm_field(hmu, gmu,      [s.vx N*s.mu_settings/noise.nu_factor s.mg_its]);

    case {'multinomial','categorical'}
        sig    = SoftMax(mu);
        gmu    = bsxfun(@times,s0,sig) - s1;
        Hmu    = zeros([d(1:3) d(4)*(d(4)+1)/2],'single');
        for l = 1:d(4)
            Hmu(:,:,:,l) = s0.*(sig(:,:,:,l) - sig(:,:,:,l).^2 + 1e-3*d(4)^2);
        end
        l = d(4);
        for l1 = 1:d(4)
            for l2 = (l1+1):d(4)
                l            = l + 1;
                Hmu(:,:,:,l) = -s0.*sig(:,:,:,l1).*sig(:,:,:,l2);
            end
        end
        noise.lam = 1;
        gmu  = gmu + spm_field('vel2mom', mu, [s.vx N*s.mu_settings/noise.nu_factor]);
        mu   = mu  - spm_field(Hmu, gmu,      [s.vx N*s.mu_settings/noise.nu_factor s.mg_its]);
    otherwise
        noise.lam   = 1;
    end

   %imagesc(mu(:,:,1,1)'); axis image xy off; drawnow
end

