function [mu,lb_pmu] = UpdateMu(mu,gmu,Hmu,N,s)
% Update the mean
% FORMAT [mu,lb_pmu] = UpdateMu(mu,gmu,Hmu,N,s)
%
% mu     - The mean
% gmu    - 1st derivatives
% Hmu    - 2nd derivatives
% N      - Number of observations
% s      - Settings. Uses s.likelihood, s.vx, s.mu_settings & s.mg_its.
% lb_pmu - Part of lower bound
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

d = [size(mu) 1 1];
d = d(1:4);

switch lower(s.likelihood)
case {'normal','gaussian','laplace','binomial','binary'}
    a0 = zeros(d,'single');
    for l=1:d(4)
        gr          = gmu(:,:,:,l) + spm_field('vel2mom', mu(:,:,:,l), [s.vx N*s.mu_settings]);
        dmu         = spm_field(Hmu(:,:,:,l), gr,                      [s.vx N*s.mu_settings+[0.01 0.01 0]*0 s.mg_its]);
        mu(:,:,:,l) = mu(:,:,:,l) - dmu;
        a0(:,:,:,l) = spm_field('vel2mom', mu(:,:,:,l),                [s.vx N*s.mu_settings]);
    end

case {'multinomial','categorical'}
    % Could be made more efficient in terms of disk I/O etc.  This representation
    % has some redundancy, which could be eliminated using an approach like the
    % one in spm12/toolbox/Shoot/spm_shoot_blur.m , which works with the subspace
    % null(ones(1,d(4))).  This approach may also increase stability as the Hessian
    % of the likelihood would not be singular.
    gr     = gmu + spm_field('vel2mom', mu, [s.vx N*s.mu_settings]);
    dmu    = spm_field(Hmu, gr,             [s.vx N*s.mu_settings+[0.01 0.01 0]*0 s.mg_its]);
    mu     = mu - dmu;
    mu     = bsxfun(@minus,mu,sum(mu,4)/size(mu,4));
    a0     = spm_field('vel2mom', mu,       [s.vx N*s.mu_settings]);
end
lb_pmu = -0.5*sum(sum(sum(sum(sum(a0.*mu)))));

