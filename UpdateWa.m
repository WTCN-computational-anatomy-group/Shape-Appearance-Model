function Wa = UpdateWa(Wa,ga,Ha,RegW,s)
% Update appearance basis functions
% FORMAT Wa = UpdateWa(Wa,ga,Ha,RegW,s)
%
% Wa   - Appearance basis functions
% ga   - Gradients of likelihood term
% Ha   - Hessian of likelihood term
% RegW - Regularisation (precision) of estimate
% s    - Settings. Uses s.likelihood, s.vx, s.a_settings & s.mg_its
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if isempty(ga)
    return;
end
if isfield(s,'verbose') && s.verbose == true, verb = true; else verb = false; end
if verb, fprintf('\nWa: '); end

spm_field('boundary',0);

switch lower(s.likelihood)
case {'normal','gaussian','binomial','binary'}
    for k=1:size(ga,5)
        for l=1:size(ga,4)
            g             = ga(:,:,:,l,k) + spm_field('vel2mom', Wa(:,:,:,l,k), double([s.vx RegW(k,k)*s.a_settings]));
            da            = spm_field(Ha(:,:,:,l,k), g,                         double([s.vx RegW(k,k)*s.a_settings+[0.01 0.1 0] s.mg_its]));
            Wa(:,:,:,l,k) = Wa(:,:,:,l,k) - s.omega*da;
        end
    end

case {'multinomial','categorical'}
    % Could be made more efficient in terms of disk I/O etc.  This representation
    % has some redundancy, which could be eliminated using an approach like the
    % one in spm12/toolbox/Shoot/spm_shoot_blur.m , which works with the subspace
    % null(ones(1,d(4))).  This approach may also increase stability as the Hessian
    % of the likelihood would not be singular.
    for k=1:size(ga,5)
        g             = ga(:,:,:,:,k) + spm_field('vel2mom', Wa(:,:,:,:,k), double([s.vx RegW(k,k)*s.a_settings]));
        da            = spm_field(Ha(:,:,:,:,k), g,                         double([s.vx RegW(k,k)*s.a_settings+[0.01 0.1 0] s.mg_its]));
        tmp           = Wa(:,:,:,:,k) - s.omega*da;
        tmp           = bsxfun(@minus,tmp,sum(tmp,4)/size(tmp,4));
        Wa(:,:,:,:,k) = tmp;
        if verb, fprintf(' %g', sqrt(sum(da(:).^2))); end
    end
end
if verb, fprintf('\n'); end

