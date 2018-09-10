function Wv = UpdateWv(Wv,gv,Hv,RegW,s)
% Update shape basis functions
% FORMAT WWv = UpdateWv(Wv,gv,Hv,RegW,s)
%
% Wv   - Shape basis functions
% gv   - Gradients of likelihood term
% Hv   - Hessian of likelihood term
% RegW - Regularisation (precision) of estimate
% s    - Settings. Uses s.vx, s.v_settings & s.mg_its
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if isempty(gv)
    return;
end
if isfield(s,'verbose') && s.verbose == true
    verb = true;
else
    verb = false;
end
if verb, fprintf('\nWv: '); end
for k=1:size(gv,5)
    spm_diffeo('boundary',0);
    g1            = spm_diffeo('vel2mom', Wv(:,:,:,:,k), [s.vx RegW(k,k)*s.v_settings]);
    g             = gv(:,:,:,:,k) + g1;
    dv0           = spm_diffeo('fmg', Hv(:,:,:,:,k), g,  [s.vx RegW(k,k)*s.v_settings s.mg_its]);
    Wv(:,:,:,:,k) = Wv(:,:,:,:,k) - s.omega*dv0;
    if verb, fprintf(' %g', sqrt(sum(dv0(:).^2))); end
end
if verb, fprintf('\n'); end

