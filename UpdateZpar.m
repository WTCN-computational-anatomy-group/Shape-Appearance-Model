function [z,S,L,omisc] = UpdateZpar(z,f,mu,Wa,Wv,A,s,noise)
% Update latent variables for one observation
% FORMAT [z,S,L,omisc] = UpdateZpar(z,f,mu,Wa,Wv,A,s,noise)
%
% z     - Vector of latent variables
% f     - This observation
% mu    - Mean
% Wa    - Appearance basis functions
% Wv    - Shape basis functions
% A     - Precision matrix of z (assumed zero mean)
% s     - Settings. Uses s.v_settings, s.nit, s.bs_args, s.vx & s.int_args
% noise - Noise precision (Gaussian model only)
%
% S     - Covariance of uncertainty of z (Laplace approximation)
% L     - Log-likelihood
% omisc - An assortment of statistics
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

verb = false;
if isfield(s,'verbose') && s.verbose>0, verb = true; end

tic
a0   = GetA0(z,Wa,mu); % Compute linear combination of appearance modes
if verb, fprintf('GetA0: %g sec\n', toc); end
tic
v0   = GetV0(z,Wv);    % Linear combination of shape components
if verb, fprintf('GetV0: %g sec\n', toc); end

subj = struct('f',f,'z',z,'a0',a0,'v0',v0,'iphi',[],'f1',[],'rho',[],'g',[],'H1',[],'nll',[],'S',[],'dz',[],'da',[],'dv',[],'stop',false);

% Approximate memory requirements per subject: f+f1+a0+da, v0+iphi+dv+rho, B
% 4*(prod(d)*4 + prod(d(1:3))*(3*3+1))
% 4*prod(d(1:3))*(d(4)*4 + 3*3+1) + 4*prod(d(1:2))*K

tic
parfor n=1:numel(subj)
    if size(Wv,5)>0  % Includes a shape model
        subj(n).iphi = GetIPhi(subj(n).v0,s);           % Generate diffeomorphism
        [f1,rho]     = Push(subj(n).f,subj(n).iphi);    % Use the resampling's adjoint operator
        subj(n).f1   = f1;
        subj(n).rho  = rho;

    else             % No shape model
        % Generate f1 and rho without warping the data
        msk = all(isfinite(subj(n).f),4);
        f1  = zeros(size(subj(n).f),'single');
        for i=1:size(f1,4)
            tmp = subj(n).f(:,:,:,i);
            tmp(~msk) = 0;
            f1(:,:,:,i) = tmp;
        end
        subj(n).f1      = f1;
        subj(n).rho     = single(msk);
        subj(n).iphi    = [];
    end

end
if verb, fprintf('Generate %d warped images: %g sec\n', numel(subj), toc); end


A   = double(A); % Use double
K   = size(A,1);

parfor n=1:numel(subj)
    subj(n).nll = -ComputeLL(subj(n).f,subj(n).iphi,subj(n).a0,s,noise) + 0.5*subj(n).z'*A*subj(n).z; % Initial objective function
end


for subit=1:s.nit

    for n=1:numel(subj)
        if ~subj(n).stop
            subj(n).g  = 0; % Gradients
            subj(n).H1 = 0; % Hessian
        end
    end

    tim = [0 0];
    for x3=1:size(Wv,3)
        % Work slice by slice to save memory
        tic;
        wa = Wa(:,:,x3,:,:);
        wv = Wv(:,:,x3,:,:);
        tim(1) = tim(1) + toc;

        tic;
        % Passing lots of data (broadcast variables) within the parfor may be slowing things down
        for n=1:numel(subj)
            if ~subj(n).stop
                [dg,dh]    = SliceComp(subj(n).z,subj(n).f1(:,:,x3,:),subj(n).rho(:,:,x3,:),subj(n).a0,x3,wa,wv,noise,s);
                subj(n).g  = subj(n).g + dg;
                subj(n).H1 = subj(n).H1 + dh;
            end
        end
        tim(2) = tim(2)+toc;
    end
    if verb, fprintf('%g sec & %g sec\n', tim); end

    for n=1:numel(subj)
        if ~subj(n).stop
            subj(n).H1  = (subj(n).H1+subj(n).H1')/2; % Ensure totally symmetric (rounding errors)

            % This should not be needed, but can be an extra layer of protection if the shooting explodes
            if any(~isfinite(subj(n).g(:))) || any(~isfinite(subj(n).H1(:)))
                subj(n).g(~isfinite(subj(n).g))   = 0;
                subj(n).H1(~isfinite(subj(n).H1)) = 0;
            end

            g  = subj(n).g  + double(A*subj(n).z); % Add prior term of gradient
            H  = subj(n).H1 + double(A);   % Add prior term of hessian

            R  = (max(diag(H))*1e-7)*eye(K); % Regularisation done in case H is singular
            H  = H+R;

            subj(n).S  = inv(H); % S encodes the uncertainty of the z estimates (Laplace approximation)
            subj(n).dz = H\g;    % Search direction for updating z
        end
    end

    tic
    da0    = GetA0({subj(:).dz},Wa);
    if verb, fprintf('GetdA0: %g sec\n', toc); end
    tic
    dv0    = GetV0({subj(:).dz},Wv);
    if verb, fprintf('GetdV0: %g sec\n', toc); end

    tic
    parfor n=1:numel(subj)
        if ~subj(n).stop
            % A Gauss-Newton update can sometimes overshoot, so a backtracking (Armijo) linesearch
            % is used to ensure the objective function improves.
            [subj(n).z,nll,misc] = LineSearch(subj(n).z,subj(n).nll,subj(n).dz,subj(n).a0,da0{n},subj(n).v0,dv0{n},subj(n).f,noise,A,s);

            % If there's an a0 field, then the linesearch has found a better solution, so save some of
            % the stuff that has been computed, ready for the next iteration.
            if isfield(misc,'a0')
                subj(n).a0   = misc.a0;
                subj(n).v0   = misc.v0;
                subj(n).iphi = misc.iphi;
                [f1,rho]     = Push(subj(n).f,subj(n).iphi);
                subj(n).f1   = f1;
                subj(n).rho  = rho;
                if abs(nll-subj(n).nll)<1e-3
                    subj(n).stop = true;
                end
                subj(n).nll = nll;
            else
                subj(n).stop = true;
            end
        end
    end
    if verb, fprintf('Linesearches: %g sec\n', toc); end
    if all(cat(1,subj.stop)), break; end
end

z = {subj.z};
S = {subj.S};

% Update the log-likelihood to include the 1/sqrt((2*pi)^d det|Sigma|) part of the Gaussian
% distribution of the prior, as well as a Laplace approximation for p(f|M) = \int_z p(f,z|M) dz
% See Section 4.4 of Bishop's book, where the Laplace approximation is described for
% use in model comparison.
L = 0;
for n=1:numel(subj)
    L = L + (-subj(n).nll + 0.5*LogDet(A) - 0.5*K*log(2*pi)...
                         -(-0.5*LogDet(subj(n).S) - 0.5*K*log(2*pi)));
end

if nargout>=4
    % If necessary, compute some extra stuff that is used for the learning part of the model

    omisc = struct('s0',0,'s1',0,'SmoSuf',0,'gmu',single(0),'Hmu',single(0));

    for n=1:numel(subj)
        % Compute sufficient statistics used for estimating the smoothness of the residuals.
        % These are used for adjusting the number of observations via a bit of random field theory.
        [SmoSuf,s0,s1] = ComputeSmoSuf(subj(n).f,subj(n).a0,subj(n).iphi,s);

        % This was originally an adjustment to the sufficient statistics for computing the noise (Gaussian npoise model).
        % The idea was to make a Bayesian estimate that accounts for the uncertainty with which
        % z is estimated.  In practice though, the effect is tiny compared to the lack of independence
        % of the neighbouring voxels.
        if false % Fix later
        switch lower(s.likelihood)
        case {'normal'}
            % Adjustment for uncertainty in z
            s1    = s1 + trace(H\H1)/noise.lam;
        end
        end

        omisc.s1     = omisc.s1 + s1;
        omisc.s0     = omisc.s0 + s0;
        omisc.SmoSuf = omisc.SmoSuf + SmoSuf;

        % Sufficient statistics that are used for recomouting the mean
        [gmu,Hmu]    = AppearanceDerivs(subj(n).f1,subj(n).rho,subj(n).a0,noise,s);
        omisc.gmu    = omisc.gmu + gmu;
        omisc.Hmu    = omisc.Hmu + Hmu;
    end
end



%==========================================================================
%
%==========================================================================
function [g,H1] = SliceComp(z,f1,rho,a0,x3,wa,wv,noise,s)

g   = 0;
H1  = 0;

d   = [size(a0) 1 1]; % Ensure dimensions are a 1x4 vector
d   = d(1:4);
Ka  = size(wa,5);    % Number of appearance modes
Kv  = size(wv,5);    % Number of shape modes

if numel(z)==Ka+Kv
    % Shape and appearance are combined and controlled together
    K      = Ka+Kv;
    Koff   = Ka;
else
    % Shape and appearance controlled by separate z elements
    K      = Ka;
    Koff   = 0;
end;


% Basis functions encoding appearance modes, or derivatives w.r.t. warps
% This encodes derivatives of f1 w.r.t z (ignoring scaling by Jacobians) - df1/dz
B   = zeros(prod([d(1:2) d(4:end)]),K);

% Fill in appearance modes
B(:,1:size(wa,5)) = reshape(wa,[prod(d([1,2,4])) size(wa,5)]);


% Fill in derivatives w.r.t. warps
tmp = zeros([d(1:2) 1 d(4:end) Kv],'single');
for l=1:d(4)
    g1 = (a0([2:end, 1],:,x3,l) - a0([end, 1:(end-1)],:,x3,l))/2;
    g2 = (a0(:,[2:end, 1],x3,l) - a0(:,[end, 1:(end-1)],x3,l))/2;
    g3 = (a0(:,:,rem(x3+1+d(3)-1,d(3))+1,l) - a0(:,:,rem(x3-1+d(3)-1,d(3))+1,l))/2;
    tmp(:,:,1,l,:) = bsxfun(@times,wv(:,:,1,1,:),-g1)...
                   + bsxfun(@times,wv(:,:,1,2,:),-g2)...
                   + bsxfun(@times,wv(:,:,1,3,:),-g3);
end
B(:,(1:Kv)+Koff) = B(:,(1:Kv)+Koff) + reshape(tmp,[prod(d([1 2 4])) Kv]);

%[a,Ha] = AppearanceDerivs(f1(:,:,x3,:),rho(:,:,x3,:),a0(:,:,x3,:),noise,s); % Compute dE/df1
[a,Ha] = AppearanceDerivs(f1,rho,a0(:,:,x3,:),noise,s); % Compute dE/df1
B      = reshape(B,[prod(d([1,2,4])),K]); % Reshape for easier matrix-vector multiplication
g      = g  + double(B'*a(:));            % dE/dz = dE/df1 * df1/dz

% Compute Hessian
switch lower(s.likelihood)
case {'multinomial','categorical'}
    % Computations are slightly different for the multinomial (categorical) noise model
    ind = Horder(d(4)); % Not all fields of Ha are saved, as it is a field of symmetric matrices (c.f. standard interview question about 1+2+3+4+...)
    B   = reshape(B,[d(1)*d(2),d(4),K]);
%% Original code didn't exploit matrix symmetry
%   for l1=1:d(4)
%       B1lt = reshape(B(:,l1,:),[d(1)*d(2),K])';
%       for l2=1:d(4)
%           B2l = reshape(B(:,l2,:),[d(1)*d(2),K]);
%           hal = Ha(:,:,1,ind(l1,l2));
%           H1  = H1 + double(B1lt*bsxfun(@times,hal(:),B2l));
%       end
%   end
    for l1=1:d(4)
        B2l  = reshape(B(:,l1,:),[d(1)*d(2),K]);
        B1lt = B2l';
        hal  = Ha(:,:,1,ind(l1,l1));
        H1   = H1 + double(B1lt*bsxfun(@times,hal(:),B2l));
        for l2=(l1+1):d(4)
            B2l = reshape(B(:,l2,:),[d(1)*d(2),K]);
            hal = Ha(:,:,1,ind(l1,l2));
            H1  = H1 + 2*double(B1lt*bsxfun(@times,hal(:),B2l));
        end
    end

otherwise
    % Simpler computations as each component can be treated separately (ie a field of diagonal matrices)
    B   = reshape(B,[d(1)*d(2),d(4),K]);
    for l=1:d(4)
        Bl  = reshape(B(:,l,:),[d(1)*d(2),K]);
        hal = Ha(:,:,1,l);
        H1  = H1 + double(Bl'*bsxfun(@times,hal(:),Bl));
    end
end


%==========================================================================
%
%==========================================================================
function [z,nll,misc] = LineSearch(oz,onll,dz,a0o,da0,v0o,dv0,f,noise,A,s)
verb   = isfield(s,'verbose') && s.verbose;
if verb, fprintf(' %g  ', onll); end

misc   = struct; % Some intermediate computations may be returned (if a better solution is found)
armijo = 1;      % Line-search parameter is decreased until objective function improves
nsubit = 6;      % Maximum number of halvings

for subit = 1:nsubit
    z    = oz  - armijo*dz;
    if numel(da0)>0
        a0   = a0o - armijo*da0;
    else
        a0   = a0o;
    end
    if numel(dv0)>0
        v0   = v0o - armijo*dv0;
        iphi = GetIPhi(v0,s); % Could achieve slight speedup by precomputing the kernel
    else
        v0   = v0o;
        iphi = [];
    end
    nll  = -ComputeLL(f,iphi,a0,s,noise) + 0.5*z'*A*z;

    if verb, fprintf(' %g', nll); end; %ShowPic(f,iphi,a0,mu,s); end

    if nll>onll+1e-3
        armijo  = armijo*0.5;
    else
        misc.a0   = a0;
        misc.v0   = v0;
        misc.iphi = iphi;
        if verb, fprintf('\n'); end
        return;
    end
end
z   = oz;
nll = onll;
if verb, fprintf('\n'); end


%==========================================================================
%
%==========================================================================
function ShowPic(f,iphi,a0,mu,s)
%return; % This function is disabled

a1  = Resamp(a0,iphi,s.bs_args);
mu1 = Resamp(mu,iphi,s.bs_args);

msk = ~isfinite(f);
ff  = f;
switch lower(s.likelihood)
case {'multinomial','categorical'}
    sig     = SoftMax(a1);
case {'binomial','binary'}
    sig = 1./(1+exp(-a1));
otherwise % case {'normal','gaussian','laplace'}
    sig = a1;
end
ff(msk) = sig(msk);

imagesc([ColourPic(ff) ColourPic(a1,s.likelihood) ColourPic(mu1,s.likelihood)]);
axis image ij off; drawnow; %set(gca,'CLim',[0 1]); drawnow


