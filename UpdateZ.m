function [z,S,L,omisc] = UpdateZ(z,f,mu,Wa,Wv,A,s,noise)
% Update latent variables for one observation
% FORMAT [z,S,L,omisc] = UpdateZ(z,f,mu,Wa,Wv,A,s,noise)
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

z   = z(:);          % Ensure z is a column vector
d   = [size(f) 1 1]; % Ensure dimensions are a 1x4 vector
d   = d(1:4);
Ka  = size(Wa,5);    % Number of appearance modes
Kv  = size(Wv,5);    % Number of shape modes

if numel(z)==Ka+Kv
    % Shape and appearance are combined and controlled together
    K      = Ka+Kv;
    Koff   = Ka;
else
    % Shape and appearance controlled by separate z elements
    K      = Ka;
    Koff   = 0;
end;

a0      = GetA0(z,Wa,mu); % Compute linear combination of appearance modes

if size(Wv,5)>0  % Includes a shape model
    v0       = GetV0(z,Wv);             % Linear combination of shape components
    iphi     = GetIPhi(v0,s);           % Generate diffeomorphism
    [f1,rho] = Push(f,iphi);            % Use the resampling's adjoint operator
    G        = CompGrads(a0,s.bs_args); % Compute spatial gradients of appearance model

else             % No shape model
    % Generate f1 and rho without warping the data
    msk     = all(isfinite(f),4);
    f1      = zeros(size(f),'single');
    for i=1:size(f1,4)
        tmp = f(:,:,:,i);
        tmp(~msk) = 0;
        f1(:,:,:,i) = tmp;
    end 
    rho     = single(msk);
    v0      = [];
    iphi    = [];
end

A   = double(A); % Use double
nll = -ComputeLL(f,iphi,a0,s,noise) + 0.5*z'*A*z; % Initial objective function

for subit=1:s.nit

    g  = zeros(K,1);   % Gradients

    if subit==1
        H1 = zeros(K); % Hessian
    end

    for x3=1:d(3)
        % Work slice by slice to save memory

        % Basis functions encoding appearance modes, or derivatives w.r.t. warps
        % This encodes derivatives of f1 w.r.t z (ignoring scaling by Jacobians) - df1/dz
        B   = zeros(prod([d(1:2) d(4:end)]),K);

        % Fill in appearance modes
        w   = Wa(:,:,x3,:,:);
        B(:,1:size(Wa,5)) = reshape(w,[prod(d([1,2,4])) size(Wa,5)]);

        % Fill in derivatives w.r.t. warps
        w   = Wv(:,:,x3,:,:);
        tmp = zeros([d(1:2) 1 d(4:end) Kv],'single');
        for l=1:d(4)
            tmp(:,:,1,l,:) = bsxfun(@times,w(:,:,1,1,:),-G{l,1}(:,:,x3))...
                           + bsxfun(@times,w(:,:,1,2,:),-G{l,2}(:,:,x3))...
                           + bsxfun(@times,w(:,:,1,3,:),-G{l,3}(:,:,x3));
        end
        B(:,(1:Kv)+Koff) = B(:,(1:Kv)+Koff) + reshape(tmp,[prod(d([1 2 4])) Kv]);

        if subit==1
            [a,Ha] = AppearanceDerivs(f1(:,:,x3,:),rho(:,:,x3,:),a0(:,:,x3,:),noise,s); % Compute dE/df1
            B      = reshape(B,[prod(d([1,2,4])),K]); % Reshape fo easier matrix-vector multiplication
            g      = g  + double(B'*a(:));            % dE/dz = dE/df1 * df1/dz

            % Compute Hessian
            switch lower(s.likelihood)
            case {'multinomial','categorical'}
                % Computations are slightly different for the multinomial (categorical) noise model
                ind = Horder(d(4)); % Not all fields of Ha are saved, as it is a field of symmetric matrices (c.f. standard interview question about 1+2+3+4+...)
                B   = reshape(B,[d(1)*d(2),d(4),K]);
%               for l1=1:d(4)
%                   B1lt = reshape(B(:,l1,:),[d(1)*d(2),K])';
%                   for l2=1:d(4)
%                       B2l = reshape(B(:,l2,:),[d(1)*d(2),K]);
%                       hal = Ha(:,:,1,ind(l1,l2));
%                       H1  = H1 + double(B1lt*bsxfun(@times,hal(:),B2l));
%                   end
%               end
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
        else
            a      = AppearanceDerivs(f1(:,:,x3,:),rho(:,:,x3,:),a0(:,:,x3,:),noise,s); % Compute dE/df1
            B      = reshape(B,[prod(d([1,2,4])),K]); % Reshape fo easier matrix-vector multiplication
            g      = g  + double(B'*a(:));            % dE/dz = dE/df1 * df1/dz
        end
    end

    if subit==1
        H1  = (H1+H1')/2; % Ensure totally symmetric (rounding errors)
    end

    % This should not be needed, but can be an extra layer of protection if the shooting explodes
    if any(~isfinite(g(:))) || any(~isfinite(H1(:)))
        g(~isfinite(g))   = 0;
        H1(~isfinite(H1)) = 0;
    end

    g      = g  + double(A*z); % Add prior term of gradient
    H      = H1 + double(A);   % Add prior term of hessian

    R      = (max(diag(H))*1e-7)*eye(K); % Regularisation done in case H is singular
    H      = H+R;

    S      = inv(H); % S encodes the uncertainty of the z estimates (Laplace approximation)
    dz     = H\g;    % Search direction for updating z
    nllo   = nll;    % Previous (negative) log-likelihood

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=S*0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A Gauss-Newton update can sometimes overshoot, so a backtracking (Armijo) linesearch
    % is used to ensure the objective function improves.
    [z,nll,misc] = LineSearch(z,nllo,dz,a0,v0,f,mu,Wa,Wv,noise,A,s);

    % If there's an a0 field, then the linesearch has found a better solution, so save some of
    % the stuff that has been computed, ready for the next iteration.
    if isfield(misc,'a0')
        a0 = misc.a0;
        G  = CompGrads(a0,s.bs_args);
        if isfield(misc,'iphi')
            v0       = misc.v0;
            iphi     = misc.iphi;
            [f1,rho] = Push(f,iphi);
        end
        if abs(nll-nllo)<1e-4, break; end;
    else
        break
    end
end

% Update the log-likelihood to include the 1/sqrt((2*pi)^d det|Sigma|) part of the Gaussian
% distribution of the prior, as well as a Laplace approximation for p(f|M) = \int_z p(f,z|M) dz
% See Section 4.4 of Bishop's book, where the Laplace approximation is described for
% use in model comparison.
L = -nll + 0.5*LogDet(A) - 0.5*K*log(2*pi)...
         -(0.5*LogDet(H) - 0.5*K*log(2*pi));

if nargout>=4
    % If necessary, compute some extra stuff that is used for the learning part of the model

    % Compute sufficient statistics used for estimating the smoothness of the residuals.
    % These are used for adjusting the number of observations via a bit of random field theory.
    [SmoSuf,s0,s1] = ComputeSmoSuf(f,a0,iphi,s);

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

    omisc.s1  = s1;
    omisc.s0  = s0;
    omisc.SmoSuf = SmoSuf;

    % Sufficient statistics that are used for recomouting the mean
    [omisc.gmu,omisc.Hmu] = AppearanceDerivs(f1,rho,a0,noise,s);
end

ShowPic(f,iphi,a0,mu,s);


%==========================================================================
%
%==========================================================================
function [z,nll,misc] = LineSearch(oz,onll,dz,a0o,v0o,f,mu,Wa,Wv,noise,A,s)
verb   = isfield(s,'verbose') && s.verbose;
if verb, fprintf(' %g  ', onll); end

misc   = struct; % Some intermediate computations may be returned (if a better solution is found)
armijo = 1;      % Line-search parameter is decreased until objective function improves
nsubit = 6;      % Maximum number of halvings

% Instead of recomputing the linear combination of basis functions at each step, it is
% more efficient (in terms of I/O) to make an initial linear combination, and rescale it.
da0    = GetA0(dz,Wa);
dv0    = GetV0(dz,Wv);

for subit = 1:nsubit
    z    = oz  - armijo*dz;
    if size(Wa,5)>0
        a0   = a0o - armijo*da0;
    else
        a0   = a0o;
    end
    if size(Wv,5)>0
        v0   = v0o - armijo*dv0;
        iphi = GetIPhi(v0,s);
    else
        v0   = v0o;
        iphi = [];
    end
    nll  = -ComputeLL(f,iphi,a0,s,noise) + 0.5*z'*A*z;

    if verb, fprintf(' %g', nll); ShowPic(f,iphi,a0,mu,s); end

    if nll>onll+1e-4
        armijo  = armijo*0.5;
    else
        misc.a0 = a0;
        misc.v0 = v0;
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


