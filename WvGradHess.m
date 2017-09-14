function [gv,Hv,nll] = WvGradHess(dat,mu,Wa,Wv,noise,s)
% Compute derivatives w.r.t. shape basis functions
% FORMAT [gv,Hv,nll] = WaGradHess(dat,mu,Wa,Wv,noise,s)
% dat   - Array of data
% mu    - Mean
% Wa    - Appearance basis functions
% Wv    - Shape basis functions
% noise - Noise information
% s     - Settings. Uses s.likelihood, s.ondisk, s.result_dir & s.result_name
%
% gv    - Gradients 
% Hv    - Hessians
% nll   - Negative log-likelihood
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if isempty(dat), gv = []; Hv = []; nll = 0; return; end

Kv = size(Wv,5);
Ka = size(Wa,5);

K  = size(dat(1).z,1);
if (Ka == Kv) && (Ka == K), Koff = 0; else Koff = Ka; end

if Kv==0, gv = zeros([Kv 0],'single'); Hv = []; nll = []; return; end

d   = [size(mu) 1 1];
d   = d(1:4);

batchsize = 1;
if isfield(s,'batchsize'), batchsize = s.batchsize; end

if isfield(s,'ondisk') && s.ondisk
    gv = file_array(fullfile(s.result_dir,[s.result_name '_gv.dat']), [d(1:3) 3 Kv],'float32',352);
    Hv = file_array(fullfile(s.result_dir,[s.result_name '_Hv.dat']), [d(1:3) 6 Kv],'float32',352);
else
    gv = zeros([d(1:3) 3 Kv],'single');
    Hv = zeros([d(1:3) 6 Kv],'single');
end

nll = 0;
% Compute 1st and 2nd derivatives w.r.t. velocities
for k=1:Kv,
    gv(:,:,:,:,k) = 0;
    Hv(:,:,:,:,k) = 0;
end

if Ka==0
    % The same gradients can be used throughout
    a0  = mu;
    Gmu = CompGrads(a0,s.bs_args);
else
    Gmu = {};
end

for n1=1:batchsize:numel(dat)
    nn    = n1:min(n1+(batchsize-1),numel(dat));
    z     = {dat(nn).z};
    S     = {dat(nn).S};
    cell1 = GetV0(z,Wv); % Replaced by gradients
    cell2 = GetA0(z,Wa,mu);  % Replaced by Hessians
    dat1  = dat(nn);

    comp_nll = nargout>=3;
    parfor n=1:numel(nn)
        iphi     = GetIPhi(cell1{n},s);
        a0       = cell2{n};
        f        = GetDat(dat1(n),s);
        if comp_nll
            nll = nll - ComputeLL(f,iphi,a0,s,noise);
        end

        [f1,rho] = Push(f,iphi);

        if Ka>0
            G   = CompGrads(a0,s.bs_args);
        else
            G   = Gmu;
        end

        g       = zeros([d(1:3),3],'single'); % First derivatives
        H       = zeros([d(1:3),6],'single'); % Second derivatives

        switch lower(s.likelihood)
        case {'normal','gaussian'}
            for l=1:d(4)
                al         = (noise.nu_factor*noise.lam(l))*(f1(:,:,:,l)-rho.*a0(:,:,:,l));
                g(:,:,:,1) = g(:,:,:,1) + al.*G{l,1};
                g(:,:,:,2) = g(:,:,:,2) + al.*G{l,2};
                g(:,:,:,3) = g(:,:,:,3) + al.*G{l,3};

                wl         = (noise.nu_factor*noise.lam(l))*rho;
                H(:,:,:,1) = H(:,:,:,1) + wl.*G{l,1}.*G{l,1};
                H(:,:,:,2) = H(:,:,:,2) + wl.*G{l,2}.*G{l,2};
                H(:,:,:,3) = H(:,:,:,3) + wl.*G{l,3}.*G{l,3};
                H(:,:,:,4) = H(:,:,:,4) + wl.*G{l,1}.*G{l,2};
                H(:,:,:,5) = H(:,:,:,5) + wl.*G{l,1}.*G{l,3};
                H(:,:,:,6) = H(:,:,:,6) + wl.*G{l,2}.*G{l,3};
            end

        case {'laplace'}
           %b   = reshape(sqrt(1./(2*noise.lam)),[1,1,1,d(4)]);
            for l=1:d(4)
%               r          = (f1(:,:,:,l)./(rho+eps) - a0(:,:,:,l))/b(l);
                r          = (f1(:,:,:,l)./(rho+eps) - a0(:,:,:,l));
                q          = 2./sqrt((noise.lam(l)/2)*r.^2+0.01);
                wl         = noise.nu_factor*noise.lam(l)*q.*rho;
                al         = wl.*(f1(:,:,:,l)./(rho+eps) - a0(:,:,:,l));
%               wt         = rho./max(abs(r),0.0001);
%               al         = r.*wt;
                g(:,:,:,1) = g(:,:,:,1) + al.*G{l,1};
                g(:,:,:,2) = g(:,:,:,2) + al.*G{l,2};
                g(:,:,:,3) = g(:,:,:,3) + al.*G{l,3};
%               wl         = wt/b(l);
                H(:,:,:,1) = H(:,:,:,1) + wl.*G{l,1}.*G{l,1};
                H(:,:,:,2) = H(:,:,:,2) + wl.*G{l,2}.*G{l,2};
                H(:,:,:,3) = H(:,:,:,3) + wl.*G{l,3}.*G{l,3};
                H(:,:,:,4) = H(:,:,:,4) + wl.*G{l,1}.*G{l,2};
                H(:,:,:,5) = H(:,:,:,5) + wl.*G{l,1}.*G{l,3};
                H(:,:,:,6) = H(:,:,:,6) + wl.*G{l,2}.*G{l,3};
            end

        case {'binomial','binary'}
            ea         = exp(a0);
            sig        = ea./(1+ea);
            a          = noise.nu_factor*(f1-sig.*rho);
            g(:,:,:,1) = a.*G{1,1};
            g(:,:,:,2) = a.*G{1,2};
            g(:,:,:,3) = a.*G{1,3};
            wt         = noise.nu_factor*rho.*(sig.*(1-sig)+1e-3);
            H(:,:,:,1) = wt.*G{1,1}.*G{1,1};
            H(:,:,:,2) = wt.*G{1,2}.*G{1,2};
            H(:,:,:,3) = wt.*G{1,3}.*G{1,3};
            H(:,:,:,4) = wt.*G{1,1}.*G{1,2};
            H(:,:,:,5) = wt.*G{1,1}.*G{1,3};
            H(:,:,:,6) = wt.*G{1,2}.*G{1,3};

        case {'multinomial','categorical'}
            sig    = SoftMax(a0);
            g      = zeros([d(1:3) 3],'single');
            for l=1:d(4)
                a          = noise.nu_factor*(f1(:,:,:,l) - sig(:,:,:,l).*rho);
                g(:,:,:,1) = g(:,:,:,1) + a.*G{l,1};
                g(:,:,:,2) = g(:,:,:,2) + a.*G{l,2};
                g(:,:,:,3) = g(:,:,:,3) + a.*G{l,3};
            end

            H      = zeros([d(1:3) 6],'single');
            for l1 = 1:d(4)
                for l2 = 1:d(4)
                    wt = -rho.*sig(:,:,:,l1).*sig(:,:,:,l2);
                    if l1==l2, wt = wt + rho.*(sig(:,:,:,l1)+1e-3); end
                    wt = noise.nu_factor*wt;
                    H(:,:,:,1) = H(:,:,:,1) + wt.*G{l1,1}.*G{l2,1};
                    H(:,:,:,2) = H(:,:,:,2) + wt.*G{l1,2}.*G{l2,2};
                    H(:,:,:,3) = H(:,:,:,3) + wt.*G{l1,3}.*G{l2,3};
                    H(:,:,:,4) = H(:,:,:,4) + wt.*G{l1,1}.*G{l2,2};
                    H(:,:,:,5) = H(:,:,:,5) + wt.*G{l1,1}.*G{l2,3};
                    H(:,:,:,6) = H(:,:,:,6) + wt.*G{l1,2}.*G{l2,3};
                end
            end
        otherwise
            error('Unknown likelihood function.');
        end

        g(~isfinite(g)) = 0;
        H(~isfinite(H)) = 0;

        cell1{n} = g;
        cell2{n} = H;
    end

    % Add appropriate amount of gradient and Hessian
    for k=1:Kv,
        g1 = single(0);
        H1 = single(0);
        for n=1:numel(cell1)
            g1  = g1 + cell1{n}*z{n}(k+Koff);
            ezz = z{n}(k+Koff)^2;% + S{n}(k+Koff,k+Koff);
            H1  = H1 + cell2{n}*ezz;
        end
        gv(:,:,:,:,k) = gv(:,:,:,:,k) + g1;
        Hv(:,:,:,:,k) = Hv(:,:,:,:,k) + H1;
    end
end

