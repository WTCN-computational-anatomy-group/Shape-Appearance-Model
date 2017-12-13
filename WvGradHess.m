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
    Gmu = CompGrads(a0);
else
    Gmu = {};
end

for n1=1:batchsize:numel(dat)
    nn    = n1:min(n1+(batchsize-1),numel(dat));
    z     = {dat(nn).z};
   %S     = {dat(nn).S};
    cell1 = GetV0(z,Wv);    % Replaced by gradients
    cell2 = GetA0(z,Wa,mu); % Replaced by Hessians
    dat1  = dat(nn);

    parfor n=1:numel(nn)
        iphi     = GetIPhi(cell1{n},s);
        a0       = cell2{n};
        f        = GetDat(dat1(n),s);
        [ll,a,Ha] = AppearanceDerivs(f,a0,iphi,noise,s);
        nll = nll - ll;

        if Ka>0
            G   = CompGrads(a0);
        else
            G   = Gmu;
        end

        g       = zeros([d(1:3),3],'single'); % First derivatives
        H       = zeros([d(1:3),6],'single'); % Second derivatives

        switch lower(s.likelihood)
        case {'normal','gaussian'}
            for l=1:d(4)
                al         =-a(:,:,:,l);
                g(:,:,:,1) = g(:,:,:,1) + al.*G{l,1};
                g(:,:,:,2) = g(:,:,:,2) + al.*G{l,2};
                g(:,:,:,3) = g(:,:,:,3) + al.*G{l,3};

                wl         = Ha(:,:,:,l);
                H(:,:,:,1) = H(:,:,:,1) + wl.*G{l,1}.*G{l,1};
                H(:,:,:,2) = H(:,:,:,2) + wl.*G{l,2}.*G{l,2};
                H(:,:,:,3) = H(:,:,:,3) + wl.*G{l,3}.*G{l,3};
                H(:,:,:,4) = H(:,:,:,4) + wl.*G{l,1}.*G{l,2};
                H(:,:,:,5) = H(:,:,:,5) + wl.*G{l,1}.*G{l,3};
                H(:,:,:,6) = H(:,:,:,6) + wl.*G{l,2}.*G{l,3};
            end

        case {'binomial','binary'}
            g(:,:,:,1) =-a.*G{1,1};
            g(:,:,:,2) =-a.*G{1,2};
            g(:,:,:,3) =-a.*G{1,3};
            wt         = Ha+1e-4;
            H(:,:,:,1) = wt.*G{1,1}.*G{1,1};
            H(:,:,:,2) = wt.*G{1,2}.*G{1,2};
            H(:,:,:,3) = wt.*G{1,3}.*G{1,3};
            H(:,:,:,4) = wt.*G{1,1}.*G{1,2};
            H(:,:,:,5) = wt.*G{1,1}.*G{1,3};
            H(:,:,:,6) = wt.*G{1,2}.*G{1,3};

        case {'multinomial','categorical'}
            g      = zeros([d(1:3) 3],'single');
            for l=1:d(4)
                al         =-a(:,:,:,l);
                g(:,:,:,1) = g(:,:,:,1) + al.*G{l,1};
                g(:,:,:,2) = g(:,:,:,2) + al.*G{l,2};
                g(:,:,:,3) = g(:,:,:,3) + al.*G{l,3};
            end

            H   = zeros([d(1:3) 6],'single');
            ind = Horder(d(4));
            for l1 = 1:d(4)
                for l2 = l1:d(4)
                    wt = Ha(:,:,:,ind(l1,l2));
                    if l1~=l2, wt = wt*2; end
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

