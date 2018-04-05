function PG(dat,s)
% Combined principal geodesic analysis and generalized(ish) PCA
% FORMAT PG(dat,s)
% dat - a data structure containing filenames etc
% s   - various settings
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

spm_field('boundary',0);
spm_diffeo('boundary',0);

if nargin<2
    s = dat;
end

if ~isfield(s,'continue') || s.continue==false,
    if ~exist('s','var'), error('No settings'); end

    PGdistribute('init',s);
    PGdistribute('share',dat);

    % Start from scratch
    [s0,s1,s2,mat] = PGdistribute('SuffStats',s);
    s.vx           = sqrt(sum(mat(1:3,1:3).^2));
    [mu,noise] = ComputeMean(s0,s1,s2,s);
    d          = [size(mu) 1 1];
    [mu_fa,Wa,Wv,WWa,WWv,WW] = CreateBases(s,mu,mat);
    K          = size(WW,1);

    PGdistribute('RandomZ',K);
    [ss.N,ss.Z,ss.ZZ,ss.S] = PGdistribute('GetZZ');
    PGdistribute('AddToZ',-ss.Z/ss.N);
    ss.ZZ     = ss.ZZ - ss.Z*ss.Z'/ss.N;
    ss.Z      = ss.Z*0;
    [U,S]     = svd(ss.ZZ);
    Rz        = sqrt(ss.N)*U/sqrtm(S);
   %Rz        = U/sqrtm(S);

    ss.ZZ     = Rz'*ss.ZZ*Rz;
    ss.S      = Rz'*ss.S *Rz;
    PGdistribute('TransfZ',Rz');
    [ss.N,ss.Z,ss.ZZ,ss.S] = PGdistribute('GetZZ');
    s.omega   = 1;
else
    % Continue from previous results
    new_s = s;
    load(fullfile(s.result_dir,['train' s.result_name '.mat']),...
        'Wa','Wv','WWa','WWv','mat','dat','ss','A','B','mu','s','noise','dat');
    d     = [size(mu) 1 1];
    if isfield(new_s,'nu_factor')
        s.nu_factor     = new_s.nu_factor;
        noise.nu_factor = s.nu_factor;
    end
    if isfield(new_s,'v_settings')
        s.v_settings  = new_s.v_settings;
    end
    if isfield(new_s,'a_settings')
        s.a_settings  = new_s.a_settings;
    end
    if isfield(new_s,'mu_settings')
        s.mu_settings = new_s.mu_settings;
    end

    % Should really include some checks here
    PGdistribute('init',s);
    [ss.N,ss.Z,ss.ZZ,ss.S] = PGdistribute('GetZZ');
    Cv    = eye(size(ss.ZZ))*max(diag(ss.ZZ))/ss.N*0.1; % Break symmetry if necessary
    PGdistribute('AddRandZ',Cv);
    [ss.N,ss.Z,ss.ZZ,ss.S] = PGdistribute('GetZZ');

    WWa   = UpdateWWa(Wa,s);
    WWv   = UpdateWWv(Wv,s);
    [Wa,Wv,WWa,WWv,ss,WW] = OrthAll(Wa,Wv,WWa,WWv,ss,s);
end

if ~isfield(s,'lambda'), s.lambda = [1 1];    end
maxit = 15;    if isfield(s,'maxit'), maxit = s.maxit; end

[A,B,lb_qA,lb_pA,ldA] = SetReg(ss.ZZ+ss.S,ss.N,s);
lb_A                  = lb_pA  - lb_qA;


if isfield(s,'debug') && s.debug
    subplot(2,2,1); image(ColourPic(mu,s.likelihood)); axis image ij off;
    drawnow
end

for iter = 1:maxit,
    fprintf('%-3d    ', iter);

    [Wa,Wv,WW,s.omega] = Mstep(mu,Wa,Wv,noise,B,ss.ZZ,A,WWa,WWv,s);

    lb_pW = -0.5*trace(B*WW); % + const

    RegZ  = double(s.lambda(1)*A + s.lambda(2)*WW);
    ss    = PGdistribute('UpdateLatentVariables',mu,Wa,Wv,noise,RegZ,s);
    lb_L  = [ss.L(1) ss.L(3)];
    noise = NoiseModel(ss,s,d);

    % Zero-mean Z and recompute the mean mu
    % fprintf(' (%g,%g) ', norm(ss.Z/ss.N),sqrt(trace(ss.ZZ)/ss.N));
    PGdistribute('AddToZ',-ss.Z/ss.N);
    ss.ZZ             = ss.ZZ - ss.Z*ss.Z'/ss.N;
    ss.Z              = ss.Z*0;

    [ss.gmu,ss.Hmu,~] = PGdistribute('MeanDerivatives',mu,Wa,Wv,noise,s);
    [mu,lb_pmu]       = UpdateMu(mu,ss.gmu,ss.Hmu,ss.N,s);
    if isfield(s,'ondisk') && s.ondisk
        mu_fa(:) = mu(:);
    end

    if isfield(s,'debug') && s.debug
        subplot(2,2,1); image(ColourPic(mu,s.likelihood)); axis image ij off;
        subplot(2,2,2); imagesc(abs(ss.ZZ)/ss.N);  colorbar; axis image; title ZZ
        subplot(2,2,3); imagesc(abs(WW));  colorbar; axis image; title WW
        subplot(2,2,4); imagesc(abs(WWv)); colorbar; axis image; title WWv
        drawnow
    end

    lb = lb_L(1) + s.lambda(1)*(lb_A + 0.5*ss.N*ldA) + s.lambda(1)*lb_pW + lb_pmu + noise.lb_lam;
    fprintf('%9.6g %7.4g    %9.6g %9.6g    %9.6g %9.6g\n',...
            (ss.L(1)+s.lambda(1)*lb_pW), s.omega, s.lambda(1)*(lb_A + 0.5*ss.N*ldA), s.lambda(1)*lb_pW + lb_pmu + noise.lb_lam, lb, lb+lb_L(2));

    WWa     = UpdateWWa(Wa,s);
    WWv     = UpdateWWv(Wv,s);
    [Wa,Wv,WWa,WWv,ss,WW] = OrthAll(Wa,Wv,WWa,WWv,ss,s);
    [A,B,lb_qA,lb_pA,ldA] = SetReg(ss.ZZ-ss.Z'*ss.Z/ss.N+ss.S,ss.N,s);
    lb_A                  = lb_pA  - lb_qA;

    RegZ   = double(s.lambda(1)*A + s.lambda(2)*WW);
    dat    = PGdistribute('Collect');
    save(fullfile(s.result_dir,['train' s.result_name '.mat']),...
        'Wa','Wv','WWa','WWv','mat','dat','ss','A','B','RegZ','mu','s','noise','dat');

    if isfield(s,'debug') && s.debug
        subplot(2,2,1); image(ColourPic(mu,s.likelihood)); axis image ij off;
        subplot(2,2,2); imagesc(abs(ss.ZZ)/ss.N);  colorbar; axis image; title ZZ
        subplot(2,2,3); imagesc(abs(WW));  colorbar; axis image; title WW
        subplot(2,2,4); imagesc(abs(WWv)); colorbar; axis image; title WWv
        drawnow
    end
end

