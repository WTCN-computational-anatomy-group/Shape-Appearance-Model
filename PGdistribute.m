function varargout = PGdistribute(opt,varargin)
% Distribute functions among workers
% FORMAT varargout = PGdistribute(opt,varargin)
%
% opt - A string, which may be one of:
%       'Init'       - Initialise the system, passing various arguments
%                      via the s data structure.
%       'Share'      - Pass data among workers.  This option would be
%                      removed for the real distributed system as the data
%                      would be held locally on each worker.
%       'Collect'    - Collect the data from each worker.  This is another
%                      option to be removed for the properly distributed
%                      version, as this data should not leave the sites.
%       'GetZZ'      - Return sufficient statistics from latent variables.
%       'TransfZ'    - Matrix pre-multiplication of latent variables.
%       'AddToZ'     - Addition of vector to latent variables.
%       'UpdateZall' - Re-estimate all latent variables and return useful
%                      suffient statistics.
%       'WvGradHess' - Return derivatives w.r.t. shape basis functions.
%       'WaGradHess' - Return derivatives w.r.t. appearance basis functions.
%       'MuGradHess' - Return derivatives w.r.t. appearance basis functions.
%       'ComputeOF'  - Return the likelihood part of objective function.
%       'RandomZ'    - Initialise latent variables randomly.
%       'AddRandZ'   - Add random nouse to latent variables.
%       'SuffStats'  - Compute sufficient statistics.
%
% varargin - See each function for descriptions of arguments.
%
% This code is intended to eventually distribute work among many different
% hospital sites. The idea is that much of the computation can be done
% using aggregates, such that the actual data from individual patients
% should not leave the sites. 
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


persistent name par NW
if isempty(name), name = 'unknown'; end
if isempty(par)
    par  = false;
    NW   = 1;
else
    if isempty(NW)
        p  = gcp;
        NW = p.NumWorkers;
    end
end

%fprintf('\n*** %s : ', opt);
%tic

drawnow

switch lower(opt)
case {'init'}
    % Initialise the system, passing various arguments
    % via the s data structure.

    s    = varargin{1};
    name = s.result_name;
    if isfield(s,'par'), par = s.par; end
    if par, p = gcp; NW = p.NumWorkers; else NW = 1; end

case {'share'}
    % Pass data among workers.  This option would be
    % removed for the real distributed system as the data
    % would be held locally on each worker.

    dat = varargin{1};
    N   = numel(dat);
    if par
        parfor nw = 1:NW
           %clear data
            ind  = floor((nw-1)*(N/NW)+1):floor(nw*N/NW);
            data = dat(ind);
            PrivateData('set',name,nw,data); % Private should mean private
        end
    else
        PrivateData('set',name,1,dat); % Private should mean private
    end

case {'collect'}
    % Collect the data from each worker.  This is another
    % option to be removed for the properly distributed
    % version, as this data should not leave the sites.

    dat_cell = cell(1,NW);
    if par
        parfor nw=1:NW
            tmp = PGworker('Collect',name, nw);
            if ~isempty(tmp)
                dat_cell{nw} = tmp(:);
            end
        end
    else
        dat_cell{1} = PGworker('Collect',name,1);
    end
    varargout{1} = cat(1,dat_cell{:});

case {'getzz'}
    N  = 0;
    Z  = 0;
    ZZ = 0;
    sS = 0;
    if par
        parfor nw=1:NW
            [n,z,zz,ss] = PGworker('getZZ', name, nw);
            if n
                N  = N  + n;
                Z  = Z + z;
                ZZ = ZZ + zz;
                sS = sS + ss;
            end
        end
    else
        [N,Z,ZZ,sS] = PGworker('getZZ', name, 1);
    end
    varargout{1} = N;
    varargout{2} = Z;
    varargout{3} = ZZ;
    varargout{4} = sS;

case {'transfz'}
    T = varargin{1};
    if par
        parfor nw=1:NW
            PGworker('transfZ',name, nw,T);
        end
    else
        PGworker('transfZ', name, 1,T);
    end

case {'addtoz'}
    zs = varargin{1};
    if par
        parfor nw=1:NW
            PGworker('AddToZ',name, nw,zs);
        end
    else
        PGworker('AddToZ', name, 1,zs);
    end

case {'updatezall'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    A      = varargin{5};
    s      = varargin{6};
    if par
        CompSmo = false;
        CompMu  = false;

        N  = 0;
        L  = 0;
        s0 = 0;
        s1 = 0;
        Z  = 0;
        ZZ = 0;
        sS = 0;
        if CompSmo
            SmoSuf = 0;
        end
        if CompMu
            gmu    = single(0);
            Hmu    = single(0);
        end

        parfor nw=1:NW
            st  = PGworker('UpdateZall',name, nw,mu,Wa,Wv,noise,A,s);
            if ~isempty(st)
                N      = N  + st.N;
                L      = L  + st.L;
                s1     = s1 + st.s1;
                s0     = s0 + st.s0;
                Z      = Z  + st.Z;
                ZZ     = ZZ + st.ZZ;
                sS     = sS + st.sS;
                if CompSmo
                    SmoSuf = SmoSuf + st.SmoSuf;
                end
                if CompMu
                    gmu    = gmu    + st.gmu;
                    Hmu    = Hmu    + st.Hmu;
                end
            end
        end
        stats = struct('N',N, 'Z',Z, 'ZZ',ZZ, 'sS',sS, 's0',s0, 's1',s1, 'L', L);
        if CompSmo
            stats.SmoSuf = SmoSuf;
        end
        if CompMu
            stats.gmu    = gmu;
            stats.Hmu    = Hmu;
        end
    else
        stats  = PGworker('UpdateZall',name, 1,mu,Wa,Wv,noise,A,s);
    end
    varargout{1} = stats;

case {'wvgradhess'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    s      = varargin{5};

    gv_cell = cell(1,NW);
    Hv_cell = cell(1,NW);
    nll     = 0;
    if par
        parfor nw=1:NW
            [gv_cell{nw},Hv_cell{nw},nll1] = PGworker('WvGradHess',name,nw, mu,Wa,Wv,noise,s);
            nll = nll + nll1;
        end
    else
        [gv_cell{1},Hv_cell{1},nll] = PGworker('WvGradHess',name,1, mu,Wa,Wv,noise,s);
    end

    % Needs to deal with possibly empty results
    rw = 0;
    for nw=1:NW
        if ~isempty(gv_cell{nw})
            rw = nw;
            break;
        end
    end
    if rw
        for nw=(rw+1):NW
            if ~isempty(gv_cell{nw})
                for k=1:size(gv_cell{rw},5)
                    gv_cell{rw}(:,:,:,:,k) = gv_cell{rw}(:,:,:,:,k) + gv_cell{nw}(:,:,:,:,k);
                    Hv_cell{rw}(:,:,:,:,k) = Hv_cell{rw}(:,:,:,:,k) + Hv_cell{nw}(:,:,:,:,k);
                end
                DeleteFA(gv_cell{nw});
                DeleteFA(Hv_cell{nw});
            end
        end
        varargout{1} = gv_cell{rw};
        varargout{2} = Hv_cell{rw};
        varargout{3} = nll;
    else
        varargout{1} = [];
        varargout{2} = [];
        varargout{3} = 0;
    end

case {'wagradhess'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    s      = varargin{5};

    ga_cell = cell(1,NW);
    Ha_cell = cell(1,NW);
    nll     = 0;
    if par
        parfor nw=1:NW
            [ga_cell{nw},Ha_cell{nw},nll1] = PGworker('WaGradHess',name,nw,mu,Wa,Wv,noise,s);
            nll = nll + nll1;
        end
    else
        [ga_cell{1},Ha_cell{1},nll] = PGworker('WaGradHess',name,1,mu,Wa,Wv,noise,s);
    end

    % Needs to deal with possibly empty results
    rw = 0;
    for nw=1:NW
        if ~isempty(ga_cell{nw})
            rw = nw;
            break;
        end
    end
    if rw
        for nw=(rw+1):NW
            if ~isempty(ga_cell{nw})
                for k=1:size(ga_cell{rw},5)
                    ga_cell{rw}(:,:,:,:,k) = ga_cell{rw}(:,:,:,:,k) + ga_cell{nw}(:,:,:,:,k);
                    Ha_cell{rw}(:,:,:,:,k) = Ha_cell{rw}(:,:,:,:,k) + Ha_cell{nw}(:,:,:,:,k);
                end
                DeleteFA(ga_cell{nw});
                DeleteFA(Ha_cell{nw});
            end
        end
        varargout{1} = ga_cell{rw};
        varargout{2} = Ha_cell{rw};
        varargout{3} = nll;
    else
        varargout{1} = [];
        varargout{2} = [];
        varargout{3} = 0;
    end

case {'mugradhess'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    s      = varargin{5};

    nll      = 0;
    gmu_cell = cell(1,NW);
    Hmu_cell = cell(1,NW);
    if par
        parfor nw=1:NW
            [gmu_cell{nw},Hmu_cell{nw},nll1] = PGworker('muGradHess',name,nw,mu,Wa,Wv,noise,s);
            nll = nll + nll1;
        end
    else
        [gmu_cell{1},Hmu_cell{1},nll] = PGworker('muGradHess',name,1,mu,Wa,Wv,noise,s);
    end

    % Needs to deal with possibly empty results
    rw = 0;
    for nw=1:NW
        if ~isempty(gmu_cell{nw})
            rw = nw;
            break;
        end
    end
    if rw
        for nw=(rw+1):NW
            if ~isempty(gmu_cell{nw})
                gmu_cell{rw}(:,:,:,:) = gmu_cell{rw}(:,:,:,:) + gmu_cell{nw}(:,:,:,:);
                Hmu_cell{rw}(:,:,:,:) = Hmu_cell{rw}(:,:,:,:) + Hmu_cell{nw}(:,:,:,:);
            end
        end
        varargout{1} = gmu_cell{rw};
        varargout{2} = Hmu_cell{rw};
        varargout{3} = nll;
    else
        varargout{1} = [];
        varargout{2} = [];
        varargout{3} = [];
    end

case {'computeof'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    s      = varargin{5};
    nll    = 0;
    if par
        parfor nw=1:NW
            nll1   = PGworker('ComputeOF',name,nw,mu,Wa,Wv,noise,s);
            nll    = nll + nll1;
        end
    else
        nll   = PGworker('ComputeOF',name,1,mu,Wa,Wv,noise,s);
    end
    varargout{1} = nll;

case {'randomz'}
    K = varargin{1};
    if par
        parfor nw=1:NW
            PGworker('randomz',name,nw,K);
        end
    else
        PGworker('randomz',name,1,K);
    end

case {'addrandz'}
    Sig = varargin{1};
    if par
        parfor nw=1:NW
            PGworker('addrandz',name,nw,Sig);
        end
    else
        PGworker('addrandz',name,1,Sig);
    end

case {'suffstats'}
    s0 = single(0);
    s1 = single(0);
    s2 = single(0);
    if par
    parfor nw=1:NW
        [t0,t1,t2,mat]   = PGworker('Suffstats',name,nw,varargin{:});
        if ~isempty(t0)
            s0 = s0 + t0;
            s1 = s1 + t1;
            s2 = s2 + t2;
        end
    end
    else
        [s0,s1,s2,mat]   = PGworker('Suffstats',name,1,varargin{:});
    end
    varargout{1} = s0;
    varargout{2} = s1;
    varargout{3} = s2;
    varargout{4} = mat;

otherwise
    error('"%s" unknown.', opt);
end

drawnow

%fprintf('%g sec\n', toc);

