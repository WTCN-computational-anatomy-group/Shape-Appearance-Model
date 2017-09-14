function varargout = PGworker(opt,name,nw,varargin)
% Receive job to do and pass results back
% FORMAT varargout = PGworker(opt,nw,varargin)
%
% opt - Option (as string)
% nw  - Worker ID
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

data = PrivateData('get',name,nw);

switch lower(opt)
case {'init'}

case {'getzz'}
    [N,Z,ZZ,sS]  = GetZZ(data);
    varargout{1} = N;
    varargout{2} = Z;
    varargout{3} = ZZ;
    varargout{4} = sS;

case {'transfz'}
    Rz     = varargin{1};
    data   = TransfZ(data,Rz);
    PrivateData('set',name,nw,data);

case {'addtoz'}
    zs     = varargin{1};
    data   = AddToZ(data,zs);
    PrivateData('set',name,nw,data);

case {'updatez'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    EA     = varargin{5};
    s      = varargin{6};
    [data,stats] = UpdateAllZ(data,mu,Wa,Wv,noise,EA,s);
    varargout{1} = stats;
    PrivateData('set',name,nw,data);

case {'wvgradhess'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    s      = varargin{5};
    s.result_name = [s.result_name '_' num2str(nw)];
    [gv,Hv,nll]  = WvGradHess(data,mu,Wa,Wv,noise,s);
    varargout{1} = gv;
    varargout{2} = Hv;
    varargout{3} = nll;

case {'wagradhess'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    s      = varargin{5};
    s.result_name = [s.result_name '_' num2str(nw)];
    [ga,Ha,nll]  = WaGradHess(data,mu,Wa,Wv,noise,s);
    varargout{1} = ga;
    varargout{2} = Ha;
    varargout{3} = nll;

case {'computeof'}
    mu     = varargin{1};
    Wa     = varargin{2};
    Wv     = varargin{3};
    noise  = varargin{4};
    s      = varargin{5};
    nll    = ComputeOF(data,mu,Wa,Wv,noise,s);
    varargout{1} = nll;

case {'randomz'}
    K      = varargin{1};
    data   = RandomZ(data,K);
    PrivateData('set',name,nw,data);

case {'addrandz'}
    Sig    = varargin{1};
    data   = AddRandZ(data,Sig);
    PrivateData('set',name,nw,data);

case {'suffstats'}
    [s0,s1,s2]  = SuffStats(data,varargin{:});
    varargout{1} = s0;
    varargout{2} = s1;
    varargout{3} = s2;

case {'collect'}
    % This would be done differently within a privacy-preserving setting
    % with this facility disabled.
    varargout{1} = data;
case {'save'}

otherwise
    error('"%s" unknown.');
end


%==========================================================================

%==========================================================================
