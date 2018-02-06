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
no   = nargout;

switch lower(opt)
case {'init'}

case {'getzz'}
    [varargout{1:no}] = GetZZ(data);

case {'transfz'}
    data = TransfZ(data,varargin{:});
    PrivateData('set',name,nw,data);

case {'addtoz'}
    data = AddToZ(data,varargin{:});
    PrivateData('set',name,nw,data);

case {'updatelatentvariables'}
    [data,varargout{1:no}] = UpdateLatentVariables(data,varargin{:});
    PrivateData('set',name,nw,data);

case {'shapederivatives'}
    s                 = varargin{5};
    s.result_name     = [s.result_name '_' num2str(nw)];
    [varargout{1:no}] = ShapeDerivatives(data,varargin{1:4},s);

case {'appearancederivatives'}
    s                 = varargin{5};
    s.result_name     = [s.result_name '_' num2str(nw)];
    [varargout{1:no}] = AppearanceDerivatives(data,varargin{1:4},s);

case {'meanderivatives'}
    [varargout{1:no}] = MeanDerivatives(data,varargin{:});

case {'computeof'}
    [varargout{1:no}] = ComputeOF(data,varargin{:});

case {'randomz'}
    data = RandomZ(data,varargin{:});
    PrivateData('set',name,nw,data);

case {'addrandz'}
    data = AddRandZ(data,varargin{:});
    PrivateData('set',name,nw,data);

case {'suffstats'}
    [varargout{1:no}] = SuffStats(data,varargin{:});

case {'collect'}
    % This would be done differently within a privacy-preserving setting
    % with this facility disabled.
    varargout{1}    = data;
otherwise
    error('"%s" unknown.');
end

