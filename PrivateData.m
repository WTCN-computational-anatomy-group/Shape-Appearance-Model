function varargout = PrivateData(code,name,nw,varargin)
% Get/get private data
% FORMAT varargout = PrivateData(code,nw,varargin)
%
% code - A string. Either 'init', 'set' or 'get'
% name - a string identifying this job
% nw   - Worker number
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

privdir = tempdir;
privnam = 'private_';
fullname = fullfile(privdir,[privnam name '_' num2str(nw) '.mat']);
switch lower(code)
case {'init','set'}
    data = varargin{1};
    save(fullname,'data')
case {'get'}
    load(fullname,'data');
    varargout{1} = data;
otherwise
    error('Unknown option.');
end


