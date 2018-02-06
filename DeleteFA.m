function DeleteFA(varargin)
% Delete file_arrays from disk
% FORMAT DeleteFA(f1 [,f2 [,f3]])
% f1, f2 etc - Possible file_arrays.
%
% For 3D data with lots of components, there are too many parameters
% in the component fields to save in memory. These are instead saved using
% the file_array data structure of SPM, which allows relatively quick I/O.
% This function deletes these file_arrays from disk if they exist.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

for i=1:nargin
    f = varargin{i};
    if isa(f,'file_array')
        delete(f.fname);
    elseif isa(f,'char')
        delete(f);
    end
end

