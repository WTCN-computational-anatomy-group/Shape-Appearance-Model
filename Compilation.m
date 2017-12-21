% Script for compiling the MATLAB code for the Shape & Appearance models.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

compile_dir = fullfile(tempdir,'compiled');

if ~exist('spm','file')
    error('SPM not on search path.');
end

if ~strcmp(spm('ver'),'SPM12')
    error('Wrong SPM version installed (should be SPM12).');
end

spm_file  = which('spm');
[pth,~,~] = fileparts(spm_file);

addpath(pth);
addpath(fullfile(pth,'toolbox','Shoot'));

[ok,message] = mkdir(compile_dir);
if ok
    fprintf('Compiled results should be in "%s"\n', compile_dir)
    mcc('-m','runPG',   '-v','-d',compile_dir);
    mcc('-m','runPGfit','-v','-d',compile_dir); 
    disp('Done.');
else
    error(message);
end
