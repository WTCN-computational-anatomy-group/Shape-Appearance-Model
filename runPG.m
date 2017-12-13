function runPG(jsn_filenames, jsn_settings)
% Run The shape-appearance model on "imported" scans
% For help, type:
% runPG --help
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if nargin>=2
    settings = spm_jsonread(jsn_settings);
    s = struct('Ka',          64,...
               'Kv',          64,...
               'linked',      true,...
               'vx',          [1.5 1.5 1.5],...
               'v_settings',  [1e-4  1e-1 2 0.25 0.5]*0.1,...
               'a_settings',  [1e-3  1e-1 0]*10,...
               'mu_settings', [1e-3  1e-1 0],...
               'mg_its',      [3 3],...
               'int_args',    8,...
               'nit',         1,...
               'maxit',       8,...
               'omega',       1.0,...
               'wt',          [1.0 1.0],...
               'nu_factor',   1,...
               'likelihood',  'multinomial',...
               'slices',      [],...
               'result_name', 'test',...
               'result_dir',  tempdir,...
               'ondisk',      true,...
               'batchsize',   4,...
               'verbose',     false);

    if isfield(settings,'K')
        s.Ka     = settings.K;
        s.Kv     = settings.K;
        s.linked = true;
    end

    fldnames = fieldnames(settings);
    for i=1:numel(fldnames)
        if isfield(s,fldnames{i})
            s.(fldnames{i}) = settings.(fldnames{i});
            s.(fldnames{i}) = s.(fldnames{i})(:)';
        end
    end
end

if nargin<1 || (ischar(jsn_filenames) && (strcmp(jsn_filenames,'--help') || strcmp(jsn_filenames,'--h')))
    show_instructions
else
    files    = spm_jsonread(jsn_filenames);
    dat      = struct('f',files);
    for n=1:numel(files)
        dat(n).f = char(dat(n).f{:});
    end

    %disp('Running PG with the following settings...');
    %disp(s)

    PG(dat,s)
end




function show_instructions
disp('Usage: runPG filenames.jsn settings.jsn');
disp(' ');
disp('Images must first be imported using SPM12''s Segment tool, giving a series of rc1*.nii and rc2*.nii files.');
disp('The filenames.jsn is a JSON file with the following general structure, where SCANXX is some filename:');
disp('[["rc1SCAN01.nii","rc2SCAN01.nii"],');
disp(' ["rc1SCAN02.nii","rc2SCAN02.nii"],');
disp(' ["rc1SCAN03.nii","rc2SCAN03.nii"],');
disp('        :                :         ');
disp(' ["rc1SCANXX.nii","rc2SCANXX.nii"]]');
disp(' ');
disp('The settings.jsn file contains various settings, which are of the following form:');
disp('{"K":64,                                    % Number of shape and appearance components');
disp(' "vx":[1.5 1.5 1.5],                        % Voxel sizes (mm - matches those of "imported" images)');
disp(' "v_settings":[1e-05,0.01,0.2,0.025,0.05],  % Registration regularisation settings.');
disp(' "a_settings":[0.01,1,0],                   % Appearance regularisation settings.');
disp(' "mu_settings":[0.001,0.1,0],               % Mean image regularisation settings.');
disp(' "maxit":8,                                 % Number of iterations.');
disp(' "result_name":"test",                      % Name of results files.');
disp([' "result_dir":"' tempdir '",                      % Name of result directory.']);
disp(' "batchsize":4}                             % Batch-size (for local parallelisation).');
disp(' ');
disp('If no settings.jsn file is provided, or if some settings are not specified,');
disp('then default values (shown in the example above) are assumed.');
disp(' ');
disp('Depending on the amount of data, execution can take a very long time.');
disp('Outputs are saved at each iteration in the specified result_dir.');
disp(' ');
disp('The latent variables of interest may be obtained via MATLAB:');
disp('    load(fullfile(settings.result_dir,[''train'' settings.result_name ''.mat'']))');
disp('    Z = cat(2,dat.z);');
disp(' ');


