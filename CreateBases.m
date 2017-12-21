function [mu_fa,Wa,Wv,WWa,WWv,WW] = CreateBases(s,mu,mat)
% Initialise basis functions
% FORMAT [mu_fa,Wa,Wv,WWa,WWv,WW] = CreateBases(s,mu,mat)
%
% s     - Various settings. Uses s.ondisk, s.Ka, s.Kv, s.linked, s.result_dir
%         s.result_name, s.a_settings, s.v_settings and s.mu_settings.
% mu    - The mean basis function
% mat   - Voxel-to-world mapping to use.
%
% mu_fa - A file array object for the mean (if required)
% Wa    - Appearance basis functions (possibly a file_array)
% Wv    - Shape basis functions (possibly a file_array)
% WWa   - Empty matrix for Wa'*La*Wa
% WWv   - Empty matrix for Wv'*Ta*Wv
% WW    - Empty matrix for combined WWa and WWv
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

Ka     = s.Ka;
Kv     = s.Kv;
if s.linked
    Koff = 0;
    Kv   = Ka;
else
    Koff = Ka;
end
K     = max(Ka,Kv+Koff);
inda  = 1:Ka;
indv  = Koff+(1:Kv);
if isempty(inda), inda = []; end
if isempty(indv), indv = []; end

d   = [size(mu) 1 1]; d  = d(1:4);
if nargin<3, mat = eye(4); end

if isfield(s,'ondisk') && s.ondisk
    % Do this if basis functions are to be written to disk as file_arrays.
    if numel(indv)>0
        % Generate empty file for shape components
        Nv0         = nifti;
        Nv0.dat     = file_array(fullfile(s.result_dir,[s.result_name '_Wv.nii']),[d(1:3) 3 numel(indv)],'float32',352);
        Nv0.mat     = mat;
        Nv0.mat0    = mat;
        Nv0.descrip = sprintf('PG velocities (%g %g %g %g %g)', s.v_settings);
        create(Nv0);
        Wv_fa       = Nv0.dat;
        Wv          = Nv0.dat;
        for k=1:numel(indv), Wv_fa(:,:,:,:,k) = 0; end
    else
        Wv          = zeros([d(1:3) 3 0],'single');
    end
    if numel(inda)>0
        % Generate empty file for appearance components
        Na0         = nifti;
        Na0.dat     = file_array(fullfile(s.result_dir,[s.result_name '_Wa.nii']),[d(1:4) numel(inda)],'float32',352);
        Na0.mat     = mat;
        Na0.mat0    = mat;
        Na0.descrip = sprintf('PG appearance (%g %g %g)', s.a_settings);
        create(Na0);
        Wa_fa       = Na0.dat;
        Wa          = Na0.dat;
        for k=1:numel(inda), Wa_fa(:,:,:,:,k) = 0; end
    else
        Wa          = zeros([d 0],'single');
    end

    % Generate file for mean
    Nmu         = nifti;
    Nmu.dat     = file_array(fullfile(s.result_dir,[s.result_name '_mu.nii']),d(1:4),'float32',352);
    Nmu.mat     = mat;
    Nmu.mat0    = mat;
    Nmu.descrip = sprintf('PG mean (%g %g %g)', s.mu_settings);
    create(Nmu);
    Nmu.dat(:,:,:,:) = mu;
    mu_fa       = Nmu.dat;
else
    % Do this if basis functions are only to be held in memory.
    mu_fa = [];
    Wa = zeros([d        numel(inda)],'single');
    Wv = zeros([d(1:3) 3 numel(indv)],'single');
end

% Create empty matrices
WWa           = zeros(numel(inda));
WWv           = zeros(numel(indv));
WW            = zeros(K);
WW(inda,inda) = WW(inda,inda) + WWa;
WW(indv,indv) = WW(indv,indv) + WWv;

