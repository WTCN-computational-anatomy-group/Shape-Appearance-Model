function [f,mat] = GetDat(dat,s)
% Read in the image data
% FORMAT f = GetDat(dat,s)
%
% dat - Structure containing various information about each image.
%       Fields for each image n are:
%       dat(n).f - Image data.
%       dat(n).z - Expectations of latent variables.
%       dat(n).S - Covariances of latent variables.
% s   - Settings. Uses s.likelihood.
%
% f   - The resulting image
% mat - Voxel to world mapping
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if ~isfield(dat,'f')
    error('No data.');
end
if nargin<2, s.likelihood = 'normal'; end
mat = eye(4);

switch class(dat.f)
case {'nifti','char'}
    if isa(dat.f, 'char')
        [~,~,ext] = fileparts(dat.f(1,:));
        switch lower(ext)
        case {'.jpg','.png','.tif'}
            f = imread(deblank(dat.f(1,:)));
            d = [size(f) 1 1];
            f = single(reshape(f,[d(1:2) 1 d(3)]));
            return
        case {'.nii','.img'}
            Nii = nifti(dat.f);
            mat = Nii.mat;
        otherwise
            error('Unknown file extension (%s).', ext);
        end
    else
        Nii = dat.f;
        mat = Nii.mat;
    end

    D = zeros(numel(Nii),6);
    for l=1:numel(Nii)
        d = [size(Nii(1).dat) 1 1 1 1 1];
        D(l,:) = d(1:6);
    end

    if any(any(diff(D(:,1:3),1,1)))
        error('Incompatible dimensions.');
    end
    d = [D(1,1:3) sum(prod(D(:,4:5),2))];
    slices = 1:d(3);
    if isfield(s,'slices') && ~isempty(s.slices)
        slices = s.slices;
        d(3)   = numel(slices);
    end
    switch s.likelihood
    case {'normal','gaussian'}
    case {'binomial','binary'}
        if d(4)~=1, error('Wrong number of volumes.'); end
    case {'multinomial','categorical'}
        d(4) = d(4) + 1;
    end


    f = zeros(d,'single');
    l = 0;
    for l1 = 1:numel(Nii)
        for l2 = 1:size(Nii(l1).dat,4)
            for l3 = 1:size(Nii(l1).dat,5)
                l = l + 1;
                f(:,:,:,l) = single(Nii(l1).dat(:,:,slices,l2,l3));
            end
        end
    end

    switch s.likelihood
    case {'multinomial','categorical'}
        sf           = sum(f,4);
        f(:,:,:,end) = max(1-sf,0);
    end
case {'double','logical'}
    f = single(dat.f);
case {'single'}
    f = dat.f;
otherwise
    error('Unrecognised data type.');
end

if size(f,4)>1
    msk = any(~isfinite(f),4);
    if any(msk(:))
        f(repmat(msk,[1 1 1 size(f,4)])) = NaN;
    end
end

