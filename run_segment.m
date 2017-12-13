function run_segment

images=spm_select(Inf,'nifti');

parfor i=1:size(images,1)
    image = deblank(images(i,:));
    fprintf('Segmenting %s... ', image);
    tic
    do_seg(image)
    fprintf('...%g s\n', toc);
end

function do_seg(image)
clear obj
obj.bb       =  NaN(2,3);
obj.bb       = [-90 -126 -72; 90 90 108];
obj.vox      =  2;
obj.cleanup  =  1;
obj.mrf      =  2;
obj.affreg   = 'mni';
obj.reg      = [0 0.001 0.5 0.05 0.2]*0.1;
obj.fwhm     = 1;
obj.samp     = 3;
obj.lkp      = [1 1 2 2 3 3 4 4 5 5 5 6 6];
obj.biasreg  = 0.001*(1/5);
obj.biasfwhm = 60;

tpmname      = fullfile(spm('dir'),'tpm','TPM.nii');
obj.lkp      = [1 1 2 2 3 3 4 4 5 5 5 6 6];

%[pth,~,~]    = fileparts(mfilename('fullpath'));
%tpmname      = fullfile(pth,'BlaiottaTPM.nii');
%obj.lkp      = [1 1 2 2 3 3 4 4 4 5 6 6 7 7 7];


obj.tpm      = spm_load_priors8(tpmname);
obj.image    = spm_vol(char(image));


M = obj.image(1).mat;
c = (obj.image(1).dim+1)/2;
obj.image(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]    = spm_maff8(obj.image(1),8,(0+1)*16,obj.tpm,[],obj.affreg); % Closer to rigid
Affine1          = Affine1*(obj.image(1).mat/M);
obj.image(1).mat = M;

% Run using the origin from the header
[Affine2,ll2]    = spm_maff8(obj.image(1),8,(0+1)*16,obj.tpm,[],obj.affreg); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, obj.Affine  = Affine1; else obj.Affine  = Affine2; end

% Initial affine registration.
obj.Affine     = spm_maff8(obj.image(1),obj.samp*2,(obj.fwhm+1)*16,obj.tpm, obj.Affine, obj.affreg); % Closer to rigid
obj.Affine     = spm_maff8(obj.image(1),obj.samp*2, obj.fwhm,      obj.tpm, obj.Affine, obj.affreg);

res = spm_preproc8(obj);

% Final iteration, so write out the required data.
required = false(max(obj.lkp),4);
required(1:3,2) = true;
spm_preproc_write8(res,required,false(1,2),false(1,2),obj.mrf,obj.cleanup,obj.bb,obj.vox,'.');



