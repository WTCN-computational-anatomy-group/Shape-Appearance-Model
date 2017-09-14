function ld = FieldLogDet(d,prm)
% Still work in progress and in need of thorough checking.
% FORMAT ld = FieldLogDet(d,prm)
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$
d = [d(:)' 1 1];
d = d(1:3);

if numel(prm)==6
    % Compute for regular fields
    spm_field('boundary',0);
    F = zeros(d(1:3),'single');
    F(1,1,1) = 1;
    F = spm_field('vel2mom',F,prm);
    F = real(fftn(F));
    if prm(4)==0
        F(1,1,1) = 0;
    end
    msk = F(:)>0;
    ld  = [-(sum(log(F(msk))) - 2*sum(msk)*sum(log(prm(1:3)))), sum(msk)];

else
    % Compute for velocity field regularisation
    spm_diffeo('boundary',0);
    F   = spm_diffeo('kernel',d,prm);
    if size(F,4) == 1,
        % The differential operator is symmetric, so the Fourier transform should be real
        F = 1./real(fftn(F));
        if prm(4)==0, F(1,1) = 0; end;
        msk = F(:)>0;
        ld  = [-(3*sum(log(F(msk))) - 2*sum(msk)*sum(log(prm(1:3)))), sum(msk)];
    else
        for j=1:3,
            for i=1:3,
                % The differential operator is symmetric, so the Fourier transform should be real
                A(:,:,:,i,j) = real(fftn(F(:,:,:,i,j)));
            end
        end
        % Compare the following with computing the determinant of a 3x3 matrix...
        dt  = A(:,:,:,1,1).*(A(:,:,:,2,2).*A(:,:,:,3,3) - A(:,:,:,2,3).*A(:,:,:,3,2)) +...
              A(:,:,:,1,2).*(A(:,:,:,2,3).*A(:,:,:,3,1) - A(:,:,:,2,1).*A(:,:,:,3,3)) +...
              A(:,:,:,1,3).*(A(:,:,:,2,1).*A(:,:,:,3,2) - A(:,:,:,2,2).*A(:,:,:,3,1));

        msk     = find(dt<=0);
        dt      = 1./dt;
        dt(msk) = 0;
        if prm(4)==0, dt(1,1,1) = 0; end;
        msk = dt(:)>0;
        ld  = [-sum(log(dt(msk))), sum(msk)];
    end
end

