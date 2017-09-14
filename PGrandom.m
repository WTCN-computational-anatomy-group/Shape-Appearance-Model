%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

d  = [size(mu) 1 1];
d  = d(1:4);
extra = '';
if exist('RegZ','var'), Li = RegZ; else Li = L; end

K  = size(Li,1);
Ka = size(Wa,5);
Kv = size(Wv,5);

if K == Ka+Kv
    linked = false;
    Koff   = Ka;
else
    linked = true;
    Koff   = 0;
end

int_args = [8 1 1];
bs_args  = [2 2 2 1 1 1];
scal = 1;
%scal = 0.5;

clear id
[id{1:3}] = ndgrid(1:d(1),1:d(2),1:d(3));

N1 = 20;
N2 = 20;
N1 = round(min(2000/d(2),N1));
N2 = round(min(2000/d(1),N2));

N1 = min(N1,4);
N2 = min(N2,16);

R   = randn(K,N1*N2);%/sqrt(0.034);

for iter=2:2
if iter==1
    extra = '';
    if exist('RegZ','var'), Li = RegZ; else Li = L; end
else
    Li    = EA;
    extra = 'noreg';
end

clf;
pic0 = [];
for n1=1:N1
    pic = [];
    for n2= 1:N2 %1:min(K,N2),
        v0  = zeros([d(1:3) 3],'single');
        a0  = zeros([d(1:4) 1],'single')+mu;
        z   = real(sqrtm(Li))\R(:,n1+(n2-1)*N1);
        for k=1:Ka,
            a0 = a0 + z(k)*Wa(:,:,:,:,k);
        end
        for k=1:Kv,
            v0 = v0 + z(k+Koff)*Wv(:,:,:,:,k);
        end
        iphi           = Shoot(v0,[s.vx s.v_settings],s.int_args);

        mu1  = zeros(d,'single');
        for l=1:d(4)
            muc          = spm_diffeo('bsplinc',single(a0(:,:,:,l)),s.bs_args);
            mu1(:,:,:,l) = spm_diffeo('bsplins',muc,iphi,s.bs_args);
        end
        pic  = [pic ColourPic(mu1,s.likelihood)];
    end
    pic0 = [pic; pic0];
    imagesc(pic0); axis image ij off; drawnow;
end
%pic0 = convn(pic0,ones(2)/4);
%pic0 = pic0(2:2:end,2:2:end);
imagesc(pic0); axis image ij off; drawnow;

imwrite(pic0,fullfile(s.result_dir,['random1_' s.result_name extra '.png']));

pic0 = [];
for n1=1:N1
    pic = [];
    for n2= 1:N2 %1:min(K,N1),
        v0  = zeros([d(1:3) 3],'single');
        a0  = zeros([d(1:4) 1],'single') + mu;
        z   =-real(Li^(-1/2))*R(:,n1+(n2-1)*N1);
        for k=1:Ka,
            a0 = a0 + z(k)*Wa(:,:,:,:,k);
        end
        for k=1:Kv,
            v0 = v0 + z(k+Koff)*Wv(:,:,:,:,k);
        end
        iphi           = Shoot(v0,[s.vx s.v_settings],s.int_args);

        mu1  = zeros(d,'single');
        for l=1:d(4)
            muc          = spm_diffeo('bsplinc',single(a0(:,:,:,l)),s.bs_args);
            mu1(:,:,:,l) = spm_diffeo('bsplins',muc,iphi,s.bs_args);
        end
        pic  = [pic ColourPic(mu1,s.likelihood)];

    end
    pic0 = [pic; pic0];
    imagesc(pic0); axis image ij off; drawnow;
end
%pic0 = convn(pic0,ones(2)/4);
%pic0 = pic0(2:2:end,2:2:end);
imagesc(pic0); axis image ij off; drawnow;

imwrite(pic0,fullfile(s.result_dir,['random2_' s.result_name extra '.png']));
end

