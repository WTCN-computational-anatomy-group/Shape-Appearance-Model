
d        = size(Wa);
d        = d(1:4);
Kv       = size(Wv,5);
Ka       = size(Wa,5);
K        = Ka+Kv;

scal     = 1;
kd       = 8;
scales   = (-4:2:4);

S = inv(sqrtm(EA));

if s.linked
    n    = 1;
    pic0 = [];
    for k=1:Kv
        pic  = [];
        for sc = scales*S(k,k),
            a0   = single(mu+sc*Wa(:,:,:,:,k));
            v0   = sc*Wv(:,:,:,:,k);
            iphi = GetIPhi(v0,s);
            a1   = Resamp(a0,iphi,s.bs_args);
            pic  = [pic ColourPic(a1,s.likelihood)];
        end
        pic0 = [pic0; pic];
        imagesc(pic0); axis image ij off; drawnow;
        if ~rem(k,kd)
            fname = fullfile(s.result_dir,[s.result_name '_Modes' num2str(n) '.png']);
            imwrite(pic0,fname);
            n    = n + 1;
            pic0 = [];
        end
    end
    if ~isempty(pic0)
        fname = fullfile(s.result_dir,[s.result_name '_Modes' num2str(n) '.png']);
        imwrite(pic0,fname);
    end

else
    n    = 1;
    pic0 = [];
    for k=1:Kv
        pic  = [];
        for sc = scales*S(k,k), 
            v0   = sc*Wv(:,:,:,:,k);
            iphi = GetIPhi(v0,s);
            mu1  = zeros(d,'single');
            for l=1:d(4)
                muc          = spm_diffeo('bsplinc',single(mu(:,:,:,l)),s.bs_args);
                mu1(:,:,:,l) = spm_diffeo('bsplins',muc,iphi,s.bs_args);
            end
            pic  = [pic ColourPic(mu1,s.likelihood)];
        end
        pic0 = [pic0; pic];
        imagesc(pic0); axis image ij off; drawnow;
        if ~rem(k,kd)
            fname = fullfile(s.result_dir,[s.result_name '_ShapeModes' num2str(n) '.png']);
            imwrite(pic0,fname);
            n    = n + 1;
            pic0 = [];
        end
    end
    if ~isempty(pic0)
        fname = fullfile(s.result_dir,[s.result_name '_ShapeModes' num2str(n) '.png']);
        imwrite(pic0,fname);
        imwrite(pic0,['PGapsep_shape_modes' num2str(n) '.png']);
    end

    n    = 1;
    pic0 = [];
    for k=1:Ka
        pic  = [];
        for sc = scales*S(k,k), %-0.2:0.1:0.2,
            mu1 = mu+sc*Wa(:,:,:,:,k);
            pic = [pic ColourPic(mu1,s.likelihood)];
        end
        pic0 = [pic0; pic];
        imagesc(pic0); axis image ij off; drawnow;
        if ~rem(k,kd)
            fname = fullfile(s.result_dir,[s.result_name '_AppModes' num2str(n) '.png']);
            imwrite(pic0,fname);
            n    = n + 1;
            pic0 = [];
        end
    end
    if ~isempty(pic0)
        fname = fullfile(s.result_dir,[s.result_name '_AppModes' num2str(n) '.png']);
        imwrite(pic0,fname);
    end
end


