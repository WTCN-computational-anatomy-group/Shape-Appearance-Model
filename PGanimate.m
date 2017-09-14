

d  = [size(mu) 1 1];
[id{1:3}] = ndgrid(1:d(1),1:d(2),1:d(3));
Ka  = size(Wa,5);
Kv  = size(Wv,5);
if Ka==Kv && size(dat(1).z,1)==Ka, Koff = 0; else, Koff = Ka; end

N     = numel(dat);
%N     = min(N,100);
fname = fullfile(s.result_dir,[s.result_name '_movie.avi']);
avi   = VideoWriter(fname);
open(avi);

ii = 1;
z1 = dat(ii).z;
f1 = GetDat(dat(ii),s);f1(~isfinite(f1))=0;
sc = cos(0:pi/8:pi)*0.5+0.5;

for jj=1:N
    z0 = z1;
    f0 = f1;
    ii = rem(ii,N)+1;
    z1 = dat(ii).z;
    f1 = GetDat(dat(ii),s);f1(~isfinite(f1))=0;

    for i=1:numel(sc)
        z    = (z0*sc(i)+z1*(1-sc(i)));
        a0   = GetA0(z,Wa,mu);
        iphi = GetIPhi(z,Wv,s);
        ff   = f0*sc(i)+f1*(1-sc(i));
        a1   = Resamp(a0,iphi,s.bs_args);
        mu1  = Resamp(mu,iphi,s.bs_args);
        pic  = [ColourPic(ff)               ColourPic(a1,s.likelihood)
                ColourPic(mu1,s.likelihood) ColourPic(a0,s.likelihood)];
        %pic = ColourPic(a1,s.likelihood);
        %image([ColourPic(ff, s.likelihood) pic]);
        image(pic); axis ij image off; drawnow
        axis image off;
        drawnow
        currFrame = struct('cdata',pic,'colormap',[]);
        writeVideo(avi,currFrame);
    end
    fprintf('.');
end

close(avi);

