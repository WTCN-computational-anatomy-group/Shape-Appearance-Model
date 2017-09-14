d  = [size(mu) 1 1];
[id{1:3}] = ndgrid(1:d(1),1:d(2),1:d(3));
Ka  = size(Wa,5);
Kv  = size(Wv,5);
K=Ka;

N     = numel(dat);
%N     = min(N,100);
fname = fullfile(s.result_dir,['modes' '_movie.avi']);
avi   = VideoWriter(fname);
open(avi);


ii = 1;
for jj=1:6
    for ph = (0:19),
        sc1  = sin(ph/20*2*pi)*4;
        sc2  = sin(ph/20*pi)*0.5+0.5;
        z     = zeros(K,1);
        z(jj) = sc1*sc2/sqrt(EA(jj,jj));
        
        if ph<=10,
           jp = rem(jj-2+6,6)+1;
        else
           jp = rem(jj,6)+1;
        end
        z(jp) = sc1*(1-sc2)/sqrt(EA(jp,jp));
        a0   = GetA0(z,Wa,mu);
        iphi = GetIPhi(z,Wv,s);
        a1   = Resamp(a0,iphi,s.bs_args);
        %mu1  = Resamp(mu,iphi,s.bs_args);
        pic  = [ColourPic(a1,s.likelihood)];
        pic(:,:,3)=0;
        image(pic); axis ij image off; drawnow
        axis image off;
        drawnow
        currFrame = struct('cdata',pic,'colormap',[]);
        writeVideo(avi,currFrame);
    end
    fprintf('.');
end

close(avi);

