
d        = size(Wa);
d        = d(1:4);
Ka       = size(Wa,5);
Kv       = size(Wv,5);

if ~isfield(dat,'z')
    dat      = PGdistribute('Collect');
end

N        = numel(dat);
Z        = cat(2,dat.z);

int_args = 8;
bs_args  = [2 2 0 1 1 1];

N1=floor(sqrt(N));
N2=floor(N/N1);

N1 = min(N1,4);
N2 = min(N2,12);

%N1=5; N2=5;

rand('seed',1);
randn('seed',1);
select = randperm(N);
select = select(1:(N1*N2));
select = sort(select);
%select(1)=1583;

%N1=1; N2=1; N=1; select=411;

pic0e = [];
pic0f = [];
pic0a = [];
pic0s = [];
pic0m = [];
for n1=1:N1
    pice = [];
    picf = [];
    pica = [];
    pics = [];
    picm = [];
    for n2= 1:N2 %1:min(K,20),
        n    = select(n2+((n1-1)*N2));
        f    = GetDat(dat(n),s);
        z    = dat(n).z;
        a0   = GetA0(z,Wa,mu);
        iphi = GetIPhi(z,Wv,s);
        a1   = Resamp(a0,iphi,s.bs_args);
        mu1  = Resamp(mu,iphi,s.bs_args);
        pice = [pice ColourPic(f)];
        picf = [picf ColourPic( a1,s.likelihood)];
        pics = [pics ColourPic(mu1,s.likelihood)];
        pica = [pica ColourPic( a0,s.likelihood)];
    end
    pic0e = [pic0e; pice];
    pic0f = [pic0f; picf];
    pic0s = [pic0s; pics];
    pic0a = [pic0a; pica];
    imagesc(pic0f); axis image ij off; drawnow;
end
%pic0 = convn(pic0,ones(2)/4);
%pic0 = pic0(2:2:end,2:2:end);
imagesc(pic0f); axis image ij off; drawnow;

fname = fullfile(s.result_dir,[s.result_name '_Orig.png']);
imwrite(pic0e,fname);

fname = fullfile(s.result_dir,[s.result_name '_FullFit.png']);
imwrite(pic0f,fname);

fname = fullfile(s.result_dir,[s.result_name '_ShapeFit.png']);
imwrite(pic0s,fname);

fname = fullfile(s.result_dir,[s.result_name '_AppearanceFit.png']);
imwrite(pic0a,fname);

