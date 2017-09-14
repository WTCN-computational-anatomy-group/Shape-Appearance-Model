%load mnist/mnist_test.mat
load mnist/mnist_hard.mat

result_dir = s.result_dir;

dat = struct('f',F);
lab = y;
d   = size(dat(1).f);
[id{1:3}] = ndgrid(1:d(1),1:d(2),1);
clear train
nufac = 0;
for c=1:10,
    tmp = [result_name num2str(Ntrain) '_' num2str(c-1)];
    result_file = fullfile(result_dir,['train' tmp '.mat']);
    train(c)    = load(result_file);
    train(c).mu = train(c).mu(:,:,:,:,:);
    train(c).Wa = train(c).Wa(:,:,:,:,:);
    train(c).Wv = train(c).Wv(:,:,:,:,:);

    if false
    th  = max(diag(train(c).ss.ZZ))*1e-3;
    ind = find(diag(train(c).ss.ZZ>th));
    tmp = inv(train(c).EA); train(c).EA = inv(tmp(ind,ind));
    tmp = inv(train(c).RegZ); train(c).RegZ = inv(tmp(ind,ind));
    train(c).Wa = train(c).Wa(:,:,:,:,ind);
    train(c).Wv = train(c).Wv(:,:,:,:,ind);
    end
    nufac = nufac + train(c).noise.nu_factor;
end
nufac = nufac/10;
for c=1:10,
    train(c).noise.nu_factor = nufac;
end

corr  = 0;
incor = 0;
nll0  = 0;
nll   = 0;
LL    = nan(numel(dat),numel(train));

for i=1:numel(dat)

    f = GetDat(dat(i));
    parfor c=1:numel(train)
        s     = train(c).s;
        s.nit = 20;
        K     = size(train(c).Wv,5);
        z     = zeros(K,1);
        for reg = 2.^(5:-1:0)
           [z,~,~] = UpdateZ(z,f,train(c).mu,train(c).Wa,train(c).Wv,reg*train(c).RegZ,s,train(c).noise);
        end
        s.nit = 100;
        [z,~,L] = UpdateZ(z,f,train(c).mu,train(c).Wa,train(c).Wv,train(c).RegZ,s,train(c).noise);
        LL(i,c) = L/train(c).noise.nu_factor*nufac;
    end

    ll       = LL(i,:) - LL(i,lab(i)+1);
    [mn,ind] = max(ll);
    corr  = corr  +  (lab(i)==ind-1);
    incor = incor + ~(lab(i)==ind-1);
    fprintf('%3d %d %d  %d  %6.4g  ', i, lab(i),ind-1, lab(i)==ind-1,100*incor/(corr+incor));
    fprintf(' %-6.4g',ll);
    fprintf('\n');
end

save(fullfile(s.result_dir,['MNIST' num2str(Ntrain) '_' result_name '_Xval.mat']));


if false

for ii=1:size(P,1)
load(deblank(P(ii,:)));
[~,guess]=max(LL,[],2);
guess=guess-1;
[~,nam] = fileparts(deblank(P(ii,:)));
fprintf('%s\t%g\n', nam, sum(guess~=lab)/numel(lab))
%disp(s)
end




end
