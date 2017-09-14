function I = Horder(d)
I = diag(1:d);
l = d;
for i1=2:d
    for i2=1:(i1-1)
        l = l + 1;
        I(i1,i2) = l;
        I(i2,i1) = l;
    end
end

