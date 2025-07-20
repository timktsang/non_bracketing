function r = drchrnd(a,b,n)
p = length(a);
r = gamrnd(repmat(a+b,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

end