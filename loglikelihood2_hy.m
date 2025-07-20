
function result = loglikelihood2_hy(para2,para2_hy)
a = log(para2(1:10).^(para2_hy(1)-1)*(gamma(10*para2_hy(1)))^(1/10)/gamma(para2_hy(1)));
b = log(para2(11:20).^(para2_hy(2)-1)*(gamma(10*para2_hy(2)))^(1/10)/gamma(para2_hy(2)));

result= sum(a(a<Inf & a>-Inf))+...
 sum(b(b<Inf & b>-Inf));

end

