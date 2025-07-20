function result = loglikelihood4_hy(para4,para4_hy)
np = 11;
a = log(para4(1:np).^(para4_hy(1)-1)*(gamma(np*para4_hy(1)))^(1/np)/gamma(para4_hy(1)));
b = log(para4((np+1):(2*np)).^(para4_hy(2)-1)*(gamma(np*para4_hy(2)))^(1/np)/gamma(para4_hy(2)));

result= sum(a(a<Inf & a>-Inf))+...
 sum(b(b<Inf & b>-Inf));
end




