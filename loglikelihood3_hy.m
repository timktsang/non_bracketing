function result = loglikelihood3_hy(para3,para3_hy)
np = 8;
a = log(para3(1:np).^(para3_hy(1)-1)*(gamma(np*para3_hy(1)))^(1/np)/gamma(para3_hy(1)));
b = log(para3((np+1):(2*np)).^(para3_hy(2)-1)*(gamma(np*para3_hy(2)))^(1/np)/gamma(para3_hy(2)));

result= sum(a(a<Inf & a>-Inf))+...
 sum(b(b<Inf & b>-Inf));

end




