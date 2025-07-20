% likelihood function 

function [result result2] = loglikelihood1(data1,para1,cum_ILILAB,cutoff)

pro = ILI_pro(data1(:,[1 2 4]),para1,cum_ILILAB,cutoff); 

h = pro.*exp(para1(7)*log2(data1(:,3)/5));
p = binopdf(data1(:,5)./data1(:,3) >= 4,1,1-exp(-h));
result = sum(log(p));
result2 = 1-exp(-h);
end



