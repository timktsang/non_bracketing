% log prior
function result = logprior3_hy(para)
%%%% model 1
lik = log(unifpdf(para(1,1),0,1000)) + log(unifpdf(para(1,2),0,1000));
result = lik;
end
