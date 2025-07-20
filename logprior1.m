% log prior
function result = logprior1(para)
%%%% model 1
lik = log(unifpdf(para(1,1),0,30)) + log(unifpdf(para(1,2),0,30)) + log(unifpdf(para(1,3),0,30)) + ...
      log(unifpdf(para(1,4),0,30)) + log(unifpdf(para(1,5),0,30)) + log(unifpdf(para(1,6),0,30)) + ...
      log(unifpdf(para(1,7),-0.001,0.001));
result = lik;
end
