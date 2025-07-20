% metropolis-hastings method
function result = metropolis1(data,old_para,new_para,cum_ILILAB,cutoff)
result = old_para;
L_new = loglikelihood1(data,new_para,cum_ILILAB,cutoff);
L_old = loglikelihood1(data,old_para,cum_ILILAB,cutoff); 
    lu = log(unifrnd(0,1,1,1));
    if lu < L_new + logprior1(new_para) - L_old - logprior1(old_para);
        result = new_para;
    end;
end
