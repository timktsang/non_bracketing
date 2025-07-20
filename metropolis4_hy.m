% metropolis-hastings method
function result = metropolis4_hy(data,old_para,new_para)
result = old_para;
L_new = loglikelihood3_hy(data,new_para);
L_old = loglikelihood3_hy(data,old_para); 
    lu = log(unifrnd(0,1,1,1));
    if lu < L_new + logprior3_hy(new_para) - L_old - logprior3_hy(old_para);
        result = new_para;
    end;
end
