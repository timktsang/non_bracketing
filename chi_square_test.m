%%% function to compute the DIC

function result = chi_square_test(data,para0,para1,para2,ILILAB)

% first compuate the number of raw data
simN = 500;
data1 = data;
exp_raw = final_size_dist(data1);
exp_raw(:,3) = 0;

% indicator of missing 
miss_1 = data(:,13) == -1 & data(:,14) ~= -1 ;
miss_2 = data(:,13) ~= -1 & data(:,14) == -1 ;
miss_1_2 = data(:,13) == -1 & data(:,14) == -1 ;
PCR = data(:,12) - data(:,11) < 10;

% declare the matrix first
y = min(data(data(:,11)>0,11)):max(data(:,12));
z = y(ones(size(data,1),1),:);
time_matrix = bsxfun(@le,data(:,11),z) + bsxfun(@le,z,data(:,12)) - 1;
hazard_comm = bsxfun(@times,time_matrix,ILILAB(y));

miss_age = data(:,5) == -1;
age_prob = [sum(data(:,5) > 0 & data(:,5) <= 18) sum(data(:,5) > 18 & data(:,5) <= 50) sum(data(:,5) > 50)];


for v = 1:500;

% impute age
data(miss_age,5) = randsample([17 39 55],sum(miss_age),1,age_prob);   
    
% impute the missing titer 
% for missing 1_2 & y_missing_1_2, draw the missing_1 from the normal dist first
data(miss_1_2,13) = 5*2.^((data(miss_1_2,5)<=18).*randsample(0:9,sum(miss_1_2),1,para0(1,1:10))' + (data(miss_1_2,5)>18).*randsample(0:9,sum(miss_1_2),1,para0(1,11:20))');

% for missing 1
at1_prob = (data(miss_1,5) <= 18)*para0(1,1:10) + (data(miss_1,5)>=18)*para0(1,11:20);
fold = max(2^-3,min(2^6,bsxfun(@rdivide,data(miss_1,14),5*2.^(0:9))));
y_temp = fold>=4;
sup_temp = exp(bsxfun(@plus,para1(1,6)*(data(miss_1,5)<=18)+para1(1,7)*(data(miss_1,5)>50),para1(1,8)*log(5*2.^(0:9))));
prob_temp = bsxfun(@times,sum(para1(1,3)*hazard_comm(miss_1,:),2),sup_temp);
y_prob=binopdf(y_temp,1,1-exp(-prob_temp));   
para2_temp1 = para2(1,1:10);
para2_temp2 = para2(1,11:20);
at2_prob = bsxfun(@times,data(miss_1,5)<=18,para2_temp1(log2(fold)+4)) + bsxfun(@times,data(miss_1,5)>18,para2_temp2(log2(fold)+4));  
prob_matrix = at1_prob.*at2_prob.*y_prob;
titer_mat = 0:9;
titer_cell = mat2cell(titer_mat(ones(sum(miss_1==1),1),:),ones(sum(miss_1==1),1));
ones_cell = mat2cell(ones(sum(miss_1==1),1),ones(sum(miss_1==1),1));
prob_cell= mat2cell(prob_matrix,ones(sum(miss_1==1),1));
data(miss_1,13) = 5*2.^cellfun(@randsample,titer_cell,ones_cell,ones_cell,prob_cell);

% for impute missing_2 & missing_1_2
miss_2_com = miss_2 | miss_1_2;
fold = min(2^6,max(2^-3,bsxfun(@rdivide,5*2.^(0:9),data(miss_2_com,13))));
y_temp = fold>=4;
prob_temp = 1-exp(-sum(para1(1,3)*hazard_comm(miss_2_com,:),2).*exp(para1(1,6)*(data(miss_2_com,5)<=18)+para1(1,7)*(data(miss_2_com,5)>50)+para1(1,8)*log(data(miss_2_com,13))));
y_prob=binopdf(y_temp,1,prob_temp(:,ones(10,1))); 
para2_temp1 = para2(1,1:10);
para2_temp2 = para2(1,11:20);
at2_prob = bsxfun(@times,data(miss_2_com,5)<=18,para2_temp1(log2(fold)+4)) + bsxfun(@times,data(miss_2_com,5)>18,para2_temp2(log2(fold)+4)); 
prob_matrix = at2_prob.*y_prob;
titer_mat = 0:9;
titer_cell = mat2cell(titer_mat(ones(sum(miss_2_com==1),1),:),ones(sum(miss_2_com==1),1));
ones_cell = mat2cell(ones(sum(miss_2_com==1),1),ones(sum(miss_2_com==1),1));
prob_cell= mat2cell(prob_matrix,ones(sum(miss_2_com==1),1));
data(miss_2_com,14) = 5*2.^cellfun(@randsample,titer_cell,ones_cell,ones_cell,prob_cell);

data(data(:,14)./data(:,13) >=4,9) = 1;
data(data(:,14)./data(:,13) < 4,9) = 0;

% impute the infection time based on cover interval
% first simulate the order of infection based on community risk only.
time_cell = mat2cell(z(data(:,9)==1,:),ones(sum(data(:,9)==1),1));
ones_cell = mat2cell(ones(sum(data(:,9)==1),1),ones(sum(data(:,9)==1),1));
prob_cell= mat2cell(hazard_comm(data(:,9)==1,:),ones(sum(data(:,9)==1),1));
data(data(:,9)==1,10) = cellfun(@randsample,time_cell,ones_cell,ones_cell,prob_cell);
% save this as the order to impute the infection
data(data(:,9) == 0,10) = -1;
data_temp = data(:,[1:10 13 14]);
temp = reshape_data(data_temp);
temp = temp(data(:,2) == 0,size(data_temp,2)*(0:max(data(:,3)))+10);
temp = temp(sum(temp,2) >= 0,:);
data(data(:,9) == 1,15) = infection_order(temp);

for u = 2:max(data(:,15));
cond = data(:,9) == 1 & data(:,15) == u & PCR ~= 1;
data_reshape = reshape_data(data);
data1 = data_reshape(cond,:);
hazard1 = hazard_comm(cond,:)*para1(1,3);
time_matrix1 = time_matrix(cond,:);
data(cond,10) = update_infection_time(data1,y,hazard1,time_matrix1,para1(1,:),size(data,2));    
end;

data(data(:,9) == 0,10) = data(data(:,9) == 0,8)+1;


exp_temp = final_size_dist(data);

exp_raw(:,3) = exp_raw(:,3) + exp_temp(:,3)/simN;
end;


exp_sim = exp_raw;
exp_sim(:,3) = 0;

for i = 1:simN;
data_sim = simdata_chisquare(data,para0,para1,ILILAB);
data_sim_final_size = final_size_dist(data_sim);
exp_sim(:,3) = exp_sim(:,3) + data_sim_final_size(:,3)/simN;
end;


chi_square_matrix = (exp_raw(:,3)-exp_sim(:,3)).^2./exp_sim(:,3);

chi_square_stat = sum(chi_square_matrix(exp_sim(:,3) >=5));

p_value = 1-chi2cdf(chi_square_stat,length(chi_square_matrix(exp_sim(:,3) >=5))-1);


result = p_value;
end
