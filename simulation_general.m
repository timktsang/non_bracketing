%%% program for real data

rand('state',sum(clock));

ILILAB = load('ILILAB_ph1.csv');
ILILAB = max(ILILAB,10^-10);
ILILAB = (ILILAB')/sum(ILILAB);

cum_ILILAB = cumsum(ILILAB);

data = load('y1_ph1.csv');

%data = data(data(:,10)~=-1&data(:,12)~=-1&data(:,14)~=-1,:);

data(data(:,7) == -1,7) = median(data(data(:,7)~=-1,7));
data(data(:,9) == -1,9) = median(data(data(:,9)~=-1,9));
data(data(:,11) == -1,11) = median(data(data(:,11)~=-1,11));
data(data(:,13) == -1,13) = median(data(data(:,13)~=-1,13));

%data = data(data(:,5) > 18,:);
data = data(data(:,10)~=-1 & data(:,12)~=-1 & data(:,14)~=-1,:);

data(:,[ 7 9 11 13]) = data(:,[ 7 9 11 13]) - 14;


y1 = load('mcmc_result_1.csv');
y2 = load('mcmc_result_2.csv');
y2_hy = load('mcmc_result_2_hy.csv');
y3 = load('mcmc_result_3.csv');
y3_hy = load('mcmc_result_3_hy.csv');
y4 = load('mcmc_result_4.csv');
y4_hy = load('mcmc_result_4_hy.csv');

startpt = 001;
endpt = 1000;

para1 = median(y1(startpt:endpt,:));
para2 = median(y2(startpt:endpt,:));
para2_hy = median(y2_hy(startpt:endpt,:));
para3 = median(y3(startpt:endpt,:));
para3_hy = median(y3_hy(startpt:endpt,:));
para4 = median(y4(startpt:endpt,:));
para4_hy = median(y4_hy(startpt:endpt,:));

cutoff = 495;

para1 = [ 1 1 1 3 3 3 0];
%for i = 1:9;

data_sim=simdata_general(data,para1,para2,para3,para4,cum_ILILAB,cutoff,1000);


tic
[y1 y2 y2_hy y3 y3_hy y4 y4_hy LL pred_obs] = mcmc(data_sim,1000,para1,para2,para2_hy,para3,para3_hy,para4,para4_hy,cum_ILILAB,cutoff);
toc

startpt = 001;
endpt = 1000;

z1 = para_summary(y1(startpt:endpt,:),7,2);
z2 = para_summary(y2(startpt:endpt,:),8,5);
z2_hy = para_summary(y2_hy(startpt:endpt,:),2,2);
z3 = para_summary(y3(startpt:endpt,:),8,5);
z3_hy = para_summary(y3_hy(startpt:endpt,:),2,2);
z4 = para_summary(y4(startpt:endpt,:),9,5);
z4_hy = para_summary(y4_hy(startpt:endpt,:),2,2);

z_pred = para_summary(pred_obs(startpt:endpt,:),12,2);

%csvwrite(['sim' num2str(i) '.csv'],z1);

%end;




% 
% 
% inc_pair = zeros(size(data_sim,1),4);
% inc_pair(:,1) = data_sim(:,8) ~= -1 & data_sim(:,10) ~= -1;
% inc_pair(:,2) = data_sim(:,10) ~= -1 & data_sim(:,12) ~= -1;
% inc_pair(:,3) = data_sim(:,12) ~= -1 & data_sim(:,14) ~= -1;
% inc_pair(:,4) = data_sim(:,10) ~= -1 & data_sim(:,14) ~= -1 & data_sim(:,12) == -1;
% 
% sum(data_sim(:,10)./data_sim(:,8) >= 4 & inc_pair(:,1) & data_sim(:,5) <= 18)
% sum(data_sim(:,10)./data_sim(:,8) >= 4 & inc_pair(:,1) & data_sim(:,5) > 18 & data_sim(:,5) <= 50)
% sum(data_sim(:,10)./data_sim(:,8) >= 4 & inc_pair(:,1) & data_sim(:,5) > 50 )
% 
% 
% 
% sum(data_sim(:,12)./data_sim(:,10) >= 4 & inc_pair(:,2) & data_sim(:,5) <= 18)
% sum(data_sim(:,12)./data_sim(:,10) >= 4 & inc_pair(:,2) & data_sim(:,5) > 18 & data_sim(:,5) <= 50)
% sum(data_sim(:,12)./data_sim(:,10) >= 4 & inc_pair(:,2) & data_sim(:,5) > 50 )
% 
% 
% sum(data_sim(:,14)./data_sim(:,12) >= 4 & inc_pair(:,3) & data_sim(:,5) <= 18)
% sum(data_sim(:,14)./data_sim(:,12) >= 4 & inc_pair(:,3) & data_sim(:,5) > 18 & data_sim(:,5) <= 50)
% sum(data_sim(:,14)./data_sim(:,12) >= 4 & inc_pair(:,3) & data_sim(:,5) > 50 )
% 
% sum(data_sim(:,14)./data_sim(:,10) >= 4 & inc_pair(:,4) & data_sim(:,5) <= 18)
% sum(data_sim(:,14)./data_sim(:,10) >= 4 & inc_pair(:,4) & data_sim(:,5) > 18 & data_sim(:,5) <= 50)
% sum(data_sim(:,14)./data_sim(:,10) >= 4 & inc_pair(:,4) & data_sim(:,5) > 50 )




