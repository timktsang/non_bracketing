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


data(:,[ 7 9 11 13]) = data(:,[ 7 9 11 13]) - 14;


para1 = [1 0.26 0.15 2.5 0.43 0.51 0 ];
para2 = ones(1,20)/10;
para2_hy = [1.01 1.01 ];
para3 = ones(1,16)/10;     
para3_hy = [1.01 1.01 ];
para4 = [1 2 3 4 5 6 7 10 22 38 2   1 2 3 4 5 6 7 10 22 38 2  ]/100;
para4_hy = [1.01 1.01 ];

cutoff = 509;

tic
[y1 y2 y2_hy y3 y3_hy y4 y4_hy LL pred_obs ] = mcmc(data,11000,para1,para2,para2_hy,para3,para3_hy,para4,para4_hy,cum_ILILAB,cutoff);
toc


startpt = 1001;
endpt = 11000;


csvwrite('mcmc_result_1.csv',y1);
csvwrite('mcmc_result_2.csv',y2);
csvwrite('mcmc_result_2_hy.csv',y2_hy);
csvwrite('mcmc_result_3.csv',y3);
csvwrite('mcmc_result_3_hy.csv',y3_hy);
csvwrite('mcmc_result_4.csv',y4);
csvwrite('mcmc_result_4_hy.csv',y4_hy);
csvwrite('mcmc_result_5.csv',pred_obs);
csvwrite('LL.csv',LL);

z1 = para_summary(y1(startpt:endpt,:),7,2)
z2 = para_summary(y2(startpt:endpt,:),8,5)
z2_hy = para_summary(y2_hy(startpt:endpt,:),2,2)
z3 = para_summary(y3(startpt:endpt,:),8,4)
z3_hy = para_summary(y3_hy(startpt:endpt,:),2,2)
z4 = para_summary(y4(startpt:endpt,:),8,6)
z4_hy = para_summary(y4_hy(startpt:endpt,:),2,2)

z_pred = para_summary(pred_obs(startpt:endpt,:),12,2)


csvwrite('mcmc_summary_1.csv',z1);
csvwrite('mcmc_summary_2.csv',z2);
csvwrite('mcmc_summary_2_hy.csv',z2_hy);
csvwrite('mcmc_summary_3.csv',z3);
csvwrite('mcmc_summary_3_hy.csv',z3_hy);
csvwrite('mcmc_summary_4.csv',z4);
csvwrite('mcmc_summary_4_hy.csv',z4_hy);
