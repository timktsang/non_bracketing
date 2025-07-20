% summary of parameter
function result = para_summary(para,k,l)
[nrow ncol ] = size(para);
mat = zeros(ncol,4);
for i = 1:ncol;
mat(i,:) = [ median(para(:,i)) quantile(para(:,i),0.025) quantile(para(:,i),0.975) sum( para(1:(length(para(:,i))-1),i) ~= para(2:length(para(:,i)),i))/(length(para)-1)];   
subplot(k,l,i),plot(para(:,i));
subplot(k,l,i+ncol),hist(para(:,i));
end;
result = mat;
end
