function f = mean_nan(data,dim)
% temp = ~isnan(data);
% data(isnan(data))=0;
% temp1 = sum(temp,dim);temp1(temp1==0)=1;
% f = sum(data,dim)./temp1;


data(isnan(data))=0;%zero padding
f = mean(data,dim);