function[data_prob, data_prob20, data_prob20_2] = fliri_klima(mm, n, data); 
%nn=number of years
%n=number of days in month. 

data_prob=zeros(n,4);
 
   for i=1:n
       data_prob(i,1)=roundn((size(find(data(i,:)<=-2),2)./mm)*100,-1);
       data_prob(i,2)=roundn(min(data(i,:)),-1);
       data_prob(i,3)=roundn(max(data(i,:)),-1);
       data_prob(i,4)=roundn(nanmean(data(i,:)),-1);
   end
   
  data_prob20=zeros(n,4);
   for i=1:n
       data_prob20(i,1)=roundn((size(find(data(i,mm-19:mm)<=-2),2)./20)*100,-1);
       data_prob20(i,2)=roundn(min(data(i,mm-19:mm)),-1);
       data_prob20(i,3)=roundn(max(data(i,mm-19:mm)),-1);
       data_prob20(i,4)=roundn(nanmean(data(i,mm-19:mm)),-1);
   end
%   
      data_prob20_2=zeros(n,4);
      
      %comment from here for Großer Arber and other short time series.
   for i=1:n
       data_prob20_2(i,1)=roundn((size(find(data(i,mm-39:mm-20)<=-2),2)./20)*100,-1);
       data_prob20_2(i,2)=roundn(min(data(i,mm-39:mm-20)),-1);
       data_prob20_2(i,3)=roundn(max(data(i,mm-39:mm-20)),-1);
        data_prob20_2(i,4)=roundn(max(data(i,mm-39:mm-20)),-1);
   end
   

