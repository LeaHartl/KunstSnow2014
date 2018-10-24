function[data_prob] = fliri_abb_paper(mm, n, data); 
%nn=number of years
%n=number of days in month. 

   
     data_prob=zeros(n,3);
   for i=1:n
       data_prob(i,1)=roundn((size(find(data(i,mm-19:mm)<=-2),2)./20)*100,-1);
       data_prob(i,2)=roundn((size(find(data(i,mm-19:mm)<=-3),2)./20)*100,-1);
       data_prob(i,3)=roundn((size(find(data(i,mm-19:mm)<=-4),2)./20)*100,-1);

   end


