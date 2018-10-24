function[data_prob20]=fliri_klima20(mm, n, data); 

  data_prob20=zeros(n,3);
   for i=1:31
       data_prob20(i,1)=roundn((size(find(data(i,mm-19:mm)<=-2),2)./20)*100,-1);
       data_prob20(i,2)=roundn(min(data(i,mm-19:mm)),-1);
       data_prob20(i,3)=roundn(max(data(i,mm-19:mm)),-1);
   end



