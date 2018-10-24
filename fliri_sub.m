function[tw_shortterm_hours, tw_count_sum, std_sum, mean_sum]=fliri_sub(tw_shortterm, a, b, n);   

%  if (a== 720 | 744)
 %a=hours per month, n= nr. of rows of snow_hours_seas
 %b=days per month
   tw_shortterm=reshape(tw_shortterm,a,n);
   
   for j=1:n
        k=1;
        for i=1:b
           tw_count_hour(i,j)=size(find(tw_shortterm(k:k+23,j)<=-2),1);
       k=k+24;
       end
   end
   
   tw_shortterm_hours(:,1)=nanmean(tw_count_hour,2);
   tw_shortterm_hours(:,2)=max(tw_count_hour,[],2);
   tw_shortterm_hours(:,3)=min(tw_count_hour,[],2);
   
   tw_count_sum=sum(tw_count_hour);
   std_sum=std(tw_count_sum);
   mean_sum=mean(tw_count_sum);
   
   
% else 
%         days_feb=daysinmonth(jahr_start+1:jahr_end,2)';
%    
% [tw_shortterm_jan_hours, tw_count_jan_sum, std_jan_sum, mean_jan_sum]=fliri_sub(tw_shortterm_jan, days_feb, n);
%    
%    
%     ii=1;
%     k=1;
%     for j=1:n
%         for i=1:days_feb(ii)
%            tw_count_feb_hour(i,j)=size(find(tw_shortterm_feb(k:k+23)<=-2),2);
%        k=k+24;
%        end
%        ii=ii+1;
%    end
%      
% 
% tw_count_feb_hour(29,:)=[];   
% 
% tw_shortterm_feb_hours(:,1)=nanmean(tw_count_feb_hour,2);
% tw_shortterm_feb_hours(:,2)=max(tw_count_feb_hour,[],2);
% tw_shortterm_feb_hours(:,3)=min(tw_count_feb_hour,[],2);
% 
% tw_count_feb_sum=sum(tw_count_feb_hour);
% std_feb_sum=std(tw_count_feb_sum);
% mean_feb_sum=mean(tw_count_feb_sum);
% 
%     end 
% end