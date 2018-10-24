function[tw_longterm_oct_seas,tw_count_oct, tw_longterm_nov_seas,tw_count_nov, tw_longterm_dec_seas,tw_count_dec, tw_longterm_jan_seas,tw_count_jan, ... 
    tw_longterm_feb_seas,tw_count_feb, tw_longterm_mar_seas,tw_count_mar, tw_longterm_apr_seas,tw_count_apr, tab_tw_count,tw_longterm_oct, tw_longterm_nov,...
    tw_longterm_dec, tw_longterm_jan, tw_longterm_feb, tw_longterm_mar, tw_longterm_apr]=monthly_tw(months, tw_daily_longterm, mm);
%Same calculation for months
   %find months in date variable
   
   longterm_oct=find(months==10);
   longterm_nov=find(months==11);
   longterm_dec=find(months==12);
   longterm_jan=find(months==1);
   longterm_feb=find(months==2);
   longterm_mar=find(months==3);
   longterm_apr=find(months==4);
   
   %Wet-bulb Temp per months
   tw_longterm_oct=tw_daily_longterm(longterm_oct);
   tw_longterm_nov=tw_daily_longterm(longterm_nov);
   tw_longterm_dec=tw_daily_longterm(longterm_dec);
   tw_longterm_jan=tw_daily_longterm(longterm_jan);
   tw_longterm_feb=tw_daily_longterm(longterm_feb);
   tw_longterm_mar=tw_daily_longterm(longterm_mar);
   tw_longterm_apr=tw_daily_longterm(longterm_apr);
   
   %differentiate by season
   
   
   tw_longterm_oct_seas=reshape(tw_longterm_oct,31,mm);
   for i=1:mm
       tw_count_oct(i)=size(find(tw_longterm_oct_seas(:,i)<=-2),1);
   end
   
   for i=1:mm;
    if length(  find(isnan(tw_longterm_oct_seas(:,i))))   >15;
        tw_count_oct(i)=NaN;
    end
end
   
   tw_count_oct=tw_count_oct';
   
   tw_longterm_nov_seas=reshape(tw_longterm_nov,30,mm);
   for i=1:mm
       tw_count_nov(i)=size(find(tw_longterm_nov_seas(:,i)<=-2),1);
   end
   
   
   for i=1:mm;
    if length(  find(isnan(tw_longterm_nov_seas(:,i))))   >15;
        tw_count_nov(i)=NaN;
    end
end
   
   tw_count_nov=tw_count_nov';
   
   tw_longterm_dec_seas=reshape(tw_longterm_dec,31,mm);
   for i=1:mm
       tw_count_dec(i)=size(find(tw_longterm_dec_seas(:,i)<=-2),1);
   end
   
   for i=1:mm;
    if length(  find(isnan(tw_longterm_dec_seas(:,i))))   >15;
        tw_count_dec(i)=NaN;
    end
end
   
   tw_count_dec=tw_count_dec';
   
   
    tw_longterm_jan_seas=reshape(tw_longterm_jan,31,mm);
   for i=1:mm
       tw_count_jan(i)=size(find(tw_longterm_jan_seas(:,i)<=-2),1);
   end
   
   for i=1:mm;
    if length(  find(isnan(tw_longterm_jan_seas(:,i))))   >15;
        tw_count_jan(i)=NaN;
    end
end
   tw_count_jan=tw_count_jan';
   
   
     tw_longterm_feb_seas=reshape(tw_longterm_feb,28,mm);
   for i=1:mm
       tw_count_feb(i)=size(find(tw_longterm_feb_seas(:,i)<=-2),1);
   end
   
   for i=1:mm;
    if length(  find(isnan(tw_longterm_feb_seas(:,i))))   >15;
        tw_count_feb(i)=NaN;
    end
end
   
   tw_count_feb=tw_count_feb';
   
%    days_feb=daysinmonth(jahr_start_2+1:jahr_end,2)';
%   
%   k=cumsum(days_feb);
%    tw_count_feb(1)=size(find(tw_longterm_feb(1:days_feb(1))<=-2),1);
%     for i=2:size(days_feb,1)
%        tw_count_feb(i)=size(find(tw_longterm_feb(k(i-1):k(i))<=-2),1);
%     end
%     
%     tw_count_feb=tw_count_feb';
  
    tw_longterm_mar_seas=reshape(tw_longterm_mar,31,mm);
   for i=1:mm
       tw_count_mar(i)=size(find(tw_longterm_mar_seas(:,i)<=-2),1);
   end
   
   for i=1:mm;
    if length(  find(isnan(tw_longterm_mar_seas(:,i))))   >15;
        tw_count_mar(i)=NaN;
    end
end
   
   tw_count_mar=tw_count_mar';
   
   %april
   %tw_longterm_apr(size(tw_longterm_apr,1)+1:mm*30)=NaN;
   
    tw_longterm_apr_seas=reshape(tw_longterm_apr,30,mm);
   for i=1:mm
       tw_count_apr(i)=size(find(tw_longterm_apr_seas(:,i)<=-2),1);
   end
   
   for i=1:mm;
    if length(  find(isnan(tw_longterm_apr_seas(:,i))))   >15;
        tw_count_apr(i)=NaN;
    end
end
   
   tw_count_apr=tw_count_apr';
   
   
   
      %Tabelle Beschneitage
   zzz=size(tw_count_oct,1);
   tab_tw_count(1,1:zzz)=tw_count_oct;
   tab_tw_count(1,zzz+1)=std(tw_count_oct);
   tab_tw_count(1,zzz+2)=mean(tw_count_oct);
   
   tab_tw_count(2,1:zzz)=tw_count_nov;
   tab_tw_count(2,zzz+1)=std(tw_count_nov);
   tab_tw_count(2,zzz+2)=mean(tw_count_nov);
   
   tab_tw_count(3,1:zzz)=tw_count_dec;
   tab_tw_count(3,zzz+1)=std(tw_count_dec);
   tab_tw_count(3,zzz+2)=mean(tw_count_dec);
   
   tab_tw_count(4,1:zzz)=tw_count_jan;
   tab_tw_count(4,zzz+1)=std(tw_count_jan);
   tab_tw_count(4,zzz+2)=mean(tw_count_jan);
   
   tab_tw_count(5,1:zzz)=tw_count_feb;
   tab_tw_count(5,zzz+1)=std(tw_count_feb);
   tab_tw_count(5,zzz+2)=mean(tw_count_feb);
   
   tab_tw_count(6,1:zzz)=tw_count_mar;
   tab_tw_count(6,zzz+1)=std(tw_count_mar);
   tab_tw_count(6,zzz+2)=mean(tw_count_mar);
   
   tab_tw_count(7,1:zzz)=tw_count_apr;
   tab_tw_count(7,zzz+1)=std(tw_count_apr);
   tab_tw_count(7,zzz+2)=mean(tw_count_apr);