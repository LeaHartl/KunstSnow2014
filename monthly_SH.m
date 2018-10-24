function[SH_longterm_oct_seas,SH_count_oct, SH_longterm_nov_seas,SH_count_nov, SH_longterm_dec_seas,SH_count_dec, SH_longterm_jan_seas,SH_count_jan, ... 
    SH_longterm_feb_seas,SH_count_feb, SH_longterm_mar_seas,SH_count_mar, SH_longterm_apr_seas,SH_count_apr, tab_SH_count]=monthly_SH(months, SH_daily_longterm, mm);
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
   SH_longterm_oct=SH_daily_longterm(longterm_oct);
   SH_longterm_nov=SH_daily_longterm(longterm_nov);
   SH_longterm_dec=SH_daily_longterm(longterm_dec);
   SH_longterm_jan=SH_daily_longterm(longterm_jan);
   SH_longterm_feb=SH_daily_longterm(longterm_feb);
   SH_longterm_mar=SH_daily_longterm(longterm_mar);
   SH_longterm_apr=SH_daily_longterm(longterm_apr);
   
   %differentiate by season
   
   
   SH_longterm_oct_seas=reshape(SH_longterm_oct,31,mm);
   for i=1:mm
       SH_count_oct(i)=size(find(SH_longterm_oct_seas(:,i)>=30),1);
   end
   
   SH_count_oct=SH_count_oct';
   
   SH_longterm_nov_seas=reshape(SH_longterm_nov,30,mm);
   for i=1:mm
       SH_count_nov(i)=size(find(SH_longterm_nov_seas(:,i)>=30),1);
   end
   SH_count_nov=SH_count_nov';
   
   SH_longterm_dec_seas=reshape(SH_longterm_dec,31,mm);
   for i=1:mm
       SH_count_dec(i)=size(find(SH_longterm_dec_seas(:,i)>=30),1);
   end
   SH_count_dec=SH_count_dec';
   
   
    SH_longterm_jan_seas=reshape(SH_longterm_jan,31,mm);
   for i=1:mm
       SH_count_jan(i)=size(find(SH_longterm_jan_seas(:,i)>=30),1);
   end
   SH_count_jan=SH_count_jan';
   
   
     SH_longterm_feb_seas=reshape(SH_longterm_feb,28,mm);
   for i=1:mm
       SH_count_feb(i)=size(find(SH_longterm_feb_seas(:,i)>=30),1);
   end
   SH_count_feb=SH_count_feb';
   
%    days_feb=daysinmonth(jahr_start_2+1:jahr_end,2)';
%   
%   k=cumsum(days_feb);
%    SH_count_feb(1)=size(find(SH_longterm_feb(1:days_feb(1))<=-2),1);
%     for i=2:size(days_feb,1)
%        SH_count_feb(i)=size(find(SH_longterm_feb(k(i-1):k(i))<=-2),1);
%     end
%     
%     SH_count_feb=SH_count_feb';
  
    SH_longterm_mar_seas=reshape(SH_longterm_mar,31,mm);
   for i=1:mm
       SH_count_mar(i)=size(find(SH_longterm_mar_seas(:,i)>=30),1);
   end
   
   SH_count_mar=SH_count_mar';
   
   %april
   %SH_longterm_apr(size(SH_longterm_apr,1)+1:mm*30)=NaN;
   
    SH_longterm_apr_seas=reshape(SH_longterm_apr,30,mm);
   for i=1:mm
       SH_count_apr(i)=size(find(SH_longterm_apr_seas(:,i)>=30),1);
   end
   
   SH_count_apr=SH_count_apr';
   
   
   
      %Tabelle Beschneitage
   zzz=size(SH_count_oct,1);
   tab_SH_count(1,1:zzz)=SH_count_oct;
   tab_SH_count(1,zzz+1)=std(SH_count_oct);
   tab_SH_count(1,zzz+2)=mean(SH_count_oct);
   
   tab_SH_count(2,1:zzz)=SH_count_nov;
   tab_SH_count(2,zzz+1)=std(SH_count_nov);
   tab_SH_count(2,zzz+2)=mean(SH_count_nov);
   
   tab_SH_count(3,1:zzz)=SH_count_dec;
   tab_SH_count(3,zzz+1)=std(SH_count_dec);
   tab_SH_count(3,zzz+2)=mean(SH_count_dec);
   
   tab_SH_count(4,1:zzz)=SH_count_jan;
   tab_SH_count(4,zzz+1)=std(SH_count_jan);
   tab_SH_count(4,zzz+2)=mean(SH_count_jan);
   
   tab_SH_count(5,1:zzz)=SH_count_feb;
   tab_SH_count(5,zzz+1)=std(SH_count_feb);
   tab_SH_count(5,zzz+2)=mean(SH_count_feb);
   
   tab_SH_count(6,1:zzz)=SH_count_mar;
   tab_SH_count(6,zzz+1)=std(SH_count_mar);
   tab_SH_count(6,zzz+2)=mean(SH_count_mar);
   
   tab_SH_count(7,1:zzz)=SH_count_apr;
   tab_SH_count(7,zzz+1)=std(SH_count_apr);
   tab_SH_count(7,zzz+2)=mean(SH_count_apr);