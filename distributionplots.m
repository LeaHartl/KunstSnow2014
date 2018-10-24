function [] = distributionplots(location, yyyy, mon, dd, tw_mean_d, t_mean_d,rel_mean_d, press, tw_longterm_oct, tw_longterm_nov, tw_longterm_dec, tw_longterm_jan, tw_longterm_feb, tw_longterm_mar, tw_longterm_apr, U)

AA=find(yyyy==1974 & mon==10 & dd==1);
AAA=find(yyyy==1994 & mon==4 & dd==30);

yyyy1=yyyy(AA:AAA);
dd1=dd(AA:AAA);
mon1=mon(AA:AAA);
tw_mean_d1=tw_mean_d(AA:AAA);
t_mean_d1=t_mean_d(AA:AAA);
rel_mean_d1=rel_mean_d(AA:AAA);

oct1= find(mon1==10);
nov1= find(mon1==11);
dec1= find(mon1==12);
jan1= find(mon1==1);
feb1= find(mon1==2);
mar1= find(mon1==3);
apr1= find(mon1==4);

AA2=find(yyyy==1994 & mon==10 & dd==1);
AAA2=find(yyyy==2014 & mon==4 & dd==30);

yyyy2=yyyy(AA2:AAA2);
dd2=dd(AA2:AAA2);
mon2=mon(AA2:AAA2);
tw_mean_d2=tw_mean_d(AA2:AAA2);
t_mean_d2=t_mean_d(AA2:AAA2);
rel_mean_d2=rel_mean_d(AA2:AAA2);

oct2= find(mon2==10);
nov2= find(mon2==11);
dec2= find(mon2==12);
jan2= find(mon2==1);
feb2= find(mon2==2);
mar2= find(mon2==3);
apr2= find(mon2==4);

ano_oct=nanmean(tw_mean_d2(oct2))-tw_mean_d2(oct2);
ano_nov=nanmean(tw_mean_d2(nov2))-tw_mean_d2(nov2);
ano_dec=nanmean(tw_mean_d2(dec2))-tw_mean_d2(dec2);
ano_jan=nanmean(tw_mean_d2(jan2))-tw_mean_d2(jan2);
ano_feb=nanmean(tw_mean_d2(feb2))-tw_mean_d2(feb2);
ano_mar=nanmean(tw_mean_d2(mar2))-tw_mean_d2(mar2);
ano_apr=nanmean(tw_mean_d2(apr2))-tw_mean_d2(apr2);


% % Distribution of TW

    [foct1, xioct1]=ksdensity(tw_mean_d1(oct1));
    tw_count_oct1=size(find(tw_mean_d1(oct1)<=-2),1);
    
    [foct2, xioct2]=ksdensity(tw_mean_d2(oct2));
    tw_count_oct2=size(find(tw_mean_d2(oct2)<=-2),1);
    
    [tw_oct3]=wetbulb(t_mean_d2(oct2)+1, rel_mean_d2(oct2), press);
    [foct3, xioct3]=ksdensity(tw_oct3);
    tw_count_oct3=size(find((tw_oct3)<=-2),1);
    
    [tw_oct4]=wetbulb(t_mean_d2(oct2)+1.8, rel_mean_d2(oct2), press);
    [foct4, xioct4]=ksdensity(tw_oct4);
    tw_count_oct4=size(find((tw_oct4)<=-2),1);
%-----------------nov-----------    
    [fnov1, xinov1]=ksdensity(tw_mean_d1(nov1));
    tw_count_nov1=size(find(tw_mean_d1(nov1)<=-2),1);
    
    [fnov2, xinov2]=ksdensity(tw_mean_d2(nov2));
    tw_count_nov2=size(find(tw_mean_d2(nov2)<=-2),1);

    [tw_nov3]=wetbulb(t_mean_d2(nov2)+1, rel_mean_d2(nov2), press);
    [fnov3, xinov3]=ksdensity(tw_nov3);
    tw_count_nov3=size(find((tw_nov3)<=-2),1);
    
    [tw_nov4]=wetbulb(t_mean_d2(nov2)+1.8, rel_mean_d2(nov2), press);
    [fnov4, xinov4]=ksdensity(tw_nov4);
    tw_count_nov4=size(find((tw_nov4)<=-2),1);  
%------------------------dec-----------
    
    [fdec1, xidec1]=ksdensity(tw_mean_d1(dec1));
    tw_count_dec1=size(find(tw_mean_d1(dec1)<=-2),1);
    
    [fdec2, xidec2]=ksdensity(tw_mean_d2(dec2));
    tw_count_dec2=size(find(tw_mean_d2(dec2)<=-2),1);
    
    [tw_dec3]=wetbulb(t_mean_d2(dec2)+1, rel_mean_d2(dec2), press);
    [fdec3, xidec3]=ksdensity(tw_dec3);
    tw_count_dec3=size(find((tw_dec3)<=-2),1);
    
    [tw_dec4]=wetbulb(t_mean_d2(dec2)+1.8, rel_mean_d2(dec2), press);
    [fdec4, xidec4]=ksdensity(tw_dec4);
    tw_count_dec4=size(find((tw_dec4)<=-2),1); 
            
%---------------jan-----------------------    
    
    
    [fjan1, xijan1]=ksdensity(tw_mean_d1(jan1));
    tw_count_jan1=size(find(tw_mean_d1(jan1)<=-2),1);
    
    [fjan2, xijan2]=ksdensity(tw_mean_d2(jan2));
    tw_count_jan2=size(find(tw_mean_d2(jan2)<=-2),1);
    
    [tw_jan3]=wetbulb(t_mean_d2(jan2)+1, rel_mean_d2(jan2), press);
    [fjan3, xijan3]=ksdensity(tw_jan3);
    tw_count_jan3=size(find((tw_jan3)<=-2),1);
    
    [tw_jan4]=wetbulb(t_mean_d2(jan2)+1.8, rel_mean_d2(jan2), press);
    [fjan4, xijan4]=ksdensity(tw_jan4);
    tw_count_jan4=size(find((tw_jan4)<=-2),1); 


%---------------feb-----------------------        
    
    [ffeb1, xifeb1]=ksdensity(tw_mean_d1(feb1));
    tw_count_feb1=size(find(tw_mean_d1(feb1)<=-2),1);

    [ffeb2, xifeb2]=ksdensity(tw_mean_d2(feb2));
    tw_count_feb2=size(find(tw_mean_d2(feb2)<=-2),1);
    
    [tw_feb3]=wetbulb(t_mean_d2(feb2)+1, rel_mean_d2(feb2), press);
    [ffeb3, xifeb3]=ksdensity(tw_feb3);
    tw_count_feb3=size(find((tw_feb3)<=-2),1);
    
    [tw_feb4]=wetbulb(t_mean_d2(feb2)+1.8, rel_mean_d2(feb2), press);
    [ffeb4, xifeb4]=ksdensity(tw_feb4);
    tw_count_feb4=size(find((tw_feb4)<=-2),1); 
%-------------------mar----------------------
    
    [fmar1, ximar1]=ksdensity(tw_mean_d1(mar1));
    tw_count_mar1=size(find(tw_mean_d1(mar1)<=-2),1);

    [fmar2, ximar2]=ksdensity(tw_mean_d2(mar2));
    tw_count_mar2=size(find(tw_mean_d2(mar2)<=-2),1);
    
    [tw_mar3]=wetbulb(t_mean_d2(mar2)+1, rel_mean_d2(mar2), press);
    [fmar3, ximar3]=ksdensity(tw_mar3);
    tw_count_mar3=size(find((tw_mar3)<=-2),1);
    
    [tw_mar4]=wetbulb(t_mean_d2(mar2)+1.8, rel_mean_d2(mar2), press);
    [fmar4, ximar4]=ksdensity(tw_mar4);
    tw_count_mar4=size(find((tw_mar4)<=-2),1); 
%---------------------apr-------------------
    
    [fapr1, xiapr1]=ksdensity(tw_mean_d1(apr1));
    tw_count_apr1=size(find(tw_mean_d1(apr1)<=-2),1);
    
    [fapr2, xiapr2]=ksdensity(tw_mean_d2(apr2));
    tw_count_apr2=size(find(tw_mean_d2(apr2)<=-2),1);
    
    [tw_apr3]=wetbulb(t_mean_d2(apr2)+1, rel_mean_d2(apr2), press);
    [fapr3, xiapr3]=ksdensity(tw_apr3);
    tw_count_apr3=size(find((tw_apr3)<=-2),1);
    
    [tw_apr4]=wetbulb(t_mean_d2(apr2)+1.8, rel_mean_d2(apr2), press);
    [fapr4, xiapr4]=ksdensity(tw_apr4);
    tw_count_apr4=size(find((tw_apr4)<=-2),1); 
%------------------------------------------------------------------

 data_all(1,1)=tw_count_oct1/20;
data_all(2,1)=tw_count_oct2/20;
data_all(3,1)=tw_count_oct3/20;
data_all(4,1)=tw_count_oct4/20;

 data_all(1,2)=tw_count_nov1/20;
data_all(2,2)=tw_count_nov2/20;
data_all(3,2)=tw_count_nov3/20;
data_all(4,2)=tw_count_nov4/20;

 data_all(1,3)=tw_count_dec1/20;
data_all(2,3)=tw_count_dec2/20;
data_all(3,3)=tw_count_dec3/20;
data_all(4,3)=tw_count_dec4/20;

 data_all(1,4)=tw_count_jan1/20;
data_all(2,4)=tw_count_jan2/20;
data_all(3,4)=tw_count_jan3/20;
data_all(4,4)=tw_count_jan4/20;

 data_all(1,5)=tw_count_feb1/20;
data_all(2,5)=tw_count_feb2/20;
data_all(3,5)=tw_count_feb3/20;
data_all(4,5)=tw_count_feb4/20;

 data_all(1,6)=tw_count_mar1/20;
data_all(2,6)=tw_count_mar2/20;
data_all(3,6)=tw_count_mar3/20;
data_all(4,6)=tw_count_mar4/20;

 data_all(1,7)=tw_count_apr1/20;
data_all(2,7)=tw_count_apr2/20;
data_all(3,7)=tw_count_apr3/20;
data_all(4,7)=tw_count_apr4/20;

filename2=['tabelle.txt'];
dlmwrite(filename2,data_all,'delimiter', ' ' , 'newline', 'pc');
    
    
   m7= figure('name','Distribution of Tw')
    set(gca,'fontsize',12)

    subplot(4,2,1)
    h=plot(xioct1,foct1, 'g',xioct2, foct2, 'b', xioct3, foct3, 'k--',  xioct4, foct4, 'r--')%
    %h=plot(xioct2, foct2, 'b', xioct3, foct3, 'k--',  xioct4, foct4, 'r--')%
    set(h, 'linewidth', 2);
    axis([(min(xijan2)-0.5) (max(xioct4)+0.5) 0 0.15])
    hline =line([-2, -2], ylim, 'Color', 'k');
    title('Okt.');
    xlabel('Feuchttemperatur (°C)');
 
    
    subplot(4,2,2)
    h=plot(xinov1,fnov1, 'g', xinov2, fnov2, 'b', xinov3, fnov3, 'k--',  xinov4, fnov4, 'r--')%
    %h=plot(xinov2, fnov2, 'b', xinov3, fnov3, 'k--',  xinov4, fnov4, 'r--')%
    set(h, 'linewidth', 2);
    axis([(min(xijan2)-0.5) (max(xioct4)+0.5) 0 0.15])
    hline =line([-2, -2], ylim, 'Color', 'k');
    title('Nov.');
    xlabel('Feuchttemperatur (°C)');

    subplot(4,2,3)
    h=plot(xidec1,fdec1, 'g', xidec2, fdec2, 'b', xidec3, fdec3, 'k--',  xidec4, fdec4, 'r--')%
    %h=plot(xidec2, fdec2, 'b', xidec3, fdec3, 'k--',  xidec4, fdec4, 'r--')%
    set(h, 'linewidth', 2);
    axis([(min(xijan2)-0.5) (max(xioct4)+0.5) 0 0.15])
    hline =line([-2, -2], ylim, 'Color', 'k');
    title('Dez.');
    xlabel('Feuchttemperatur (°C)');
  
    
    subplot(4,2,4)
    h=plot(xijan1,fjan1, 'g', xijan2, fjan2, 'b', xijan3, fjan3, 'k--',  xijan4, fjan4, 'r--')%
    %h=plot(xijan2, fjan2, 'b', xijan3, fjan3, 'k--',  xijan4, fjan4, 'r--')%
    set(h, 'linewidth', 2);
    axis([(min(xijan2)-0.5) (max(xioct4)+0.5) 0 0.15])
    hline =line([-2, -2], ylim, 'Color', 'k');
    title('Jan.');
    xlabel('Feuchttemperatur (°C)');

    
    subplot(4,2,5)
    h=plot(xifeb1,ffeb1, 'g', xifeb2, ffeb2, 'b', xifeb3, ffeb3, 'k--',  xifeb4, ffeb4, 'r--')%
   % h=plot(xifeb2, ffeb2, 'b', xifeb3, ffeb3, 'k--',  xifeb4, ffeb4, 'r--')%
    set(h, 'linewidth', 2);
    axis([(min(xijan2)-0.5) (max(xioct4)+0.5) 0 0.15])
    hline =line([-2, -2], ylim, 'Color', 'k');
    title('Feb.');
    xlabel('Feuchttemperatur (°C)');
  
    
    subplot(4,2,6)
    h=plot(ximar1,fmar1, 'g', ximar2, fmar2, 'b', ximar3, fmar3, 'k--',  ximar4, fmar4, 'r--')%
   %h=plot(ximar2, fmar2, 'b', ximar3, fmar3, 'k--',  ximar4, fmar4, 'r--')%
   set(h, 'linewidth', 2);
    axis([(min(xijan2)-0.5) (max(xioct4)+0.5) 0 0.15])
    hline =line([-2, -2], ylim, 'Color', 'k');
    title('Mar.');
    xlabel('Feuchttemperatur (°C)');
    %ylabel('relative Häufigkeit');
    
    subplot(4,2,7)
    h=plot(xiapr1,fapr1, 'g', xiapr2, fapr2, 'b', xiapr3, fapr3, 'k--',  xiapr4, fapr4, 'r--')%
    %h=plot(xiapr2, fapr2, 'b', xiapr3, fapr3, 'k--',  xiapr4, fapr4, 'r--.')%
    set(h, 'linewidth', 2);
    axis([(min(xijan2)-0.5) (max(xioct4)+0.5) 0 0.15])
    hline =line([-2, -2], ylim, 'Color', 'k');
    title('Apr.');
    xlabel('Feuchttemperatur (°C)');
    %ylabel('relative Häufigkeit');
    
    legend('1974-94', '1994-2014', 'Zeitraum 1994-2014 plus mögliche Klimaänderung bis 2030', 'Zeitraum 1994-2014 plus mögliche Klimaänderung bis 2050');
    %legend( '1994-2014', '2030', '2050');
saveas(m7,  'fig_dist.jpg');
saveas(m7,  'fig_dist.fig');
    
    %----------------------------------------------------------
%     %----------------------------------------------------------
%     
%     file1=fopen('braun_zukunft.txt');
% 
% C=textscan(file1, '%f %f', 'headerlines',0);
% fclose(file1)
% 
% future_time=C{1};
% tw_fut=C{2};
% a=  size(tw_fut) 
% 
% [yy_f, mon_f, dd_f]=datevec(future_time);
% 
% 
% BB=find(yy_f==2029 & mon_f==10 & dd_f==1)
% BBB=find(yy_f==2049 & mon_f==4 & dd_f==30)
% 
% yy_f=yy_f(BB:BBB);
% dd_f=dd_f(BB:BBB);
% mon_f=mon_f(BB:BBB);
% tw_fut=tw_fut(BB:BBB);
% 
% 
% 
% oct_f= find(mon_f==10);
% nov_f= find(mon_f==11);
% dec_f= find(mon_f==12);
% jan_f= find(mon_f==1);
% feb_f= find(mon_f==2);
% mar_f= find(mon_f==3);
% apr_f= find(mon_f==4);
% 
% 
% % % Distribution of TW
% nanmean(tw_mean_d2(oct2))
% mean(tw_fut(oct_f))
% 
% 
%     [foct_f, xioctf]=ksdensity(tw_fut(oct_f)-ano_oct);
%     tw_count_oct_f=size(find((tw_fut(oct_f)-ano_oct)<=-2),1);
% %-----------------nov-----------    
%     [fnov_f, xinovf]=ksdensity(tw_fut(nov_f));
%     tw_count_nov_f=size(find((tw_fut(nov_f))<=-2),1);  
% %------------------------dec-----------
%     [fdec_f, xidecf]=ksdensity(tw_fut(dec_f));
%     tw_count_dec_f=size(find((tw_fut(dec_f))<=-2),1);
% %---------------jan-----------------------    
%     [fjan_f, xijanf]=ksdensity(tw_fut(jan_f));
%     tw_count_jan_f=size(find((tw_fut(jan_f))<=-2),1); 
% %---------------feb-----------------------        
%     [ffeb_f, xifebf]=ksdensity(tw_fut(feb_f));
%     tw_count_feb_f=size(find((tw_fut(feb_f))<=-2),1);
% %-------------------mar----------------------
%     [fmar_f, ximarf]=ksdensity(tw_fut(mar_f));
%     tw_count_mar_f=size(find((tw_fut(mar_f))<=-2),1);
% %---------------------apr-------------------
%     [fapr_f, xiaprf]=ksdensity(tw_fut(apr_f));
%     tw_count_apr_f=size(find((tw_fut(apr_f))<=-2),1);
%     
    
% figure('name','Distribution of Tw IPCC')
%     set(gca,'fontsize',12)
%    plot(xioct1,foct1, 'r', xioct2, foct2, 'b', xioctf, foct_f, 'k-.')%  
% %     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Okt.');
%     
%     subplot(4,2,1)
%     plot(xioct1,foct1, 'r', xioct2, foct2, 'b', xioctf, foct_f, 'k-.')%  
% %     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Okt.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,2)
%     plot(xinov1,fnov1, 'r', xinov2, fnov2, 'b', xinovf, fnov_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Nov.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,3)
%     plot(xidec1,fdec1, 'r', xidec2, fdec2, 'b', xidecf, fdec_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Dez.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,4)
%     plot(xijan1,fjan1, 'r', xijan2, fjan2, 'b', xijanf, fjan_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Jan.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,5)
%     plot(xifeb1,ffeb1, 'r', xifeb2, ffeb2, 'b', xifebf, ffeb_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Feb.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,6)
%     plot(ximar1,fmar1, 'r', ximar2, fmar2, 'b', ximarf, fmar_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Mar.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,7)
%     plot(xiapr1,fapr1, 'r', xiapr2, fapr2, 'b', xiaprf, fapr_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Apr.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     legend('1974-94', '1994-2014', '2030-2050');  