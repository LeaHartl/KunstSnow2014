%Einlesen der MATLAB-Formatierten DWD Meteodaten aus files
clear all
close all

%***************
%EINSTELLUNGEN
%***************
%==========================================================================
%Tageswerte Schneehöhe 

%ohne februar 29.

location='Zugspitze';
alti=2964; %Seehöhe
 jahr_start=1901;
 jahr_end=2014;
 seasNr=jahr_end-jahr_start;
 dir='snow';
 %--------------------------------------------------!!!!!!!-------------       
     
%==========================================================================
% Daten einlesen
file1=fopen('snow\ZugspitzeSnow.out');
C=textscan(file1, '%f %f %f %f', 'headerlines',0);
fclose(file1);

%**************************************
%PROGRAM SECTION -> nicht verändern!
%**************************************

yyyy=C{1};
mon=C{2};
dd=C{3};
SH=C{4};

%replace 9999 values with NaNs
bad=find(SH==9999);
SH(bad)=NaN;

%delete feb 29
Q=find(mon==2 & dd==29);
yyyy(Q)=[];
mon(Q)=[];
dd(Q)=[];
SH(Q)=[];

%create datevetors and find specific dates to be used later
HH=zeros(size(yyyy));
mins=zeros(size(yyyy));
sec=zeros(size(yyyy));
datumvec=[yyyy, mon, dd,HH, mins, sec];
datestring=datestr(datumvec, 'yyyy-mm-dd');
date_winter=datenum(datestring, 'yyyy-mm-dd');
date_daily=date_winter;
  
pp=find(yyyy==1960 & mon==10 & dd==1);
pp2=find(yyyy==1991 & mon==4 & dd==30);

pp_a=find(unique(yyyy)==1960);
pp2_a=find(unique(yyyy)==1990);

AA=find(yyyy==1951 & mon==10 & dd==1);
AAA=find(yyyy==1952 & mon==4 & dd==30);
U=unique(yyyy);

    %find 1st of January of each year for xticklabeling
    j=1;
    for i=jahr_start:jahr_end-1
      jan_1(j)=datenum(i+1,1,1,0,0,0);
      j=j+1;
    end
     
    %Create date for seasonal mean temps (=15th January)
    %find 15 January of each year for xticklabeling.
    j=1;
    for i=jahr_start:jahr_end-1;
      jan_15(j)=datenum(i+1,1,15,0,0,0);
      j=j+1;
    end
    jan_15=jan_15';
    




%---------------------------Calculate seasonal means-------------------
%----------------------------------------------------------------------
k=212; %Number of days per winter season (Oct-Apr)
   
%Calculate seasonal mean of SH (1st of Oct-30th of April)

SH_winter = reshape(SH, k, seasNr);
  
season_mean_SH=nanmean(SH_winter,1);

%count number of days per season where SH is more than 30.
    for i=1:seasNr;
     count_SH30(i)=length(find(SH_winter(:,i)>=30));
    end
%calculate standard deviation of count variables.    
    count_SH30_std=std(count_SH30);

%Calculate climatological reference period 60/61-90/91 out of seasonal means   
    clim_ref_SH=nanmean(SH_winter(pp:pp2));
     
    clim_ref_count_SH30=round(nanmean(count_SH30(pp_a:pp2_a)));
    
%Abweichung der saisonalen Mittel vom klimatologischen Mittel berechnen
    diff_clim_SH=season_mean_SH-clim_ref_SH;
    diff_clim_count_SH30=count_SH30-clim_ref_count_SH30;
      
%Calculate 10-season running-mean
 
 [diff_clim_SH_10]=running_10(diff_clim_SH); 
 
 [diff_clim_count_SH30_10]=running_10(diff_clim_count_SH30);   
     
%Trend/Rauschverhältnis berechnen, linearen Trend berechnen
      
 %signifikanten Trend finden; Trend ist dann signifikant, wenn Trend/Rauschverhältnis >1.64 (90% Confidence-Level)
jan_15=jan_15' ;

[Tr_max, trend, xxx, y]=trend_rausch(diff_clim_SH, jan_15);

[Tr_max2, trend2, xxxcount2, ycount2]=trend_rausch(count_SH30, jan_15);


%---------------------------section seasonal means over-------------------
%-------------------------------------------------------------------------

%---------------------------Calulate monthly means - same calculations as for seasonal means but for single months -------------------

  [SH_longterm_oct_seas,SH_count_oct, SH_longterm_nov_seas,SH_count_nov, SH_longterm_dec_seas,SH_count_dec, SH_longterm_jan_seas,SH_count_jan, ... 
    SH_longterm_feb_seas,SH_count_feb, SH_longterm_mar_seas,SH_count_mar, SH_longterm_apr_seas,SH_count_apr, tab_SH_count]=monthly_SH(mon, SH, seasNr);

%for monthly means. can also do calculations for number of days with SH
%above 30cm or other value. see commented code blocks below.
   %Calculate climatological reference (1960-90)
   clim_ref_SH_oct=round(nanmean(SH_longterm_oct_seas(pp_a:pp2_a)));
   clim_ref_SH_nov=round(nanmean(SH_longterm_nov_seas(pp_a:pp2_a)));
   clim_ref_SH_dec=round(nanmean(SH_longterm_dec_seas(pp_a:pp2_a)));
   clim_ref_SH_jan=round(nanmean(SH_longterm_jan_seas(pp_a:pp2_a)));
   clim_ref_SH_feb=round(nanmean(SH_longterm_feb_seas(pp_a:pp2_a)));
   clim_ref_SH_mar=round(nanmean(SH_longterm_mar_seas(pp_a:pp2_a)));
   clim_ref_SH_apr=round(nanmean(SH_longterm_apr_seas(pp_a:pp2_a)));
   
   %Abweichung d monatsmittel vom monatsmittel der periode d. langjährigen
   %mittels (1960-1990)
   diff_clim_SH_oct=nanmean(SH_longterm_oct_seas-clim_ref_SH_oct);
   diff_clim_SH_nov=nanmean(SH_longterm_nov_seas-clim_ref_SH_nov);
   diff_clim_SH_dec=nanmean(SH_longterm_dec_seas-clim_ref_SH_dec);
   diff_clim_SH_jan=nanmean(SH_longterm_jan_seas-clim_ref_SH_jan);
   diff_clim_SH_feb=nanmean(SH_longterm_feb_seas-clim_ref_SH_feb);
   diff_clim_SH_mar=nanmean(SH_longterm_mar_seas-clim_ref_SH_mar);
   diff_clim_SH_apr=nanmean(SH_longterm_apr_seas-clim_ref_SH_apr);
   
   SH_oct=nanmean(SH_longterm_oct_seas);
   SH_nov=nanmean(SH_longterm_nov_seas);
   SH_dec=nanmean(SH_longterm_dec_seas);
   SH_jan=nanmean(SH_longterm_jan_seas);
   SH_feb=nanmean(SH_longterm_feb_seas);
   SH_mar=nanmean(SH_longterm_mar_seas);
   SH_apr=nanmean(SH_longterm_apr_seas);


   [diff_clim_SH_oct_10]=running_10(diff_clim_SH_oct); 
   [diff_clim_SH_nov_10]=running_10(diff_clim_SH_nov); 
   [diff_clim_SH_dec_10]=running_10(diff_clim_SH_dec); 
   [diff_clim_SH_jan_10]=running_10(diff_clim_SH_jan); 
   [diff_clim_SH_feb_10]=running_10(diff_clim_SH_feb); 
   [diff_clim_SH_mar_10]=running_10(diff_clim_SH_mar); 
   [diff_clim_SH_apr_10]=running_10(diff_clim_SH_apr); 
   
   %Trend 
%jan_15=jan_15';
[Tr_SH_oct_max, trend_SH_oct, xxx8, y8]=trend_rausch(diff_clim_SH_oct, jan_15);
[Tr_SH_nov_max, trend_SH_nov, xxx9, y9]=trend_rausch(diff_clim_SH_nov, jan_15);
[Tr_SH_dec_max, trend_SH_dec, xxx10, y10]=trend_rausch(diff_clim_SH_dec, jan_15);
[Tr_SH_jan_max, trend_SH_jan, xxx11, y11]=trend_rausch(diff_clim_SH_jan, jan_15);
[Tr_SH_feb_max, trend_SH_feb, xxx12, y12]=trend_rausch(diff_clim_SH_feb, jan_15);
[Tr_SH_mar_max, trend_SH_mar, xxx13, y13]=trend_rausch(diff_clim_SH_mar, jan_15);
[Tr_SH_apr_max, trend_SH_apr, xxx14, y14]=trend_rausch(diff_clim_SH_apr, jan_15);
   
      
 

%for number of days with SH > 30cm per month
%    %Calculate climatological reference (1960-90)
%    clim_ref_count_SH_oct=round(nanmean(SH_count_oct(pp_a:pp2_a)));
%    clim_ref_count_SH_nov=round(nanmean(SH_count_nov(pp_a:pp2_a)));
%    clim_ref_count_SH_dec=round(nanmean(SH_count_dec(pp_a:pp2_a)));
%    clim_ref_count_SH_jan=round(nanmean(SH_count_jan(pp_a:pp2_a)));
%    clim_ref_count_SH_feb=round(nanmean(SH_count_feb(pp_a:pp2_a)));
%    clim_ref_count_SH_mar=round(nanmean(SH_count_mar(pp_a:pp2_a)));
%    clim_ref_count_SH_apr=round(nanmean(SH_count_apr(pp_a:pp2_a)));
%    
%    %Abweichung vom langjährigen Mittel
%    diff_clim_count_SH_oct=SH_count_oct-clim_ref_count_SH_oct;
%    diff_clim_count_SH_nov=SH_count_nov-clim_ref_count_SH_nov;
%    diff_clim_count_SH_dec=SH_count_dec-clim_ref_count_SH_dec;
%    diff_clim_count_SH_jan=SH_count_jan-clim_ref_count_SH_jan;
%    diff_clim_count_SH_feb=SH_count_feb-clim_ref_count_SH_feb;
%    diff_clim_count_SH_mar=SH_count_mar-clim_ref_count_SH_mar;
%    diff_clim_count_SH_apr=SH_count_apr-clim_ref_count_SH_apr;
% 
%    [diff_clim_count_SH_oct_10]=running_10(diff_clim_count_SH_oct); 
%    [diff_clim_count_SH_nov_10]=running_10(diff_clim_count_SH_nov); 
%    [diff_clim_count_SH_dec_10]=running_10(diff_clim_count_SH_dec); 
%    [diff_clim_count_SH_jan_10]=running_10(diff_clim_count_SH_jan); 
%    [diff_clim_count_SH_feb_10]=running_10(diff_clim_count_SH_feb); 
%    [diff_clim_count_SH_mar_10]=running_10(diff_clim_count_SH_mar); 
%    [diff_clim_count_SH_apr_10]=running_10(diff_clim_count_SH_apr); 
 
% %Trend 
% jan_15=jan_15';
% [Tr_count_SH_oct_max, trend_count_SH_oct, xxx8, y8]=trend_rausch(diff_clim_count_SH_oct, jan_15);
% [Tr_count_SH_nov_max, trend_count_SH_nov, xxx9, y9]=trend_rausch(diff_clim_count_SH_nov, jan_15);
% [Tr_count_SH_dec_max, trend_count_SH_dec, xxx10, y10]=trend_rausch(diff_clim_count_SH_dec, jan_15);
% [Tr_count_SH_jan_max, trend_count_SH_jan, xxx11, y11]=trend_rausch(diff_clim_count_SH_jan, jan_15);
% [Tr_count_SH_feb_max, trend_count_SH_feb, xxx12, y12]=trend_rausch(diff_clim_count_SH_feb, jan_15);
% [Tr_count_SH_mar_max, trend_count_SH_mar, xxx13, y13]=trend_rausch(diff_clim_count_SH_mar, jan_15);
% [Tr_count_SH_apr_max, trend_count_SH_apr, xxx14, y14]=trend_rausch(diff_clim_count_SH_apr, jan_15);
%    
%      
% %Fliri-Probability Plot; Probability to have artificial snow conditions on a specific day. Call fliri_klima subfunction. Input: seasNr, number of
%    %days per month, SH_longterm_seas :SH_count_oct_prob
%    
% %Shorter Period (1987/88-2006/07) last 20 seasons to detect changes in probability : SH_count_oct_prob_last20
%    
% [SH_count_oct_prob, SH_count_oct_prob_last20, SH_count_oct_prob_last20_2]=fliri_snow(seasNr, 31, SH_longterm_oct_seas); 
% [SH_count_nov_prob, SH_count_nov_prob_last20, SH_count_nov_prob_last20_2]=fliri_snow(seasNr, 30, SH_longterm_nov_seas); 
% [SH_count_dec_prob, SH_count_dec_prob_last20, SH_count_dec_prob_last20_2]=fliri_snow(seasNr, 31, SH_longterm_dec_seas); 
% [SH_count_jan_prob, SH_count_jan_prob_last20, SH_count_jan_prob_last20_2]=fliri_snow(seasNr, 31, SH_longterm_jan_seas); 
% [SH_count_feb_prob, SH_count_feb_prob_last20, SH_count_feb_prob_last20_2]=fliri_snow(seasNr, 28, SH_longterm_feb_seas); 
% [SH_count_mar_prob, SH_count_mar_prob_last20, SH_count_mar_prob_last20_2]=fliri_snow(seasNr, 31, SH_longterm_mar_seas); 
% [SH_count_apr_prob, SH_count_apr_prob_last20, SH_count_apr_prob_last20_2]=fliri_snow(seasNr, 30, SH_longterm_apr_seas); 
% 
%    
% %   Merge monthly probabilities to yearly plot
%   SH_prob_year(1:31,1:3)=SH_count_oct_prob;
%   SH_prob_year(32:61,1:3)=SH_count_nov_prob;
%   SH_prob_year(62:92,1:3)=SH_count_dec_prob;
%   SH_prob_year(93:123,1:3)=SH_count_jan_prob;
%   SH_prob_year(124:151,1:3)=SH_count_feb_prob;
%   SH_prob_year(152:182,1:3)=SH_count_mar_prob;
%   SH_prob_year(183:212,1:3)=SH_count_apr_prob; 
%   
%  date_year=date_winter(1:212); %dates for one season for plot (takes first season in time series but year is not shown in plot)
% 
%   
%   %Merge monthly probabilities to yearly plot
%   SH_prob_year_last20(1:31,1:3)=SH_count_oct_prob_last20;
%   SH_prob_year_last20(32:61,1:3)=SH_count_nov_prob_last20;
%   SH_prob_year_last20(62:92,1:3)=SH_count_dec_prob_last20;
%   SH_prob_year_last20(93:123,1:3)=SH_count_jan_prob_last20;
%   SH_prob_year_last20(124:151,1:3)=SH_count_feb_prob_last20;
%   SH_prob_year_last20(152:182,1:3)=SH_count_mar_prob_last20;
%   SH_prob_year_last20(183:212,1:3)=SH_count_apr_prob_last20; 
%     
%   %Merge monthly probabilities to yearly plot
%   SH_prob_year_last20_2(1:31,1:3)=SH_count_oct_prob_last20_2;
%   SH_prob_year_last20_2(32:61,1:3)=SH_count_nov_prob_last20_2;
%   SH_prob_year_last20_2(62:92,1:3)=SH_count_dec_prob_last20_2;
%   SH_prob_year_last20_2(93:123,1:3)=SH_count_jan_prob_last20_2;
%   SH_prob_year_last20_2(124:151,1:3)=SH_count_feb_prob_last20_2;
%   SH_prob_year_last20_2(152:182,1:3)=SH_count_mar_prob_last20_2;
%   SH_prob_year_last20_2(183:212,1:3)=SH_count_apr_prob_last20_2; 
%   
  %WRITE OUTPUT FOR TABLES
 
%Write files with number of snowdays per year and per year and month.
%Insert desired file name. 

% jahr=U;
% jahr(end)=[];
% jahr=jahr';
% data(:,1)=jahr;
% data(:,2)=count_SH30; 
% data(:,3)=season_mean_SH;
% 
%  filename1=['mehr als 30cm und mittel',location, '.txt'];
%  dlmwrite(filename1,data,'delimiter', ' ', 'newline', 'pc');
 
% data_all(:,1)=SH_count_oct;
% data_all(:,2)=SH_count_nov;
% data_all(:,3)=SH_count_dec;
% data_all(:,4)=SH_count_jan;
% data_all(:,5)=SH_count_feb;
% data_all(:,6)=SH_count_mar;
% data_all(:,7)=SH_count_apr;
%   
% 
% last20_2=nanmean(SH_prob_year_last20_2)
% last20=nanmean(SH_prob_year_last20)
% tr=trend
% 
% filename2=['mehr als 30cm',location, 'all','.txt'];
% dlmwrite(filename2,data_all,'delimiter', ' ' , 'newline', 'pc');

 
%  
%   %PLOTS
%    
%    Trendanalyse SH
   m= figure('name','Trend Schneehöhe');
   
    subplot(2,1,1)
    set(gca,'fontsize',14)
    plot(jan_15,season_mean_SH,'r','linewidth',2)
    grid on
    axis tight
    ylabel('Schneehöhe (cm)')
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    title(['Saisonale Schneehöhe ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
    %if nofile~=1
%    line([datum(1) datum(1)],[min(season_mean) max(season_mean)],'linewidth',2,'color',[0.7 0.7 0.7])
    %else
    %  line([datum_clim(1) datum_clim(1)],[min(season_mean) max(season_mean)],'linewidth',2,'color',[0.7 0.7 0.7])  
    %end
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    hold on
    bar(jan_15,diff_clim_SH,'b')
    plot(jan_15,diff_clim_SH_10,'k','linewidth',2)
    %if Tr_max>1.64
    size(y)
    size(xxx)
    plot(jan_15,y(:,xxx),'-.r','linewidth',2)
    %end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Anomalie d. Schneehöhe (cm)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH)) round(max(diff_clim_SH))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl+2,['Stabw.: ',num2str(roundn(std(diff_clim_SH),-1)),' cm'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+500,ymittl+1.5,['Klimatol. Mittel (1960-90): ',num2str(roundn(clim_ref_SH,-1)),' cm'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl+2,['T/R.: ',num2str(roundn((Tr_max),-2)),' (',num2str(datestr(jan_15(xxx),'yyyy')),'/',num2str(datestr(jan_15(xxx+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_max>1.64
    text(xpos(1)+15500,ymittl+2,['Linearer Trend: ',num2str(roundn(trend,-1)),'cm/a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',14)
    title(['Abweichung der saisonalen Schneehöhe ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.) vom langjährigen Mittel (1960-90), ',num2str(location)],'fontsize',12,'fontweight','bold')
    %line([jan_15(47)-175 jan_15(47)-175],[min(diff_clim_t) max(diff_clim_t)+2],'linewidth',2,'color',[0.7 0.7 0.7])
    
saveas(m, 'fig1.jpg');
saveas(m, 'fig1.fig');



    
m2=    figure('name','Trend number of days monthly');
    subplot(3,1,1)
    set(gca,'fontsize',14)
    %plot(jan_15,SH_count_oct,'r',jan_15,SH_count_nov,'-.k')
    plot(jan_15,SH_oct,'r',jan_15,SH_nov,'-.k')
    grid off
    axis tight
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    legend('Oct','Nov','location','best')
    xlabel('Jahre')
    %ylabel('Anomalie d. Tagesanzahl mit Tagesmittel Tf \leq -2°C')
    title('Vorsaison (Oktober & November)')
    
    subplot(3,1,2)
    set(gca,'fontsize',14)
    %plot(jan_15,SH_count_dec,'r',jan_15,SH_count_jan,'-.k',jan_15,SH_count_feb,':b')
    plot(jan_15,SH_dec,'r',jan_15,SH_jan,'-.k',jan_15,SH_feb,':b')
    grid off
    axis tight
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    legend('Dez','Jan','Feb','location','best')
    xlabel('Jahre')
    %ylabel('Tagesanzahl mit SH >= 30cm')
    ylabel('Schneehöhe - Monatsmittel (cm)')
    title('Hauptsaison (Dezember - Februar)')
    
    subplot(3,1,3)
    set(gca,'fontsize',14)
%     plot(jan_15,SH_count_mar,'r',jan_15,SH_count_apr,'-.k')
    plot(jan_15,SH_mar,'r',jan_15,SH_apr,'-.k')
    grid off
    axis tight
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    legend('Mar','Apr','location','best')
    xlabel('Jahre')
   % ylabel('Tagesanzahl mit Tagesmittel Tf  \leq -2°C')
    title('Nachsaison (März - April)')
saveas(m2, 'fig6.jpg');
saveas(m2, 'fig6.fig');



    
m3=    figure('name','Anomalie d. Anzahl an Tagen nach Monat: Vorsaison');
%    Trend number of days
    
    subplot(2,1,1)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_count_SH_oct,'b')
    plot(jan_15,diff_clim_count_SH_oct_10,'k','linewidth',2)
    if Tr_count_SH_oct_max>1.64
    plot(jan_15,y8(:,xxx8),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Oktober (Tage)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_oct)) round(max(diff_clim_count_SH_oct))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_oct),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_oct_max,-2)),' (',num2str(datestr(jan_15(xxx8),'yyyy')),'/',num2str(datestr(jan_15(xxx8+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_count_SH_oct_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_oct,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',14)
     title('Vorsaison (Oktober & November) Anomalie d. Anzahl an Tagen mit SH \geq 30cm')
    
    subplot(2,1,2)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_count_SH_nov,'b')
    plot(jan_15,diff_clim_count_SH_nov_10,'k','linewidth',2)
    if Tr_count_SH_nov_max>1.64
    plot(jan_15,y9(:,xxx9),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('November (Tage)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_nov)) round(max(diff_clim_count_SH_nov))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_nov),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_nov_max,-2)),' (',num2str(datestr(jan_15(xxx9),'yyyy')),'/',num2str(datestr(jan_15(xxx9+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_count_SH_nov_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_nov,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    
saveas(m3, 'fig7.jpg');
saveas(m3, 'fig7.fig');
    
%     

%für schneehöhenanomalie, nicht anzahl an tagen. 
m3=    figure('name','Anomalie d. Anzahl an Tagen nach Monat: Vorsaison');
%    Trend number of days
    
    subplot(2,1,1)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_SH_oct,'b')
    plot(jan_15,diff_clim_SH_oct_10,'k','linewidth',2)
    if Tr_SH_oct_max>1.64
    plot(jan_15,y8(:,xxx8),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Oktober (cm)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_oct)) round(max(diff_clim_SH_oct))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Kilm. Mittel.: ',num2str(clim_ref_SH_oct)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_oct_max,-2)),' (',num2str(datestr(jan_15(xxx8),'yyyy')),'/',num2str(datestr(jan_15(xxx8+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_SH_oct_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_oct,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',14)
     title('Vorsaison (Oktober & November) Anomalie d. Schneehöhe im Vergleich zum langjährigen Mittel.')
    
    subplot(2,1,2)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_SH_nov,'b')
    plot(jan_15,diff_clim_SH_nov_10,'k','linewidth',2)
    if Tr_SH_nov_max>1.64
    plot(jan_15,y9(:,xxx9),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('November (cm)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_nov)) round(max(diff_clim_SH_nov))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Klim. Mittel: ',num2str(clim_ref_SH_nov)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_nov_max,-2)),' (',num2str(datestr(jan_15(xxx9),'yyyy')),'/',num2str(datestr(jan_15(xxx9+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_SH_nov_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_nov,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    
saveas(m3, 'fig7.jpg');
saveas(m3, 'fig7.fig');
    




    
    
 m4=   figure('name','Anomalie d. Anzahl an Tagen nach Monat: Hauptsaison');
  %  Trend number of days
    
    subplot(3,1,1)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_count_SH_dec,'b')
    plot(jan_15,diff_clim_count_SH_dec_10,'k','linewidth',2)
    if Tr_count_SH_dec_max>1.64
    plot(jan_15,y10(:,xxx10),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Dezember')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_dec)) round(max(diff_clim_count_SH_dec))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_dec),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_dec_max,-2)),' (',num2str(datestr(jan_15(xxx10),'yyyy')),'/',num2str(datestr(jan_15(xxx10+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_count_SH_dec_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_dec,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
   end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    title('Hauptsaison (Dezember, Januar, Februar) Anomalie d. Anzahl an Tagen mit SH \geq 30cm');
    
    
    subplot(3,1,2)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_count_SH_jan,'b')
    plot(jan_15,diff_clim_count_SH_jan_10,'k','linewidth',2)
   if Tr_count_SH_jan_max>1.64
    plot(jan_15,y11(:,xxx11),'-.r','linewidth',2)
   end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Januar')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_jan)) round(max(diff_clim_count_SH_jan))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_jan),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_jan_max,-2)),' (',num2str(datestr(jan_15(xxx11),'yyyy')),'/',num2str(datestr(jan_15(xxx11+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_count_SH_jan_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_jan,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',14)
        
    subplot(3,1,3)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_count_SH_feb,'b')
    plot(jan_15,diff_clim_count_SH_feb_10,'k','linewidth',2)
    if Tr_count_SH_feb_max>1.64
    plot(jan_15,y12(:,xxx12),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Februar')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_feb)) round(max(diff_clim_count_SH_feb))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_feb),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_feb_max,-2)),' (',num2str(datestr(jan_15(xxx12),'yyyy')),'/',num2str(datestr(jan_15(xxx12+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
   if Tr_count_SH_feb_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_feb,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
   end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    
    saveas(m4, 'fig8.jpg');
saveas(m4, 'fig8.fig');
    
    

 m4=   figure('name','Anomalie d. Anzahl an Tagen nach Monat: Hauptsaison');
  %  Trend number of days
    
    subplot(3,1,1)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_SH_dec,'b')
    plot(jan_15,diff_clim_SH_dec_10,'k','linewidth',2)
    if Tr_SH_dec_max>1.64
    plot(jan_15,y10(:,xxx10),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Dezember (cm)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_dec)) round(max(diff_clim_SH_dec))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Klim. Mittel.: ',num2str(clim_ref_SH_dec)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_dec_max,-2)),' (',num2str(datestr(jan_15(xxx10),'yyyy')),'/',num2str(datestr(jan_15(xxx10+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_SH_dec_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_dec,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
   end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    title('Hauptsaison (Dezember, Januar, Februar) Anomalie d. Schneehöhe im Vergleich zum langjährigen Mittel');
   
    
    subplot(3,1,2)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_SH_jan,'b')
    plot(jan_15,diff_clim_SH_jan_10,'k','linewidth',2)
   if Tr_SH_jan_max>1.64
    plot(jan_15,y11(:,xxx11),'-.r','linewidth',2)
   end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Januar (cm)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_jan)) round(max(diff_clim_SH_jan))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Klim. Mittel.: ',num2str(clim_ref_SH_jan)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_jan_max,-2)),' (',num2str(datestr(jan_15(xxx11),'yyyy')),'/',num2str(datestr(jan_15(xxx11+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_SH_jan_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_jan,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',14)
        
    subplot(3,1,3)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_SH_feb,'b')
    plot(jan_15,diff_clim_SH_feb_10,'k','linewidth',2)
    if Tr_SH_feb_max>1.64
    plot(jan_15,y12(:,xxx12),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Februar (cm)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_feb)) round(max(diff_clim_SH_feb))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Klim. Mittel.: ',num2str(clim_ref_SH_feb)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_feb_max,-2)),' (',num2str(datestr(jan_15(xxx12),'yyyy')),'/',num2str(datestr(jan_15(xxx12+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
   if Tr_SH_feb_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_feb,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
   end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    
    saveas(m4, 'fig8.jpg');
saveas(m4, 'fig8.fig');






    
  m5=  figure('name','Anomalie d. Anzahl an Tagen nach Monat: Nachsaison');
   
    subplot(2,1,1)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_count_SH_mar,'b')
    plot(jan_15,diff_clim_count_SH_mar_10,'k','linewidth',2)
    if Tr_count_SH_mar_max>1.64
    plot(jan_15,y13(:,xxx13),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('März')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_mar)) round(max(diff_clim_count_SH_mar))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_mar),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_mar_max,-2)),' (',num2str(datestr(jan_15(xxx13),'yyyy')),'/',num2str(datestr(jan_15(xxx13+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_count_SH_mar_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_mar,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    title('Nachsaison (März & April)Anomalie d. Anzahl an Tagen mit SH \geq 30cm')
    
    subplot(2,1,2)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_count_SH_apr,'b')
    plot(jan_15,diff_clim_count_SH_apr_10,'k','linewidth',2)
    if Tr_count_SH_apr_max>1.64
    plot(jan_15,y14(:,xxx14),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('April')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_apr)) round(max(diff_clim_count_SH_apr))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_apr),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_apr_max,-2)),' (',num2str(datestr(jan_15(xxx14),'yyyy')),'/',num2str(datestr(jan_15(xxx14+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_count_SH_apr_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_apr,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    
    
    saveas(m5, 'fig9.jpg');
saveas(m5, 'fig9.fig');


    
  m5=  figure('name','Anomalie d. Anzahl an Tagen nach Monat: Nachsaison');
   
    subplot(2,1,1)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_SH_mar,'b')
    plot(jan_15,diff_clim_SH_mar_10,'k','linewidth',2)
    if Tr_SH_mar_max>1.64
    plot(jan_15,y13(:,xxx13),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('März (cm)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_mar)) round(max(diff_clim_SH_mar))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Klim. Mittel: ',num2str(clim_ref_SH_mar)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_mar_max,-2)),' (',num2str(datestr(jan_15(xxx13),'yyyy')),'/',num2str(datestr(jan_15(xxx13+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_SH_mar_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_mar,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    title('Nachsaison (März & April) Anomalie d. Schneehöhe im Vergleich zum langjährigen Mittel')
    
    subplot(2,1,2)
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    bar(jan_15,diff_clim_SH_apr,'b')
    plot(jan_15,diff_clim_SH_apr_10,'k','linewidth',2)
    if Tr_SH_apr_max>1.64
    plot(jan_15,y14(:,xxx14),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('April (cm)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_apr)) round(max(diff_clim_SH_apr))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl,['Klim. Mittel: ',num2str(clim_ref_SH_mar)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_apr_max,-2)),' (',num2str(datestr(jan_15(xxx14),'yyyy')),'/',num2str(datestr(jan_15(xxx14+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_SH_apr_max>1.64
    text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_apr,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',12,'fontweight','bold')
    
    
    saveas(m5, 'fig9.jpg');
saveas(m5, 'fig9.fig');


dateplot=[date_daily(AA):date_daily(AAA-1)];    

 m6=   figure('name','Probability Year');
    [haxes,hline1,hline2]=plotyy(date_year,SH_prob_year(:,1),date_year,SH_prob_year(:,2),'area','plot');
    grid off
    set(hline1,'FaceColor',[0.8 0.8 0.8])
    set(hline2,'Color','b')
    axes(haxes(1))
    set(gca,'fontsize',14)
    ylabel('Wahrscheinlichkeit (%)')
    set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
    set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
    datetick('x','mmm','keeplimits','keepticks')
    xlabel('Monate')
    axis tight
    axes(haxes(2))
    set(gca,'fontsize',14,'YColor','k')
    ylabel('Extremwerte Schneehöhe (cm)')
    set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
    set(gca,'YTick',[-25 -20 -15 -10 -5 0 5])
    datetick('x','mmm','keeplimits','keepticks')
    xlabel('Monate')
    hold on
    plot(date_year,SH_prob_year(:,3),'r')
    title(['Wahrscheinlichkeit f. Schneehöhe \geq 30 cm, und Extremwerte (basierend auf ',num2str(jahr_start),' - ',num2str(jahr_end),' ) ,',num2str(location),' ',num2str(alti),' m']);
    axis tight
    hline=refline(0,-2);
    set(hline, 'Color', 'k', 'LineWidth', 2);
    
    saveas(m6, 'fig10.jpg');
saveas(m6, 'fig10.fig');
    
    
    
  m7=  figure('name','Probability last 20 years');
    [haxes,hline1,hline2]=plotyy(date_year,SH_prob_year_last20(:,1),date_year,SH_prob_year_last20(:,2),'area','plot');
    grid off
    set(hline1,'FaceColor',[0.8 0.8 0.8])
    set(hline2,'Color','b')
    axes(haxes(1))
    set(gca,'fontsize',14)
    ylabel('Wahrscheinlichkeit (%)')
    set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
    set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
    datetick('x','mmm','keeplimits','keepticks')
    xlabel('Monate')
    axis tight
    axes(haxes(2))
    set(gca,'fontsize',14,'YColor','k')
    ylabel('Extremwerte Schneehöhe (cm)')
    set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
    set(gca,'YTick',[-25 -20 -15 -10 -5 0 5])
    datetick('x','mmm','keeplimits','keepticks')
    xlabel('Monate')
    hold on
    plot(date_year,SH_prob_year_last20(:,3),'r')
    title(['Wahrscheinlichkeit f. Tagesmittelwert Schneehöhe \geq 30 cm, und Extremwerte (basierend auf 1993-2014), ', num2str(location),' ',num2str(alti)])
    axis tight
    hline=refline(0,-2);
    set(hline, 'Color', 'k', 'LineWidth', 2);
    
        saveas(m7, 'fig11.jpg');
saveas(m7, 'fig11.fig');
    
    
  m8=  figure('name','Comparison Probabilities last 20 years/20 years before');
    area(date_year,SH_prob_year_last20(:,1),'facecolor','g')
    hold on
    area(date_year,SH_prob_year_last20_2(:,1))
    grid on
    ylabel('Wahrscheinlichkeit')
    xlabel('Monate')
    axis tight
    set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
    datetick('x','mmm','keeplimits','keepticks')
    jahr_20=[num2str(jahr_end-20),'/ ',num2str(jahr_end-1919),'-',num2str(jahr_end-1),'/ ',num2str(jahr_end-2000)];
    jahr_40=[num2str(jahr_end-40),'/ ',num2str(jahr_end-1939),'-',num2str(jahr_end-21),'/ ',num2str(jahr_end-1940)];
    %legend('1987/88-2006/07','1967/68-1986/87')
    legend(jahr_20,jahr_40)
    
        saveas(m8, 'fig12.jpg');
saveas(m8, 'fig12.fig');
% % % % % % end
% % % % 
% % % % toc
% % % % t=toc;
% % % % 
% % % %display(t)
% % % 
% % % 
% % % 
