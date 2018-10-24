%Einlesen der MATLAB-Formatierten DWD Meteodaten aus files
clear all
close all

%***************
%EINSTELLUNGEN
%***************
%==========================================================================
%Tageswerte Schneehöhe 

%ohne februar 29.

location='Gr. Arber';
alti=1446; %Seehöhe
 jahr_start=1983;
 jahr_end=2014;
 seasNr=jahr_end-jahr_start;
 dir='snow';
 %--------------------------------------------------!!!!!!!-------------       
     
%==========================================================================
% Daten einlesen
file1=fopen('snow\GrArberSnow.out');
C=textscan(file1, '%f %f %f %f', 'headerlines',0);
fclose(file1);

file2=fopen('snow\twprobyearGrosser Arber.out');
D=textscan(file2, '%f %f %f', 'headerlines',0);
fclose(file2);

tw_prob_year(:,1)=D{1};
tw_prob_year(:,2)=D{2};
tw_prob_year(:,3)=D{3};

file3=fopen('snow\twprobyearlast20Grosser Arber.out');
DD=textscan(file3, '%f %f %f', 'headerlines',0);
fclose(file3);

tw_prob_year_last20(:,1)=DD{1};
tw_prob_year_last20(:,2)=DD{2};
tw_prob_year_last20(:,3)=DD{3};

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
  


AA=find(yyyy==1998 & mon==10 & dd==1);
AAA=find(yyyy==1999 & mon==4 & dd==30);
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
size(SH)
k
seasNr
SH_winter = reshape(SH, k, seasNr);
  
season_mean_SH=nanmean(SH_winter,1);


%count number of days per season where SH is more than 30.
    for i=1:seasNr;
     count_SH30(i)=length(find(SH_winter(:,i)>=30));
    end
    
    thirthy=zeros(size(SH_winter));
    abc=find(SH>=30);
    thirty(abc)=1;
    
    
     for i=1:seasNr;
    if length(  find(isnan(SH_winter(:,i))))   >106
        season_mean_SH(i)=NaN;
    end
end
 
     fixNaN=find(isnan(season_mean_SH));

     
% diff_clim_SH(fixNaN)=NaN;
count_SH30(fixNaN)=NaN;
% diff_clim_count_SH30(fixNaN)=NaN;
% diff_clim_count_SH30_10(fixNaN)=NaN;
    
%calculate standard deviation of count variables.    
    count_SH30_std=nanstd(count_SH30);

% %Calculate climatological reference period 60/61-90/91 out of seasonal means   
%     clim_ref_SH=nanmean(SH_winter(pp:pp2));
%      
%     clim_ref_count_SH30=round(nanmean(count_SH30(pp_a:pp2_a)));
%     
% %Abweichung der saisonalen Mittel vom klimatologischen Mittel berechnen
%     diff_clim_SH=season_mean_SH-clim_ref_SH;
%     diff_clim_count_SH30=count_SH30-clim_ref_count_SH30;
%       
% %Calculate 10-season running-mean
%  
%  [diff_clim_SH_10]=running_10(diff_clim_SH); 
%  
%  [diff_clim_count_SH30_10]=running_10(diff_clim_count_SH30);   
     

%Trend/Rauschverhältnis berechnen, linearen Trend berechnen

      
 %signifikanten Trend finden; Trend ist dann signifikant, wenn Trend/Rauschverhältnis >1.64 (90% Confidence-Level)
jan_15=jan_15' ;
% 
% [Tr_max, trend, xxx, y]=trend_rausch(diff_clim_SH, jan_15);
% 
% [Tr_max2, trend2, xxxcount2, ycount2]=trend_rausch(count_SH30, jan_15);


%---------------------------section seasonal means over-------------------
%-------------------------------------------------------------------------

%---------------------------Calulate monthly means - same calculations as for seasonal means but for single months -------------------

  [SH_longterm_oct_seas,SH_count_oct, SH_longterm_nov_seas,SH_count_nov, SH_longterm_dec_seas,SH_count_dec, SH_longterm_jan_seas,SH_count_jan, ... 
    SH_longterm_feb_seas,SH_count_feb, SH_longterm_mar_seas,SH_count_mar, SH_longterm_apr_seas,SH_count_apr, tab_SH_count]=monthly_SH(mon, SH, seasNr);

% %for monthly means. can also do calculations for number of days with SH
% %above 30cm or other value. see commented code blocks below.
%    %Calculate climatological reference (1960-90)
%    clim_ref_SH_oct=round(nanmean(SH_longterm_oct_seas(pp_a:pp2_a)));
%    clim_ref_SH_nov=round(nanmean(SH_longterm_nov_seas(pp_a:pp2_a)));
%    clim_ref_SH_dec=round(nanmean(SH_longterm_dec_seas(pp_a:pp2_a)));
%    clim_ref_SH_jan=round(nanmean(SH_longterm_jan_seas(pp_a:pp2_a)));
%    clim_ref_SH_feb=round(nanmean(SH_longterm_feb_seas(pp_a:pp2_a)));
%    clim_ref_SH_mar=round(nanmean(SH_longterm_mar_seas(pp_a:pp2_a)));
%    clim_ref_SH_apr=round(nanmean(SH_longterm_apr_seas(pp_a:pp2_a)));
   
%    %Abweichung d monatsmittel vom monatsmittel der periode d. langjährigen
%    %mittels (1960-1990)
%    diff_clim_SH_oct=nanmean(SH_longterm_oct_seas-clim_ref_SH_oct);
%    diff_clim_SH_nov=nanmean(SH_longterm_nov_seas-clim_ref_SH_nov);
%    diff_clim_SH_dec=nanmean(SH_longterm_dec_seas-clim_ref_SH_dec);
%    diff_clim_SH_jan=nanmean(SH_longterm_jan_seas-clim_ref_SH_jan);
%    diff_clim_SH_feb=nanmean(SH_longterm_feb_seas-clim_ref_SH_feb);
%    diff_clim_SH_mar=nanmean(SH_longterm_mar_seas-clim_ref_SH_mar);
%    diff_clim_SH_apr=nanmean(SH_longterm_apr_seas-clim_ref_SH_apr);
%    
%    SH_oct=nanmean(SH_longterm_oct_seas);
%    SH_nov=nanmean(SH_longterm_nov_seas);
%    SH_dec=nanmean(SH_longterm_dec_seas);
%    SH_jan=nanmean(SH_longterm_jan_seas);
%    SH_feb=nanmean(SH_longterm_feb_seas);
%    SH_mar=nanmean(SH_longterm_mar_seas);
%    SH_apr=nanmean(SH_longterm_apr_seas);


%    [diff_clim_SH_oct_10]=running_10(diff_clim_SH_oct); 
%    [diff_clim_SH_nov_10]=running_10(diff_clim_SH_nov); 
%    [diff_clim_SH_dec_10]=running_10(diff_clim_SH_dec); 
%    [diff_clim_SH_jan_10]=running_10(diff_clim_SH_jan); 
%    [diff_clim_SH_feb_10]=running_10(diff_clim_SH_feb); 
%    [diff_clim_SH_mar_10]=running_10(diff_clim_SH_mar); 
%    [diff_clim_SH_apr_10]=running_10(diff_clim_SH_apr); 
%    
%    %Trend 
% %jan_15=jan_15';
% [Tr_SH_oct_max, trend_SH_oct, xxx8, y8]=trend_rausch(diff_clim_SH_oct, jan_15);
% [Tr_SH_nov_max, trend_SH_nov, xxx9, y9]=trend_rausch(diff_clim_SH_nov, jan_15);
% [Tr_SH_dec_max, trend_SH_dec, xxx10, y10]=trend_rausch(diff_clim_SH_dec, jan_15);
% [Tr_SH_jan_max, trend_SH_jan, xxx11, y11]=trend_rausch(diff_clim_SH_jan, jan_15);
% [Tr_SH_feb_max, trend_SH_feb, xxx12, y12]=trend_rausch(diff_clim_SH_feb, jan_15);
% [Tr_SH_mar_max, trend_SH_mar, xxx13, y13]=trend_rausch(diff_clim_SH_mar, jan_15);
% [Tr_SH_apr_max, trend_SH_apr, xxx14, y14]=trend_rausch(diff_clim_SH_apr, jan_15);
   
      
 

% %for number of days with SH > 30cm per month
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
 
   
     
%Fliri-Probability Plot; Probability to have artificial snow conditions on a specific day. Call fliri_klima subfunction. Input: seasNr, number of
   %days per month, SH_longterm_seas :SH_count_oct_prob
   
%Shorter Period (1987/88-2006/07) last 20 seasons to detect changes in probability : SH_count_oct_prob_last20
   
[SH_count_oct_prob]=fliri_snow(seasNr, 31, SH_longterm_oct_seas); 
[SH_count_nov_prob]=fliri_snow(seasNr, 30, SH_longterm_nov_seas); 
[SH_count_dec_prob]=fliri_snow(seasNr, 31, SH_longterm_dec_seas); 
[SH_count_jan_prob]=fliri_snow(seasNr, 31, SH_longterm_jan_seas); 
[SH_count_feb_prob]=fliri_snow(seasNr, 28, SH_longterm_feb_seas); 
[SH_count_mar_prob]=fliri_snow(seasNr, 31, SH_longterm_mar_seas); 
[SH_count_apr_prob]=fliri_snow(seasNr, 30, SH_longterm_apr_seas); 

   
%   Merge monthly probabilities to yearly plot
  SH_prob_year(1:31,1:3)=SH_count_oct_prob;
  SH_prob_year(32:61,1:3)=SH_count_nov_prob;
  SH_prob_year(62:92,1:3)=SH_count_dec_prob;
  SH_prob_year(93:123,1:3)=SH_count_jan_prob;
  SH_prob_year(124:151,1:3)=SH_count_feb_prob;
  SH_prob_year(152:182,1:3)=SH_count_mar_prob;
  SH_prob_year(183:212,1:3)=SH_count_apr_prob; 
  
 date_year=date_winter(1:212); %dates for one season for plot (takes first season in time series but year is not shown in plot)

  
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
  
%   %WRITE OUTPUT FOR TABLES
%  
%Write files with number of snowdays per year and per year and month.
%Insert desired file name. 

jahr=U;
jahr(end)=[];
jahr=jahr';
data(:,1)=jahr;
data(:,2)=count_SH30; 
data(:,3)=season_mean_SH;

% filename1=['mehr als 30cm und mittel',location, '.txt'];
% dlmwrite(filename1,data,'delimiter', ' ', 'newline', 'pc');
% 
% data_all(:,1)=SH_count_oct;
% data_all(:,2)=SH_count_nov;
% data_all(:,3)=SH_count_dec;
% data_all(:,4)=SH_count_jan;
% data_all(:,5)=SH_count_feb;
% data_all(:,6)=SH_count_mar;
% data_all(:,7)=SH_count_apr;
% 
% 
% filename2=['mehr als 30cm',location, 'all','.txt'];
% dlmwrite(filename2,data_all,'delimiter', ' ' , 'newline', 'pc');


length(isnan(SH(i,:)))
%----------------------------------------------------------------------


%---------------------------section seasonal means over-------------------
%  
   m123= figure('name','Schneehöhe');
   
    subplot(2,1,1)
    set(gca,'fontsize',14)
    plot(jan_15,season_mean_SH,'r','linewidth',2)
    grid on
    axis tight
    ylabel('Schneehöhe (cm)')
    set(gca,'XTick',jan_15(1:10:end))
    set(gca,'XTickLabel',jan_15(1:10:end))
     datetick('x','keeplimits','keepticks')
    title(['Saisonale Schneehöhe, ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    bar(jan_15,count_SH30,'b')
    grid on
    axis tight
    ylabel('Anzahl an Tagen mit SH  \geq 30 cm')
    set(gca,'XTick',jan_15(1:10:end))
    set(gca,'XTickLabel',jan_15(1:10:end))
      datetick('x','keeplimits','keepticks')
    title(['Tage mit SH über 30 cm, ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
  %   %PLOTS
saveas(m123, 'fig3Snow.jpg');
saveas(m123, 'fig3Snow.fig');
    
   
% 
  dateplot=[date_daily(AA):date_daily(AAA)];    
% 
 m6=   figure('name','Probability Year');
    %[haxes,hline1,hline2]=plotyy(date_year,SH_prob_year(:,1),date_year,SH_prob_year(:,2),'area','plot');
    har1= area(date_year, SH_prob_year(:,1))
    set(har1, 'FaceColor', [1 105/255 180/255])
    hold on
    
    har2= area(date_year, tw_prob_year(:,1))
    set(har2, 'FaceColor', [0.2 0.2 0.2])
    child2=get(har2, 'Children');
    alpha(.5)
    grid off
    %axes(haxes(1))
    set(gca,'fontsize',14)
    ylabel('Wahrscheinlichkeit (%)')
    set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
    set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
    datetick('x','mmm','keeplimits','keepticks')
    xlabel('Monate')
    axis tight
   
    legend('Schneehöhe \geq 30 cm, ganze Zeitreihe', 'Tagesmittel Tf \leq -2°C, ganze Zeitreihe')
    title(['Wahrscheinlichkeit f. Schneehöhe \geq 30 cm, und Tagesmittel Tf \leq -2°C',' (', num2str(location),' ',num2str(alti), 'm)']);
    
saveas(m6, 'fig2Snow.jpg');
saveas(m6, 'fig2Snow.fig');
    


%CONTOUR Plot Schneehöhe
 a3=dateplot;
  size(SH_winter);
  mn=min(SH)
  mx=max(SH)
  map=zeros(mx, 3);
    
  map(1,1)=192/255;
  map(1,2)=192/255;
  map(1,3)=192/255;
  
  map(2:9, 1)=1;
  map(2:9, 2)=0;
  map(2:9, 3)=0;
  
  map(10:19, 1)=1;
  map(10:19, 2)=(140/255);
  map(10:19, 3)=0;
  
  map(20:29, 1)=1;
  map(20:29, 2)=(215/255);
  map(20:29, 3)=0;
  
  if mx<64
  map(30:mx, 1)=0;
  map(30:mx, 2)=1;
  map(30:mx, 3)=1;
  else
  map(30:64, 1)=0;
  map(30:64, 2)=1;
  map(30:64, 3)=1;

  map(65:100, 1)=0;
  map(65:100, 2)=191/255;
  map(65:100, 3)=1;
  
  if mx<199
      
    map(101:mx, 1)=123/255;
    map(101:mx, 2)=104/255;
    map(101:mx, 3)=238/255;
  
  else
    map(101:199, 1)=123/255;
    map(101:199, 2)=104/255;
    map(101:199, 3)=238/255;
  

        if mx<299
              map(200:mx, 1)=30/255;
              map(200:mx, 2)=144/255;
              map(200:mx, 3)=1;          
        else      
             map(200:299, 1)=30/255;
             map(200:299, 2)=144/255;
             map(200:299, 3)=1;
             
             if mx<399
                 
             map(300:mx, 1)=0;
             map(300:mx, 2)=0;
             map(300:mx, 3)=1;
             
             else
             map(300:399, 1)=0;
             map(300:399, 2)=0;
             map(300:399, 3)=1;
             
             if mx<499
            map(400:mx, 1)=25/255;
            map(400:mx, 2)=25/255;
            map(400:mx, 3)=112/255;
            
             else
            map(400:499, 1)=25/255;
            map(400:499, 2)=25/255;
            map(400:499, 3)=112/255;
             end
             end
        end
  end
  end
  
%   map(30:60, 1)=linspace(135/255,30/255,31);
%   map(30:60, 2)=linspace(206/255,144/255,31);
%   map(30:60, 3)=linspace(250/255,1,31);
  
  
size(a3)
size(SH_winter)
    %Contour Plot
     m1= figure('name','Entwicklung Tf');
     [C,h]=contourf(jan_15, a3, SH_winter,100);
     colormap(map)
     set(h,'LineStyle','none');
     set(gca,'XTick',jan_15(5:10:end))
     set(gca,'XTickLabel',jan_15(5:10:end))
     datetick('x','keeplimits','keepticks')
     datetick('y','keeplimits')     

     set(gca,'fontsize',14)
     axis tight
     xlabel('Jahre')
     grid on
     grid minor
     h=colorbar
     set(gca, 'CLim', [mn, mx])
     set(h, 'Xtick', [mn, mx])
     title(['Schneehöhe (cm) (', num2str(location),', ',num2str(alti),' m)']);
     
saveas(m1, 'fig1Snow.jpg');
saveas(m1, 'fig1Snow.fig');

