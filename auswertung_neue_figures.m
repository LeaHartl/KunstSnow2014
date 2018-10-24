%Einlesen der MATLAB-Formatierten DWD Meteodaten aus files
clear all
close all

%***************
%EINSTELLUNGEN
%***************
%==========================================================================
%Tagesmitteldaten T und RH

%ohne februar 29.

%ohne februar 29.

location='Sonnblick';
alti=3109; %Seehöhe
 jahr_start=1948;
 jahr_end=2014;
 seasNr=jahr_end-jahr_start;
 %--------------------------------------------------!!!!!!!-------------       
     
%==========================================================================
% Daten einlesen
file1=fopen('Daten_ok\ZAMG\SonnblickTageMittel1948_2014.txt');
C=textscan(file1, '%f %f %f %f %f', 'headerlines',0);
fclose(file1);
%**************************************
%PROGRAM SECTION -> nicht verändern!
%**************************************

%calculate mean air pressure using simple barometric formula with altitude (for psychrometric formula) 
press=1013.25.*exp((alti.*(-1))./7290);

yyyy=C{1};
mon=C{2};
dd=C{3};
%careful here
rel_mean_d=C{4};
t_mean_d=C{5};


% t_mean_d_2030=t_mean_d+1;
% t_mean_d_2050=t_mean_d+1.8;

%function[tw_prob_year_last20, date_year] = klima_Tage(press, yyyy, mon, dd, rel_mean_d, t_mean_d); 

Q=find(mon==2 & dd==29);
yyyy(Q)=[];
mon(Q)=[];
dd(Q)=[];
rel_mean_d(Q)=[];
t_mean_d(Q)=[];
size(t_mean_d)
asdf=2






HH=zeros(size(yyyy));
mins=zeros(size(yyyy));
sec=zeros(size(yyyy));
datumvec=[yyyy, mon, dd,HH, mins, sec];
datestring=datestr(datumvec, 'yyyy-mm-dd');
date_winter=datenum(datestring, 'yyyy-mm-dd');
date_daily=date_winter;
  asdf=3
pp=find(yyyy==1961 & mon==10 & dd==1);
pp2=find(yyyy==1991 & mon==4 & dd==30);

pp_a=find(unique(yyyy)==1961);
pp2_a=find(unique(yyyy)==1990);

AA=find(yyyy==1971 & mon==10 & dd==1);
AAA=find(yyyy==1972 & mon==4 & dd==30);
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
    

asdf=1
%Calculate longterm wet-bulb temperatur - select Temperature mean and rh mean to use when calculating Tw
  [tw_mean_d]=wetbulb(t_mean_d, rel_mean_d, press);  
 
%if wetblub takes too long, write output to file once and then read tw from file. 
  data_1(:,1)=tw_mean_d;
 filename1 =['tw_tage.out'];
 %File schreiben
  dlmwrite(filename1,data_1,' ')
  
  
fixNaN1=find(isnan(t_mean_d)& isnan(rel_mean_d));
tw_mean_d(fixNaN1)=NaN;


%   file3=fopen(filename21);
%  E=textscan(file3, '%f', 'headerlines',0);
%  fclose(file3);
%  tw_mean_d=E{1};
%  
% % % Distribution of TW
% % 
%     figure('name','Distribution of Tw')
%     set(gca,'fontsize',12)
%     hist(tw_mean_d,-25:0.5:max(tw_mean_d));
%     grid on
%     axis tight   

%---------------------------Calculate seasonal means-------------------
%----------------------------------------------------------------------
k=212; %Number of days per winter season (Oct-Apr)
    
%Calculate seasonal mean t, rh and tw (1st of Oct-30th of April)
size(t_mean_d)
seasNr
t_winter = reshape(t_mean_d, k, seasNr);
 
season_mean_t=nanmean(t_winter,1); 
 

rh_winter = reshape(rel_mean_d, k, seasNr);
 
season_mean_rh=nanmean(rh_winter, 1);

tw_winter = reshape(tw_mean_d, k, seasNr);
 
season_mean_tw=nanmean(tw_winter, 1);





%set seasonal value to NaN if more than half the days are nan.
for i=1:seasNr;
    if length(  find(isnan(tw_winter(:,i))))   >106;
        season_mean_tw(i)=NaN;
    end
end
season_mean_tw    ;
%count number of days per season where tw is less than -2, -3, -4 and greater than -2.
    for i=1:seasNr;
     count_tw_2(i)=length(find(tw_winter(:,i)<=-2));
     count_tw_3(i)=length(find(tw_winter(:,i)<=-3));
     count_tw_4(i)=length(find(tw_winter(:,i)<=-4));
     count_tw_5(i)=length(find(tw_winter(:,i)>-2));
    end
        fixNaN=find(isnan(season_mean_tw));
    count_tw_2(fixNaN)=NaN;
    count_tw_3(fixNaN)=NaN;
    count_tw_4(fixNaN)=NaN;
    count_tw_5(fixNaN)=NaN;
  
%calculate standard deviation of count variables.    
    count_tw_2_std=std(count_tw_2);
    count_tw_3_std=std(count_tw_3);
    count_tw_4_std=std(count_tw_4);
    count_tw_5_std=std(count_tw_5);
    
%Calculate climatological reference period 60/61-90/91 out of seasonal means   
    clim_ref_t=nanmean(t_winter(pp:pp2));
    clim_ref_rh=nanmean(rh_winter(pp:pp2));
    clim_ref_tw=nanmean(tw_winter(pp:pp2));
    
    clim_ref_count_tw_2=round(nanmean(count_tw_2(pp_a:pp2_a)));
    clim_ref_count_tw_3=round(nanmean(count_tw_3(pp_a:pp2_a)));
    clim_ref_count_tw_4=round(nanmean(count_tw_4(pp_a:pp2_a)));
    clim_ref_count_tw_5=round(nanmean(count_tw_5(pp_a:pp2_a)));
     
%Abweichung der saisonalen Mittel vom klimatologischen Mittel berechnen
    diff_clim_t=season_mean_t-clim_ref_t;
    diff_clim_rh=season_mean_rh-clim_ref_rh;
    diff_clim_tw=season_mean_tw-clim_ref_tw;
    
    diff_clim_count_tw_2=count_tw_2-clim_ref_count_tw_2;
    diff_clim_count_tw_3=count_tw_3-clim_ref_count_tw_3;
    diff_clim_count_tw_4=count_tw_4-clim_ref_count_tw_4; 
    diff_clim_count_tw_5=count_tw_5-clim_ref_count_tw_5; 
     
%Calculate 10-season running-mean
 
 [diff_clim_t_10]=running_10(diff_clim_t); 
 [diff_clim_rh_10]=running_10(diff_clim_rh); 
 [diff_clim_tw_10]=running_10(diff_clim_tw); 
 
 [diff_clim_count_tw_2_10]=running_10(diff_clim_count_tw_2); 
 [diff_clim_count_tw_3_10]=running_10(diff_clim_count_tw_3); 
 [diff_clim_count_tw_4_10]=running_10(diff_clim_count_tw_4); 
 [diff_clim_count_tw_5_10]=running_10(diff_clim_count_tw_5); 
 
     
%Trend/Rauschverhältnis berechnen, linearen Trend berechnen
      
 %signifikanten Trend finden; Trend ist dann signifikant, wenn Trend/Rauschverhältnis >1.64 (90% Confidence-Level)
jan_15=jan_15' ;

[Tr_max2, trend2, xxxcount2, ycount2]=trend_rausch(count_tw_2, jan_15);

[Tr_max, trend, xxx, y]=trend_rausch(diff_clim_t, jan_15);
[Tr_rh_max, trend_rh, xxx2, y2]=trend_rausch(diff_clim_rh, jan_15);
[Tr_tw_max, trend_tw, xxx3, y3]=trend_rausch(diff_clim_tw, jan_15);
    Tr_tw_max
   % y3
    %mkhgvbgj
[Tr_count_tw_2_max, trend_count_tw_2, xxx4, y4]=trend_rausch(diff_clim_count_tw_2, jan_15);
[Tr_count_tw_3_max, trend_count_tw_3, xxx5, y5]=trend_rausch(diff_clim_count_tw_3, jan_15);
[Tr_count_tw_4_max, trend_count_tw_4, xxx6, y6]=trend_rausch(diff_clim_count_tw_4, jan_15);
[Tr_count_tw_5_max, trend_count_tw_5, xxx7, y7]=trend_rausch(diff_clim_count_tw_5, jan_15);

        
%CONTOUR Plot mit Entwicklung Anzahl Tage/Feuchtemperatur in Jahren
max(tw_mean_d);
tw_mean_d;
min(tw_mean_d);
zzz=round(abs(round(max(tw_mean_d)))+abs(round(min(tw_mean_d))))+1; 

tw_freq=zeros(seasNr,zzz);
   
for i=1:seasNr-1
    tw_freq(i,:)=hist(tw_winter(:,i),zzz);
end

step=abs(round(max(tw_mean_d))-round(min(tw_mean_d)))./(size(tw_freq,2)-1);
step=1;
a3=round(min(tw_mean_d)):step:round(max(tw_mean_d));
[XII,YII]=meshgrid(jan_15,a3);
tw_freq=rot90(tw_freq);
ZII=griddata(jan_15,a3,tw_freq,XII,YII);
ZII=flipud(ZII);
     
%---------------------------section seasonal means over-------------------
%-------------------------------------------------------------------------

%---------------------------Calulate monthly means - same calculations as for seasonal means but for single months -------------------

  [tw_longterm_oct_seas,tw_count_oct, tw_longterm_nov_seas,tw_count_nov, tw_longterm_dec_seas,tw_count_dec, tw_longterm_jan_seas,tw_count_jan, ... 
    tw_longterm_feb_seas,tw_count_feb, tw_longterm_mar_seas,tw_count_mar, tw_longterm_apr_seas,tw_count_apr, tab_tw_count,tw_longterm_oct, tw_longterm_nov,...
    tw_longterm_dec, tw_longterm_jan, tw_longterm_feb, tw_longterm_mar, tw_longterm_apr]=monthly_tw(mon, tw_mean_d, seasNr);

%--------------------------------------------------------------------
%   plots und berechnung verteilungsfkt

 %     distributionplots(location, yyyy, mon, dd, tw_mean_d, t_mean_d,rel_mean_d, press, tw_longterm_oct, tw_longterm_nov, tw_longterm_dec, tw_longterm_jan, tw_longterm_feb, tw_longterm_mar, tw_longterm_apr, U)
 % index1(location, yyyy, mon, dd, tw_mean_d, t_mean_d,rel_mean_d, press, tw_longterm_oct, tw_longterm_nov, tw_longterm_dec, tw_longterm_jan, tw_longterm_feb, tw_longterm_mar, tw_longterm_apr, U)
 
%     distributionplots_kurz(location, yyyy, mon, dd, tw_mean_d, t_mean_d,rel_mean_d, press, tw_longterm_oct, tw_longterm_nov, tw_longterm_dec, tw_longterm_jan, tw_longterm_feb, tw_longterm_mar, tw_longterm_apr, U)
 %----------------------------------------

   %Calculate climatological reference (1961-90)
   clim_ref_count_tw_oct=round(nanmean(tw_count_oct(pp_a:pp2_a)));
   clim_ref_count_tw_nov=round(nanmean(tw_count_nov(pp_a:pp2_a)));
   clim_ref_count_tw_dec=round(nanmean(tw_count_dec(pp_a:pp2_a)));
   clim_ref_count_tw_jan=round(nanmean(tw_count_jan(pp_a:pp2_a)));
   clim_ref_count_tw_feb=round(nanmean(tw_count_feb(pp_a:pp2_a)));
   clim_ref_count_tw_mar=round(nanmean(tw_count_mar(pp_a:pp2_a)));
   clim_ref_count_tw_apr=round(nanmean(tw_count_apr(pp_a:pp2_a)));
   
   %Abweichung vom langjährigen Mittel
   diff_clim_count_tw_oct=tw_count_oct-clim_ref_count_tw_oct;
   diff_clim_count_tw_nov=tw_count_nov-clim_ref_count_tw_nov;
   diff_clim_count_tw_dec=tw_count_dec-clim_ref_count_tw_dec;
   diff_clim_count_tw_jan=tw_count_jan-clim_ref_count_tw_jan;
   diff_clim_count_tw_feb=tw_count_feb-clim_ref_count_tw_feb;
   diff_clim_count_tw_mar=tw_count_mar-clim_ref_count_tw_mar;
   diff_clim_count_tw_apr=tw_count_apr-clim_ref_count_tw_apr;

   [diff_clim_count_tw_oct_10]=running_10(diff_clim_count_tw_oct); 
   [diff_clim_count_tw_nov_10]=running_10(diff_clim_count_tw_nov); 
   [diff_clim_count_tw_dec_10]=running_10(diff_clim_count_tw_dec); 
   [diff_clim_count_tw_jan_10]=running_10(diff_clim_count_tw_jan); 
   [diff_clim_count_tw_feb_10]=running_10(diff_clim_count_tw_feb); 
   [diff_clim_count_tw_mar_10]=running_10(diff_clim_count_tw_mar) ;
   [diff_clim_count_tw_apr_10]=running_10(diff_clim_count_tw_apr); 
 
%Trend 
jan_15=jan_15';
[Tr_count_tw_oct_max, trend_count_tw_oct, xxx8, y8]=trend_rausch(diff_clim_count_tw_oct, jan_15);
[Tr_count_tw_nov_max, trend_count_tw_nov, xxx9, y9]=trend_rausch(diff_clim_count_tw_nov, jan_15);
[Tr_count_tw_dec_max, trend_count_tw_dec, xxx10, y10]=trend_rausch(diff_clim_count_tw_dec, jan_15);

[Tr_count_tw_jan_max, trend_count_tw_jan, xxx11, y11]=trend_rausch(diff_clim_count_tw_jan, jan_15);
[Tr_count_tw_feb_max, trend_count_tw_feb, xxx12, y12]=trend_rausch(diff_clim_count_tw_feb, jan_15);
[Tr_count_tw_mar_max, trend_count_tw_mar, xxx13, y13]=trend_rausch(diff_clim_count_tw_mar, jan_15);
[Tr_count_tw_apr_max, trend_count_tw_apr, xxx14, y14]=trend_rausch(diff_clim_count_tw_apr, jan_15);
   
     
%Fliri-Probability Plot; Probability to have artificial snow conditions on a specific day. Call fliri_klima subfunction. Input: seasNr, number of
   %days per month, tw_longterm_seas :tw_count_oct_prob
   
%Shorter Period (1987/88-2006/07) last 20 seasons to detect changes in probability : tw_count_oct_prob_last20
   
% [tw_count_oct_prob, tw_count_oct_prob_last20, tw_count_oct_prob_last20_2]=fliri_klima(seasNr, 31, tw_longterm_oct_seas); 
% [tw_count_nov_prob, tw_count_nov_prob_last20, tw_count_nov_prob_last20_2]=fliri_klima(seasNr, 30, tw_longterm_nov_seas); 
% [tw_count_dec_prob, tw_count_dec_prob_last20, tw_count_dec_prob_last20_2]=fliri_klima(seasNr, 31, tw_longterm_dec_seas); 
% [tw_count_jan_prob, tw_count_jan_prob_last20, tw_count_jan_prob_last20_2]=fliri_klima(seasNr, 31, tw_longterm_jan_seas); 
% [tw_count_feb_prob, tw_count_feb_prob_last20, tw_count_feb_prob_last20_2]=fliri_klima(seasNr, 28, tw_longterm_feb_seas); 
% [tw_count_mar_prob, tw_count_mar_prob_last20, tw_count_mar_prob_last20_2]=fliri_klima(seasNr, 31, tw_longterm_mar_seas); 
% [tw_count_apr_prob, tw_count_apr_prob_last20, tw_count_apr_prob_last20_2]=fliri_klima(seasNr, 30, tw_longterm_apr_seas); 

   
% %   Merge monthly probabilities to yearly plot
%   tw_prob_year(1:31,1:4)=tw_count_oct_prob;
%   tw_prob_year(32:61,1:4)=tw_count_nov_prob;
%   tw_prob_year(62:92,1:4)=tw_count_dec_prob;
%   tw_prob_year(93:123,1:4)=tw_count_jan_prob;
%   tw_prob_year(124:151,1:4)=tw_count_feb_prob;
%   tw_prob_year(152:182,1:4)=tw_count_mar_prob;
%   tw_prob_year(183:212,1:4)=tw_count_apr_prob; 
%   
%  date_year=date_winter(1:212); %dates for one season for plot (takes first season in time series but year is not shown in plot)
% 
%   
%   %Merge monthly probabilities to yearly plot
%   tw_prob_year_last20(1:31,1:4)=tw_count_oct_prob_last20;
%   tw_prob_year_last20(32:61,1:4)=tw_count_nov_prob_last20;
%   tw_prob_year_last20(62:92,1:4)=tw_count_dec_prob_last20;
%   tw_prob_year_last20(93:123,1:4)=tw_count_jan_prob_last20;
%   tw_prob_year_last20(124:151,1:4)=tw_count_feb_prob_last20;
%   tw_prob_year_last20(152:182,1:4)=tw_count_mar_prob_last20;
%   tw_prob_year_last20(183:212,1:4)=tw_count_apr_prob_last20; 
%     
%   %Merge monthly probabilities to yearly plot
%   tw_prob_year_last20_2(1:31,1:4)=tw_count_oct_prob_last20_2;
%   tw_prob_year_last20_2(32:61,1:4)=tw_count_nov_prob_last20_2;
%   tw_prob_year_last20_2(62:92,1:4)=tw_count_dec_prob_last20_2;
%   tw_prob_year_last20_2(93:123,1:4)=tw_count_jan_prob_last20_2;
%   tw_prob_year_last20_2(124:151,1:4)=tw_count_feb_prob_last20_2;
%   tw_prob_year_last20_2(152:182,1:4)=tw_count_mar_prob_last20_2;
%   tw_prob_year_last20_2(183:212,1:4)=tw_count_apr_prob_last20_2; 
  
  

  %-----------------------------------------------------------------------------
%[central_year_sorted, window_with, window_with_sorted, sen_mk_sorted, sen_mk_sorted_2, sen_mk_sorted_3, sen_mk_sorted_4,...
 %   sen_mk_sorted_oct, sen_mk_sorted_nov, sen_mk_sorted_dec , sen_mk_sorted_jan, sen_mk_sorted_feb, sen_mk_sorted_mar, sen_mk_sorted_apr,  min_trend_abs,  max_trend_abs, sig_mk_4_sorted]...
%    =
trend_version3(location, alti,diff_clim_t, diff_clim_rh, diff_clim_tw,diff_clim_count_tw_2 , jan_15, diff_clim_count_tw_oct,diff_clim_count_tw_nov,...
    diff_clim_count_tw_dec,diff_clim_count_tw_jan,diff_clim_count_tw_feb,diff_clim_count_tw_mar,diff_clim_count_tw_apr );    %, window_with_sorted, sen_mk_sorted, sen_mk_sorted_2,...
    %sen_mk_sorted_3, sen_mk_sorted_4, sen_mk_sorted_oct, sen_mk_sorted_nov, sen_mk_sorted_dec , sen_mk_sorted_jan, sen_mk_sorted_feb, sen_mk_sorted_mar, sen_mk_sorted_apr]...
     
   %diff_clim_t
  
  
  %------------------------------------------------------------------------------
  
  
  
  
  
% %     last20mean=nanmean(tw_prob_year_last20(:,1))
% %   last20_2mean=nanmean(tw_prob_year_last20_2(:,1))
% %   
% %   filename22=['twprobyear',location,'.out'];
% %   filename222=['twprobyearlast20',location,'.out'];
% %   
% %    dlmwrite(filename22,tw_prob_year,'delimiter', ' ', 'newline', 'pc');
% %    dlmwrite(filename222,tw_prob_year_last20,'delimiter', ' ', 'newline', 'pc');
% %    
% % %   %WRITE OUTPUT FOR TABLES
% % % [maxT, ImaxT]=max(diff_clim_t);
% % % [miniT, IminT]=min(diff_clim_t);
% % % date_maxT=datestr(jan_15(ImaxT));
% % % date_minT=datestr(jan_15(IminT));
% % % 
% % % [maxRH, ImaxRH]=max(diff_clim_rh);
% % % [minRH, IminRH]=min(diff_clim_rh);
% % % date_maxRH=datestr(jan_15(ImaxRH));
% % % date_minRH=datestr(jan_15(IminRH));
% % % 
% % % [maxTW, ImaxTW]=max(diff_clim_tw);
% % % [minTW, IminTW]=min(diff_clim_tw);
% % % date_maxTW=datestr(jan_15(ImaxTW));
% % % date_minTW=datestr(jan_15(IminTW));
% % % 
% % % [maxtage, Itagemax]=max(count_tw_2);
% % % [mintage, Itagemin]=min(count_tw_2);
% % % datemaxTage=datestr(jan_15(Itagemax))
% % % dateminTage=datestr(jan_15(Itagemin))
% % %   
% % % 
% % % datestr(jan_15(xxx4),'yyyy');
% % % datestr(jan_15(xxx4+1),'yy');
% % % datestr(jan_15(end-1),'yyyy');
% % % datestr(jan_15(end));
% % % Tr_count_tw_2_max;
% %   
% % 
% % 
% % % %Write files with number of snowdays per year and per year and month.
% % % %Insert desired file name. 
% % jahr=U;
% % jahr(1)=[];
% % jahr=jahr';
% % data=[jahr; count_tw_2];   
% % filename1=['schneitage',location, '.txt'];
% % dlmwrite(filename1,data,'delimiter', ' ', 'newline', 'pc');
% % 
% % data_all(:,1)=tw_count_oct;
% % data_all(:,2)=tw_count_nov;
% % data_all(:,3)=tw_count_dec;
% % data_all(:,4)=tw_count_jan;
% % data_all(:,5)=tw_count_feb;
% % data_all(:,6)=tw_count_mar;
% % data_all(:,7)=tw_count_apr;
% %   
% % 
% % last20_2=nanmean(tw_prob_year_last20_2)
% % last20=nanmean(tw_prob_year_last20)
% % tr=trend
% % 
% % filename2=['schneitage',location, 'all','.txt'];
% % dlmwrite(filename2,data_all,'delimiter', ' ' , 'newline', 'pc');
% % 
% % meanvalue(:,1)=season_mean_t;
% % meanvalue(:,2)=season_mean_rh;
% % meanvalue(:,3)=season_mean_tw;
% % 
% % filename1=['means_',location,'.txt'];
% % dlmwrite(filename1,meanvalue,'delimiter', ' ' , 'newline', 'pc');
% % 
% % %  PLOTS
% % 
% % 
% % 
% %    
% % %   Trendanalyse Lufttemperatur und RH
% % 
%  %plot_trendT(season_mean_t, jan_15, jahr_start, jahr_end, location, alti, diff_clim_t, diff_clim_t_10, xxx, Tr_max, clim_ref_t, trend, y); 
%  %plot_trendRH(season_mean_rh, jan_15, jahr_start, jahr_end, location, alti, diff_clim_rh, diff_clim_rh_10, xxx2, Tr_rh_max, clim_ref_rh, trend_rh, y2); 
% % plot_trendTW(season_mean_tw, jan_15, jahr_start, jahr_end, location, alti, diff_clim_tw, diff_clim_tw_10, xxx3, Tr_tw_max, clim_ref_tw, trend_tw, y3, ...
% %     count_tw_2, count_tw_3, count_tw_4, count_tw_5,diff_clim_count_tw_2,diff_clim_count_tw_2_10, y4, xxx4, Tr_count_tw_2_max, trend_count_tw_2 );
% % 
% % % 
% % %          
% %     %Contour Plot
% %    m1= figure('name','Entwicklung Tf');
% %     [C,h]=contourf(jan_15,a3,ZII,50);
% %     colormap(jet)
% %     set(h,'LineStyle','none');
% %     set(gca,'Layer','bottom')
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     colorbar
% %     title('Entwicklung der Feuchttemperatur')
% %     set(gca,'fontsize',14)
% %     axis tight
% %     xlabel('Jahre')
% %     ylabel('Feuchttemperatur (°C)')
% %     hold on
% %     bar(jan_15,diff_clim_tw,'k')
% % saveas(m1, 'fig5.jpg');
% % saveas(m1, 'fig5.fig');
% % 
% %     
% % m2=    figure('name','Trend number of days monthly');
% %     subplot(3,1,1)
% %     set(gca,'fontsize',14)
% %     plot(jan_15,tw_count_oct,'r',jan_15,tw_count_nov,'-.k')
% %     grid off
% %     axis tight
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     legend('Oct','Nov','location','best')
% %     xlabel('Jahre')
% %     %ylabel('Anomalie d. Tagesanzahl mit Tagesmittel Tf \leq -2°C')
% %     title('Vorsaison (Oktober & November)')
% %     
% %     subplot(3,1,2)
% %     set(gca,'fontsize',14)
% %     plot(jan_15,tw_count_dec,'r',jan_15,tw_count_jan,'-.k',jan_15,tw_count_feb,':b')
% %     grid off
% %     axis tight
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     legend('Dez','Jan','Feb','location','best')
% %     xlabel('Jahre')
% %     ylabel('Tagesanzahl mit Tagesmittel Tf  \leq -2°C')
% %     title('Hauptsaison (Dezember - Februar)')
% %     
% %     subplot(3,1,3)
% %     set(gca,'fontsize',14)
% %     plot(jan_15,tw_count_mar,'r',jan_15,tw_count_apr,'-.k')
% %     grid off
% %     axis tight
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     legend('Mar','Apr','location','best')
% %     xlabel('Jahre')
% %    % ylabel('Tagesanzahl mit Tagesmittel Tf  \leq -2°C')
% %     title('Nachsaison (März - April)')
% % saveas(m2, 'fig6.jpg');
% % saveas(m2, 'fig6.fig');
% % 
% % 
% % 
% %     
% % m3=    figure('name','Anomalie d. Anzahl an Tagen nach Monat: Vorsaison');
% %     %Trend number of days
% %     
% %     subplot(2,1,1)
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     hold on
% %     bar(jan_15,diff_clim_count_tw_oct,'b')
% %     plot(jan_15,diff_clim_count_tw_oct_10,'k','linewidth',2)
% %     if Tr_count_tw_oct_max>1.64
% %     plot(jan_15,y8(:,xxx8),'-.r','linewidth',2)
% %     end
% %     axis tight
% %     grid off
% %     xlabel('Jahre')
% %     ylabel('Anomalie d. Anzahl an Tagen im Oktober (Tagesmittel Tf  \leq -2°C)')
% %     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_oct)) round(max(diff_clim_count_tw_oct))])
% %     xpos=get(gca,'xlim');
% %     xmittl=(xpos(1)+xpos(2))./2;
% %     ypos=get(gca,'ylim');
% %     ymittl=(ypos(1)+ypos(2))./2;
% %     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(nanstd(diff_clim_count_tw_oct),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
% %     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_tw_oct_max,-2)),' (',num2str(datestr(jan_15(xxx8),'yyyy')),'/',num2str(datestr(jan_15(xxx8+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     if Tr_count_tw_oct_max>1.64
% %     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_tw_oct,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     end
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     set(gca,'fontsize',14)
% %      title('Vorsaison (Oktober & November)')
% %     
% %     subplot(2,1,2)
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     hold on
% %     bar(jan_15,diff_clim_count_tw_nov,'b')
% %     plot(jan_15,diff_clim_count_tw_nov_10,'k','linewidth',2)
% %     if Tr_count_tw_nov_max>1.64
% %     plot(jan_15,y9(:,xxx9),'-.r','linewidth',2)
% %     end
% %     axis tight
% %     grid off
% %     xlabel('Jahre')
% %     ylabel('Anomalie d. Anzahl an Tagen im November (Tagesmittel Tf  \leq -2°C)')
% %     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_nov)) round(max(diff_clim_count_tw_nov))])
% %     xpos=get(gca,'xlim');
% %     xmittl=(xpos(1)+xpos(2))./2;
% %     ypos=get(gca,'ylim');
% %     ymittl=(ypos(1)+ypos(2))./2;
% %     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(nanstd(diff_clim_count_tw_nov),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
% %     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_tw_nov_max,-2)),' (',num2str(datestr(jan_15(xxx9),'yyyy')),'/',num2str(datestr(jan_15(xxx9+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     if Tr_count_tw_nov_max>1.64
% %     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_tw_nov,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     end
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     
% % saveas(m3, 'fig7.jpg');
% % saveas(m3, 'fig7.fig');
% %     
% %     
% %     
% %     
% %  m4=   figure('name','Anomalie d. Anzahl an Tagen nach Monat: Hauptsaison');
% %     %Trend number of days
% %     
% %     subplot(3,1,1)
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     hold on
% %     bar(jan_15,diff_clim_count_tw_dec,'b')
% %     plot(jan_15,diff_clim_count_tw_dec_10,'k','linewidth',2)
% %     if Tr_count_tw_dec_max>1.64
% %     plot(jan_15,y10(:,xxx10),'-.r','linewidth',2)
% %     end
% %     axis tight
% %     grid off
% %     xlabel('Jahre')
% %     ylabel('Dezember (Tagesmittel Tf \leq -2°C)')
% %     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_dec)) round(max(diff_clim_count_tw_dec))])
% %     xpos=get(gca,'xlim');
% %     xmittl=(xpos(1)+xpos(2))./2;
% %     ypos=get(gca,'ylim');
% %     ymittl=(ypos(1)+ypos(2))./2;
% %     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(nanstd(diff_clim_count_tw_dec),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
% %     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_tw_dec_max,-2)),' (',num2str(datestr(jan_15(xxx10),'yyyy')),'/',num2str(datestr(jan_15(xxx10+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     %if Tr_count_tw_dec_max>1.64
% %     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_tw_dec,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %    % end
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     title('Hauptsaison (Dezember, Januar, Februar)');
% %      %title('Hauptsaison (Dezember, Januar) keine Anomalie im Januar');
% %     
% %     subplot(3,1,2)
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     hold on
% %     bar(jan_15,diff_clim_count_tw_jan,'b')
% %     plot(jan_15,diff_clim_count_tw_jan_10,'k','linewidth',2)
% %    % if Tr_count_tw_jan_max>1.64
% %     plot(jan_15,y11(:,xxx11),'-.r','linewidth',2)
% %  %   end
% %     axis tight
% %     grid off
% %     xlabel('Jahre')
% %     ylabel('Januar (Tagesmittel Tf \leq -2°C)')
% %     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_jan)) round(max(diff_clim_count_tw_jan))])
% %     xpos=get(gca,'xlim');
% %     xmittl=(xpos(1)+xpos(2))./2;
% %     ypos=get(gca,'ylim');
% %     ymittl=(ypos(1)+ypos(2))./2;
% %     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(nanstd(diff_clim_count_tw_jan),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
% % %     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_tw_jan_max,-2)),' (',num2str(datestr(jan_15(xxx11),'yyyy')),'/',num2str(datestr(jan_15(xxx11+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     if Tr_count_tw_jan_max>1.64
% %     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_tw_jan,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     end
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     set(gca,'fontsize',14)
% %         
% %     subplot(3,1,3)
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     hold on
% %     bar(jan_15,diff_clim_count_tw_feb,'b')
% %     plot(jan_15,diff_clim_count_tw_feb_10,'k','linewidth',2)
% %     if Tr_count_tw_feb_max>1.64
% %     plot(jan_15,y12(:,xxx12),'-.r','linewidth',2)
% %     end
% %     axis tight
% %     grid off
% %     xlabel('Jahre')
% %     ylabel('Februar (Tagesmittel Tf \leq -2°C)')
% %     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_feb)) round(max(diff_clim_count_tw_feb))])
% %     xpos=get(gca,'xlim');
% %     xmittl=(xpos(1)+xpos(2))./2;
% %     ypos=get(gca,'ylim');
% %     ymittl=(ypos(1)+ypos(2))./2;
% %     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(nanstd(diff_clim_count_tw_feb),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
% % %     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_tw_feb_max,-2)),' (',num2str(datestr(jan_15(xxx12),'yyyy')),'/',num2str(datestr(jan_15(xxx12+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %    % if Tr_count_tw_feb_max>1.64
% %     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_tw_feb,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %   %  end
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     
% %     saveas(m4, 'fig8.jpg');
% % saveas(m4, 'fig8.fig');
% %     
% %     
% %     
% %     
% %   m5=  figure('name','Anomalie d. Anzahl an Tagen nach Monat: Nachsaison');
% %    
% %     subplot(2,1,1)
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     hold on
% %     bar(jan_15,diff_clim_count_tw_mar,'b')
% %     plot(jan_15,diff_clim_count_tw_mar_10,'k','linewidth',2)
% %     if Tr_count_tw_mar_max>1.64
% %     plot(jan_15,y13(:,xxx13),'-.r','linewidth',2)
% %     end
% %     axis tight
% %     grid off
% %     xlabel('Jahre')
% %     ylabel('Anomalie d. Anzahl an Tagen im März (Tagesmittel Tf \leq -2°C)')
% %     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_mar)) round(max(diff_clim_count_tw_mar))])
% %     xpos=get(gca,'xlim');
% %     xmittl=(xpos(1)+xpos(2))./2;
% %     ypos=get(gca,'ylim');
% %     ymittl=(ypos(1)+ypos(2))./2;
% %     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(nanstd(diff_clim_count_tw_mar),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
% %     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_tw_mar_max,-2)),' (',num2str(datestr(jan_15(xxx13),'yyyy')),'/',num2str(datestr(jan_15(xxx13+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     if Tr_count_tw_mar_max>1.64
% %     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_tw_mar,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     end
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     title('Nachsaison (März & April)')
% %     
% %     subplot(2,1,2)
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     hold on
% %     bar(jan_15,diff_clim_count_tw_apr,'b')
% %     plot(jan_15,diff_clim_count_tw_apr_10,'k','linewidth',2)
% %     if Tr_count_tw_apr_max>1.64
% %     plot(jan_15,y14(:,xxx14),'-.r','linewidth',2)
% %     end
% %     axis tight
% %     grid off
% %     xlabel('Jahre')
% %     ylabel('Anomalie d. Anzahl an Tagen im April (Tagesmittel Tf \leq -2°C)')
% %     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_apr)) round(max(diff_clim_count_tw_apr))])
% %     xpos=get(gca,'xlim');
% %     xmittl=(xpos(1)+xpos(2))./2;
% %     ypos=get(gca,'ylim');
% %     ymittl=(ypos(1)+ypos(2))./2;
% %     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(nanstd(diff_clim_count_tw_apr),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
% %     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_tw_apr_max,-2)),' (',num2str(datestr(jan_15(xxx14),'yyyy')),'/',num2str(datestr(jan_15(xxx14+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     if Tr_count_tw_apr_max>1.64
% %     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_tw_apr,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
% %     end
% %     set(gca,'XTick',jan_15(2:10:end))
% %     set(gca,'XTickLabel',jan_15(2:10:end))
% %     datetick('x','keeplimits','keepticks')
% %     set(gca,'fontsize',12,'fontweight','bold')
% %     
% %     
% %     saveas(m5, 'fig9.jpg');
% % saveas(m5, 'fig9.fig');
% 
% % % 
%  dateplot=[date_daily(AA):date_daily(AAA-1)];    
% 
% 
% %  m6=   figure('name','Probability Year');
% %     [haxes,hline1,hline2]=plotyy(date_year,tw_prob_year(:,1),date_year,tw_prob_year(:,2),'area','plot');
% %     grid off
% %     set(hline1,'FaceColor',[0.8 0.8 0.8])
% %     set(hline2,'Color','b')
% %     axes(haxes(1))
% %     set(gca,'fontsize',14)
% %     ylabel('Wahrscheinlichkeit (%)')
% %     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
% %     set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
% %     set(haxes(2), 'Box', 'Off');
% %     set(haxes(1),'Box','Off') 
% %     datetick('x','mmm','keeplimits','keepticks')
% %     xlabel('Monate')
% %     axis tight
% %     axes(haxes(2))
% %     set(gca,'fontsize',14,'YColor','k')
% %     ylabel('Extremwerte Feuchttemperatur (°C)')
% %     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
% %     set(gca,'YTick',[-25 -20 -15 -10 -5 0 5])
% %     datetick('x','mmm','keeplimits','keepticks')
% %     xlabel('Monate')
% %     hold on
% %     plot(date_year,tw_prob_year(:,3),'r')
% %     hold on
% %     plot(date_year,tw_prob_year(:,4),'k')
% %     title(['Wahrscheinlichkeit f. Tagesmittelwert Feuchttemp. \leq -2°C, Mittelwert und Extremwerte (basierend auf ',num2str(jahr_start),' - ',num2str(jahr_end),' ) ,',num2str(location),' ',num2str(alti),' m']);
% %     axis tight
% %     hline=refline(0,-2);
% %     set(hline, 'Color', 'k', 'LineWidth', 2);
% %     
% %     saveas(m6, 'fig10.jpg');
% % saveas(m6, 'fig10.fig');
% % % 
% % 
% % 
% % %     
% % %     
% %     
%   m7=  figure('name','Probability last 20 years');
%     [haxes,hline1,hline2]=plotyy(date_year,tw_prob_year_last20(:,1),date_year,tw_prob_year_last20(:,2),'area','plot');
%     grid off
%     set(hline1,'FaceColor',[0.8 0.8 0.8])
%     set(hline2,'Color','b')
%     axes(haxes(1))
%     set(gca,'fontsize',14)
%     ylabel('Wahrscheinlichkeit (%)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
%         set(haxes(2), 'Box', 'Off');
%     set(haxes(1),'Box','Off') 
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
%     axis tight
%     axes(haxes(2))
%     set(gca,'fontsize',14,'YColor','k')
%     ylabel('Extremwerte Feuchttemperatur (°C)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     set(gca,'YTick',[-25 -20 -15 -10 -5 0 5])
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
%     hold on
%     plot(date_year,tw_prob_year_last20(:,3),'r')
%         hold on
%     plot(date_year,tw_prob_year_last20(:,4),'k')
%     title(['Wahrscheinlichkeit f. Tagesmittelwert Feuchttemp. \leq -2°C, Mittelwert und Extremwerte (basierend auf 1993-2014), ', num2str(location),' ',num2str(alti)])
%     axis tight
%     hline=refline(0,-2);
%     set(hline, 'Color', 'k', 'LineWidth', 2);
%     
%         saveas(m7, 'fig11.jpg');
% saveas(m7, 'fig11.fig');

%     
%         saveas(m7, 'fig11.jpg');
% saveas(m7, 'fig11.fig');
%     
% %     
% %   m8=  figure('name','Comparison Probabilities last 20 years/20 years before');
% %     area(date_year,tw_prob_year_last20(:,1),'facecolor','g')
% %     hold on
% %     area(date_year,tw_prob_year_last20_2(:,1))
% %     grid on
% %     ylabel('Wahrscheinlichkeit')
% %     xlabel('Monate')
% %     axis tight
% %     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
% %     datetick('x','mmm','keeplimits','keepticks')
% %     jahr_20=[num2str(jahr_end-20),'/ ',num2str(jahr_end-1919),'-',num2str(jahr_end-1),'/ ',num2str(jahr_end-2000)];
% %     jahr_40=[num2str(jahr_end-40),'/ ',num2str(jahr_end-1939),'-',num2str(jahr_end-21),'/ ',num2str(jahr_end-1940)];
% %     %legend('1987/88-2006/07','1967/68-1986/87')
% %     legend(jahr_20,jahr_40)
% %     
% %         saveas(m8, 'fig12.jpg');
% % saveas(m8, 'fig12.fig');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % % % end
% % % 
% % % toc
% % % t=toc;
% % % 
% % %display(t)
% % 
% % 
% % 
