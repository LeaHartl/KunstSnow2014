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

location='VillAlpe';
alti=2140; %Seehöhe
 jahr_start=1948;
 jahr_end=2014;
 seasNr=jahr_end-jahr_start;
 %--------------------------------------------------!!!!!!!-------------       
     
%==========================================================================
% Daten einlesen
file1=fopen('Daten_ok\ZAMG\VillacherAlpeTageMittel1948_2014.txt');
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
pp=find(yyyy==1960 & mon==10 & dd==1);
pp2=find(yyyy==1991 & mon==4 & dd==30);

pp_a=find(unique(yyyy)==1960);
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


 

     
%Fliri-Probability Plot; Probability to have artificial snow conditions on a specific day. Call fliri_klima subfunction. Input: seasNr, number of
   %days per month, tw_longterm_seas :tw_count_oct_prob
   
%Shorter Period (1987/88-2006/07) last 20 seasons to detect changes in probability : tw_count_oct_prob_last20
   
[data_oct_prob]=fliri_abb_paper(seasNr, 31, tw_longterm_oct_seas);
[data_nov_prob]=fliri_abb_paper(seasNr, 30, tw_longterm_nov_seas); 
[data_dec_prob]=fliri_abb_paper(seasNr, 31, tw_longterm_dec_seas); 
[data_jan_prob]=fliri_abb_paper(seasNr, 31, tw_longterm_jan_seas); 
[data_feb_prob]=fliri_abb_paper(seasNr, 28, tw_longterm_feb_seas); 
[data_mar_prob]=fliri_abb_paper(seasNr, 31, tw_longterm_mar_seas); 
[data_apr_prob]=fliri_abb_paper(seasNr, 30, tw_longterm_apr_seas); 


%   Merge monthly probabilities to yearly plot
  tw_prob_year(1:31,1:3)=data_oct_prob;
  tw_prob_year(32:61,1:3)=data_nov_prob;
  tw_prob_year(62:92,1:3)=data_dec_prob;
  tw_prob_year(93:123,1:3)=data_jan_prob;
  tw_prob_year(124:151,1:3)=data_feb_prob;
  tw_prob_year(152:182,1:3)=data_mar_prob;
  tw_prob_year(183:212,1:3)=data_apr_prob; 
  
 date_year=date_winter(1:212); %dates for one season for plot (takes first season in time series but year is not shown in plot)

 abc=find(date_year==datenum(jahr_start, 12, 15, 0,0,0));
 dec15=date_year(abc);

 sum2=sum(tw_prob_year(1:abc,1));
 sum3=sum(tw_prob_year(1:abc,2));
 sum4=sum(tw_prob_year(1:abc,3));
 
 SI4=sum2-sum4;
 
data(:,1)=sum2;
data(:,2)=sum3;
data(:,3)=sum4;
data(:,4)=SI4;

filename1=['index2\SI_',location,'.txt'];
dlmwrite(filename1,data,'delimiter', ' ' , 'newline', 'pc');
 

  m7_1=  figure('name','Probability last 20 years');
    hline1=plot(date_year,tw_prob_year(:,1));
   hold on
    hline2=plot(date_year,tw_prob_year(:,2), 'k');
   hold on
    hline3=plot(date_year,tw_prob_year(:,3), 'r');
    %set(hline1,'FaceColor',[0.8 0.8 0.8])
%     set(hline2,'Color','b')
%    axes(haxes(1))
    set(gca,'fontsize',14)
    ylabel('Wahrscheinlichkeit (%)')
    set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
    set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
    datetick('x','mmm','keeplimits','keepticks')
    %    set(haxes(2), 'Box', 'Off');
    title(['Wahrscheinlichkeit f. Tagesmittelwert Feuchttemp. \leq -2°C, Mittelwert und Extremwerte (basierend auf 1993-2014), ', num2str(location),' ',num2str(alti)])
    axis tight
    hline=line([dec15 dec15],[0 100])
    set(hline, 'Color', 'k', 'LineWidth', 2);
    
saveas(m7_1, ['index2\', location,'.jpg']);
saveas(m7_1, ['index2\',location,'.fig']);
saveas(m7_1, ['index2\',location,'.eps']);
