%Einlesen der MATLAB-Formatierten ZAMG Meteodaten aus files
clear all
close all
%tic%measure time matlab needs to calculate (initiate timer)

%***************
%EINSTELLUNGEN
%***************
%==========================================================================
%%JEDER TAG IM FILE MUSS 24 STUNDEN EINTRÄGE HABEN; KEINE DATENLÜCKEN!!


%FEBRUAR 29 WIRD NICHT BERÜCKSICHTIGT. IN EINLESEFILE LÖSCHEN.

% %Ort,Dateiordner, Dateinamen und Start& Endjahr eingeben

location='Loferer Alm';
alti=1620; %Seehöhe
 jahr_start=1994;
 jahr_end=2014;
 seasNr=jahr_end-jahr_start;

%==========================================================================

file1=fopen('Daten_ok\ZAMG\fix_Loferer Alm.txt');
C=textscan(file1, '%f %f %f %f %f %f', 'headerlines',0);
fclose(file1);

% %    ADJUST FOR START/END YEAR!!
%     label={'' '' '' '' '' '' '' '' '' 1960 ...
%      '' '' '' '' '' '' '' '' '' 1970 '' '' '' '' '' '' '' '' '' 1980 '' '' '' '' '' '' '' '' '' 1990 '' '' '' '' '' '' '' '' '' 2000 ...
%      '' '' '' '' '' '' '' '' '' 2010 '' '' '' ''};
% 
% % Wendelstein
%      label={ '' '' '' '' 1960 ...
%      '' '' '' '' '' '' '' '' '' 1970 '' '' '' '' '' '' '' '' '' 1980 '' '' '' '' '' '' '' '' '' 1990 '' '' '' '' '' '' '' '' '' 2000 ...
%      '' '' '' '' '' '' '' '' '' 2010 '' ''};
%1987: 
%      label={ '' '' '' '' '' 1990 '' '' '' '' '' '' '' '' '' 2000 ...
%         '' '' '' '' '' '' '' '' '' 2010 '' '' '' ''};

      label={ '' '' '' '' '' 2000 '' '' '' '' '' '' '' '' '' 2010 '' '' '' ''};
% label={ '' '' 2010 '' '' '' ''};

%short: 
%     label={ '' '' '' '' '' '' 2000 '' '' '' '' '' '' 2007  };
%**************************************
%PROGRAM SECTION -> nicht verändern!
%**************************************
seasNr=jahr_end-jahr_start;
%calculate mean air pressure using simple barometric formula with altitude (for psychrometric formula) 
press=1013.25.*exp((alti.*(-1))./7290);


yyyy=C{1};
mon=C{2};
dd=C{3};
hour=C{4};
rh=C{5};
tl2m=C{6};



% XYX=find(rh<-900);
% rh(XYX)=NaN;
% 
% XYY=find(tl2m<-900);
% tl2m(XYY)=NaN;


Q=find(mon==2 & dd==29);
yyyy(Q)=[];
mon(Q)=[];
dd(Q)=[];
hour(Q)=[];
rh(Q)=[];
tl2m(Q)=[];

minute=zeros(size(dd));
sec=zeros(size(dd));

%size(find(mon==3))
datumvec=[yyyy, mon, dd, hour, minute, sec];
datestring=datestr(datumvec, 'yyyy-mm-dd HH') ;
datum=datenum(datestring, 'yyyy-mm-dd HH');
  
    %find 1st of January of each year for xticklabeling
    j=1;
    for i=jahr_start:jahr_end-1
    jan_1(j)=datenum(i+1,1,1,0,0,0);
    j=j+1;
    end
   
     %Calculate seasonal mean temperature (1st of Oct-30th of April), calls
     %subfunction seas_mean.

seasNr
size(tl2m)
[season_mean]=seas_mean(seasNr, tl2m);

    
    %Create date for seasonal mean temps (=15th January)
    %find 15 January of each year 
    j=1;
    for i=jahr_start:jahr_end-1;
    jan_15(j)=datenum(i+1,1,15,0,0,0);
    j=j+1;
    end
    
% %  %------------------------------------------------------------------------------------------------------- 
%         %avoid calling wetbulb sub funct to save time during trials.
% %     
%  %   Calculate Wet-Bulb Temperature (calls subfunction wetbulb. Takes a long time!!)    
% %    %  select Temperature mean and rh mean to use when calculating Tw
%    rh2=rh;
%    [tw2]=wetbulb(tl2m, rh2, press);
%    
%  data_2(:,1)=tw2;

% Filenamen 

 filename22 =['tw', location, '.out'];
% 
% % File schreiben
%     dlmwrite(filename22,data_2,' ');

% % % % % 
 file3=fopen(filename22);
 E=textscan(file3, '%f', 'headerlines',0);
 fclose(file3);
 tw2=E{1};
% % %---------------------------------------------------------------------------------------------------

     
%Calculate seasonal mean Tw (1st of Oct-30th of April)
 
[season_mean_tw]=seas_mean(seasNr, tw2);
    

%Call subfunction prep_tw_plot to prepare color coded Tw plot
%(tw_snow/nosnow for limitis -1.5, -2, -3, -4)&
%Calculate total snowmaking hours (snow_hours)&
%Calculate seasonal snowmaking hours for Tw limit=-1.5°C, -2°C, -3°C and
%-4°C (snow_hours_seas, array wit values for the four limits.)
       
[tw_snow, tw_nosnow, tw_snow_2, tw_nosnow_2, tw_snow_3, tw_nosnow_3,tw_snow_4, tw_nosnow_4, snow_hours, snow_hours_seas]=prep_tw_plot(tw2, seasNr, tl2m);

%Calculate standard deviation to detect year-to-year variability
stwd=std(snow_hours_seas,0,2);
        
    %Time distribution of snowmaking hours
    
    %Considering the entire period
     %hours     
     j=1;
    for i=0:23
    tt=find(hour==i);
    xx=tw_snow(tt);
    yy=isnan(xx);
    yy(yy==1)=[];
    zz1(j)=length(yy);
    j=j+1;
    
   clear tt xx yy
    end
  
   
    clear j
    zzz=sum(zz1);
    zz_rel=(zz1./zzz).*100;
    
    zz_day(1)=sum(zz1(10:18));
    zz_day(2)=sum(zz1(1:9))+sum(zz1(19:24));
    
      %months
    j=1;
    for i=10:12
    tt=find(mon==i);
    xx=tw_snow(tt);
    yy=isnan(xx);
    yy(yy==1)=[];
    zz_2(j)=length(yy);
    j=j+1;
    clear tt xx yy
    end
    for i=1:4
    tt=find(mon==i);
    xx=tw_snow(tt);
    yy=isnan(xx);
    yy(yy==1)=[];
    zz_2(j)=length(yy);
    j=j+1;
    clear tt xx yy
    end
   
 %FLIRI PLOT  with hourly values
   %find months in date variable
%    ss=datevec(datum);
%    ss2=ss(:,2);
   shortterm_oct=find(mon==10);
   shortterm_nov=find(mon==11);
   shortterm_dec=find(mon==12);
   shortterm_jan=find(mon==1);
   shortterm_feb=find(mon==2);
   shortterm_mar=find(mon==3);
   shortterm_apr=find(mon==4);
   
   %Wet-bulb Temp per months
   %[m,n]=size(snow_hours_seas);
  
   tw_shortterm_oct=tw2(shortterm_oct);
   tw_shortterm_nov=tw2(shortterm_nov);
   tw_shortterm_dec=tw2(shortterm_dec);
   tw_shortterm_jan=tw2(shortterm_jan);
   tw_shortterm_feb=tw2(shortterm_feb);
   tw_shortterm_mar=tw2(shortterm_mar);
   tw_shortterm_apr=tw2(shortterm_apr);
   
        clear kk i ii j 
 %call subfunction fliri_sub for each month of winter season. pass number
 %of hours and days per month to subfunc. (720 hours for months with 30 days, 744 for 31
 %days)
 n=seasNr;
 size(tw_shortterm_oct)
%Oct
[tw_shortterm_oct_hours, tw_count_oct_sum, std_oct_sum, mean_oct_sum]=fliri_sub(tw_shortterm_oct, 744, 31, n);
   
%Nov
[tw_shortterm_nov_hours, tw_count_nov_sum, std_nov_sum, mean_nov_sum]=fliri_sub(tw_shortterm_nov, 720, 30, n);

%Dez
[tw_shortterm_dec_hours, tw_count_dec_sum, std_dec_sum, mean_dec_sum]=fliri_sub(tw_shortterm_dec, 744, 31, n);

%Jan
[tw_shortterm_jan_hours, tw_count_jan_sum, std_jan_sum, mean_jan_sum]=fliri_sub(tw_shortterm_jan, 744, 31, n);

%Feb
[tw_shortterm_feb_hours, tw_count_feb_sum, std_feb_sum, mean_feb_sum]=fliri_sub(tw_shortterm_feb, 672, 28, n);

%Mar
[tw_shortterm_mar_hours, tw_count_mar_sum, std_mar_sum, mean_mar_sum]=fliri_sub(tw_shortterm_mar, 744, 31, n);    

%Apr
[tw_shortterm_apr_hours, tw_count_apr_sum, std_apr_sum, mean_apr_sum]=fliri_sub(tw_shortterm_apr, 720, 30, n);   


  tw_shortterm_year(1:31,1:3)=tw_shortterm_oct_hours;
  tw_shortterm_year(32:61,1:3)=tw_shortterm_nov_hours;
  tw_shortterm_year(62:92,1:3)=tw_shortterm_dec_hours;
  tw_shortterm_year(93:123,1:3)=tw_shortterm_jan_hours;
  tw_shortterm_year(124:151,1:3)=tw_shortterm_feb_hours;
  tw_shortterm_year(152:182,1:3)=tw_shortterm_mar_hours;
  tw_shortterm_year(183:212,1:3)=tw_shortterm_apr_hours;

%-------------------------------------------------------------------------------------------------------------------------   
  %$Tabelle Auswertung Beschneistunden monatlich
  zz2=size(tw_count_apr_sum,2);
  tab_tw_count(1,1:zz2)=tw_count_oct_sum;
  tab_tw_count(1,zz2+1)=std(tw_count_oct_sum);
  tab_tw_count(1,zz2+2)=mean(tw_count_oct_sum);
  tab_tw_count(2,1:zz2)=tw_count_nov_sum;
  tab_tw_count(2,zz2+1)=std(tw_count_nov_sum);
  tab_tw_count(2,zz2+2)=mean(tw_count_nov_sum);
  tab_tw_count(3,1:zz2)=tw_count_dec_sum;
  tab_tw_count(3,zz2+1)=std(tw_count_dec_sum);
  tab_tw_count(3,zz2+2)=mean(tw_count_dec_sum);
  tab_tw_count(4,1:zz2)=tw_count_jan_sum;
  tab_tw_count(4,zz2+1)=std(tw_count_jan_sum);
  tab_tw_count(4,zz2+2)=mean(tw_count_jan_sum);
  tab_tw_count(5,1:zz2)=tw_count_feb_sum;
  tab_tw_count(5,zz2+1)=std(tw_count_feb_sum);
  tab_tw_count(5,zz2+2)=mean(tw_count_feb_sum);
  tab_tw_count(6,1:zz2)=tw_count_mar_sum;
  tab_tw_count(6,zz2+1)=std(tw_count_mar_sum);
  tab_tw_count(6,zz2+2)=mean(tw_count_mar_sum);
  tab_tw_count(7,1:zz2)=tw_count_apr_sum;
  tab_tw_count(7,zz2+1)=std(tw_count_apr_sum);
  tab_tw_count(7,zz2+2)=mean(tw_count_apr_sum);
   
    % Propeller gun
    
    %Calculate resulting potential snow production in (m³/h) for propeller guns with
    %mean relationship of snow gun producers
    snow_mass_hour=-4.83.*tw_snow+3.94;
    
    %set upper limit for snow production to 72 m³/h (due to max. water
    %pressure and nb. of cones), lower limit is already considered in
    %variable 'tw_snow' (NaN's) 
    snow_mass_hour(snow_mass_hour>72)=72;
    
     %Calculate seasonal snowmaking potential (m³/season) for Tw
     %limit=-2°C, call subfunct. seas_mean    
    hoursperseason=212*24;

season_mean1=reshape(snow_mass_hour, hoursperseason, seasNr);

snow_mass_seas=nansum(season_mean1);
stwd_2=std(snow_mass_seas,0,2)


    
    %Lanze
    
    snow_mass_hour_lan=-3.94.*tw_snow-4.23;
    %set upper limit for snow production to 51 m³/h (due to max. water
    %pressure and nb. of cones), lower limit is already considered in
    %variable 'tw_snow' (NaN's) 
    snow_mass_hour_lan(snow_mass_hour_lan>51)=51;
    
     %Calculate seasonal snowmaking potential (m³/season) for Tw limit=-2°C
  season_mean11=reshape(snow_mass_hour_lan, hoursperseason, seasNr);

snow_mass_seas_lan=nansum(season_mean11);
stwd_3=std(snow_mass_seas_lan,0,2)
% 
% %Trends in Schneleistung---------------------------------------
% %Klimat. Referenz
% climref_snowmassProp=nansum(snow_mass_hour(pp:pp2))/30;
% climref_snowmassLan=nansum(snow_mass_hour_lan(pp:pp2))/30;
%     
% %Abweichung von der Referenz
% diff_clim_prop=snow_mass_seas-climref_snowmassProp;
% diff_clim_lan=snow_mass_seas_lan-climref_snowmassLan;
% 
% %Trend Rausch Verhältnis
% [Tr_maxProp, trendProp, xxxProp, yProp]=trend_rausch(diff_clim_prop, jan_15);
% [Tr_maxLan, trendLan, xxxLan, yLan]=trend_rausch(diff_clim_lan, jan_15);
% 
% %-------------------------------------------------------------

%     %Compare daily means to hourly means
%     
%     test23=find(tw2_daily<=-2);
%     test23_2=date_daily(test23);
%     
%     for i=1:size(test23_2,1)
%     cc(i)=find(datum==test23_2(i));
%     end
%     cc=cc';
%   
%     for i=1:size(cc,1);
%     tw_dist(i,:)=tw2(cc(i):cc(i)+23);
%     end
%     tw_dist=rot90(tw_dist);
%     tw_dist=flipud(tw_dist);
%     [m,n]=size(tw_dist);
%     tw_dist_2=reshape(tw_dist,m*n,1);
%     clear cc
%     
%     test24=find(tw2_daily<=-2&tw2_daily>-3);
%     test24_2=date_daily(test24);
%     
%     for i=1:size(test24_2,1)
%     cc(i)=find(datum==test24_2(i));
%     end
%     cc=cc';
%     
%     
%     for i=1:size(cc,1)
%     tw_dist_3(i,:)=tw2(cc(i):cc(i)+23); 
%     end
%     tw_dist_3=rot90(tw_dist_3);
%     tw_dist_3=flipud(tw_dist_3);
%     [m,n]=size(tw_dist_3);
%     tw_dist_4=reshape(tw_dist_3,m*n,1);
    

%Berechnung mittlerer Tagesgang
tw22=tw2;
tw2_dailyamp=reshape(tw22,24,length(tw22)./24);
tw2_dailyamp_mean=nanmean(tw2_dailyamp,2);

%Unterscheidung Beschneiung ja/nein
tw_snow_3=tw2_dailyamp_mean;
tw_nosnow_3=tw2_dailyamp_mean;

tw_snow_3(tw_snow_3>=-2)=NaN;
tw_nosnow_3(tw_nosnow_3<-2)=NaN;

%Umrechnung in Schneileistung
tw2_dailyamp_prop=-4.83.*tw2_dailyamp_mean+3.94;
tw2_dailyamp_lan=-3.94.*tw2_dailyamp_mean+4.23;
  tw2_dailyamp_lan(tw2_dailyamp_lan<12.11)=0;
  tw2_dailyamp_prop(tw2_dailyamp_prop<13.6)=0;
%  tw2_dailyamp_lan(tw2_dailyamp_lan<0)=0;
%  tw2_dailyamp_prop(tw2_dailyamp_prop<0)=0;
tw2_dailyamp_lan(tw2_dailyamp_lan>51)=51;
tw2_dailyamp_prop(tw2_dailyamp_prop>72)=72;

meanProp=nanmean(tw2_dailyamp_prop)
meanLan=nanmean(tw2_dailyamp_lan)

%Prozentueller Anteil Schneileistung in 6 Stunden Teilen
tw2_dailyamp_rel(1)=roundn((sum(tw2_dailyamp_prop(1:6))./sum(tw2_dailyamp_prop)).*100,0);
tw2_dailyamp_rel(2)=roundn((sum(tw2_dailyamp_prop(7:12))./sum(tw2_dailyamp_prop)).*100,0);
tw2_dailyamp_rel(3)=roundn((sum(tw2_dailyamp_prop(13:18))./sum(tw2_dailyamp_prop)).*100,0);
tw2_dailyamp_rel(4)=roundn((sum(tw2_dailyamp_prop(19:24))./sum(tw2_dailyamp_prop)).*100,0);

%Dasselbe monatsweise


[m,n]=size(snow_hours_seas);
% Q=reshape(tw_shortterm_oct, length(tw_shortterm_oct)/2, n);
% tw2_dailyamp_oct=nanmean(Q,2);
% %tw2_dailyamp_oct=nanmean(tw_shortterm_oct,2);  %this line replaces the
% %three above lines in orig. code. 


%oct
[tw2_dailyamp_oct_2, tw2_dailyamp_prop_oct, tw2_dailyamp_lan_oct, tw2_dailyamp_rel_oct, snow_mass_prop_oct, snow_mass_lan_oct, snow_mass_prop_oct_sum, snow_mass_prop_oct_sum_std, snow_mass_prop_oct_sum_mean, snow_mass_lan_oct_sum, snow_mass_lan_oct_sum_std, snow_mass_lan_oct_sum_mean]=tagesgang(tw_shortterm_oct, 31, n);

snow_mass_prop_oct=reshape(snow_mass_prop_oct, 24*31, seasNr);

snow_mass_prop_oct1=nansum(snow_mass_prop_oct);

snow_mass_lan_oct=reshape(snow_mass_lan_oct, 24*31, seasNr);

snow_mass_lan_oct1=nansum(snow_mass_lan_oct);

%nov
[tw2_dailyamp_nov_2, tw2_dailyamp_prop_nov, tw2_dailyamp_lan_nov, tw2_dailyamp_rel_nov, snow_mass_prop_nov, snow_mass_lan_nov, snow_mass_prop_nov_sum, snow_mass_prop_nov_sum_std, snow_mass_prop_nov_sum_mean, snow_mass_lan_nov_sum, snow_mass_lan_nov_sum_std, snow_mass_lan_nov_sum_mean]=tagesgang(tw_shortterm_nov, 30, n);

snow_mass_prop_nov=reshape(snow_mass_prop_nov, 24*30, seasNr);

snow_mass_prop_nov1=nansum(snow_mass_prop_nov);

snow_mass_lan_nov=reshape(snow_mass_lan_nov, 24*30, seasNr);

snow_mass_lan_nov1=nansum(snow_mass_lan_nov);

%dec
[tw2_dailyamp_dec_2, tw2_dailyamp_prop_dec, tw2_dailyamp_lan_dec, tw2_dailyamp_rel_dec, snow_mass_prop_dec, snow_mass_lan_dec, snow_mass_prop_dec_sum, snow_mass_prop_dec_sum_std, snow_mass_prop_dec_sum_mean, snow_mass_lan_dec_sum, snow_mass_lan_dec_sum_std, snow_mass_lan_dec_sum_mean]=tagesgang(tw_shortterm_dec, 31, n);

snow_mass_prop_dec=reshape(snow_mass_prop_dec, 24*31, seasNr);

snow_mass_prop_dec1=nansum(snow_mass_prop_dec);

snow_mass_lan_dec=reshape(snow_mass_lan_dec, 24*31, seasNr);

snow_mass_lan_dec1=nansum(snow_mass_lan_dec);

%jan
[tw2_dailyamp_jan_2, tw2_dailyamp_prop_jan, tw2_dailyamp_lan_jan, tw2_dailyamp_rel_jan, snow_mass_prop_jan, snow_mass_lan_jan, snow_mass_prop_jan_sum, snow_mass_prop_jan_sum_std, snow_mass_prop_jan_sum_mean, snow_mass_lan_jan_sum, snow_mass_lan_jan_sum_std, snow_mass_lan_jan_sum_mean]=tagesgang(tw_shortterm_jan, 31, n);

snow_mass_prop_jan=reshape(snow_mass_prop_jan, 24*31, seasNr);

snow_mass_prop_jan1=nansum(snow_mass_prop_jan);

snow_mass_lan_jan=reshape(snow_mass_lan_jan, 24*31, seasNr);

snow_mass_lan_jan1=nansum(snow_mass_lan_jan);

%feb
[tw2_dailyamp_feb_2, tw2_dailyamp_prop_feb, tw2_dailyamp_lan_feb, tw2_dailyamp_rel_feb, snow_mass_prop_feb, snow_mass_lan_feb, snow_mass_prop_feb_sum, snow_mass_prop_feb_sum_std, snow_mass_prop_feb_sum_mean, snow_mass_lan_feb_sum, snow_mass_lan_feb_sum_std, snow_mass_lan_feb_sum_mean]=tagesgang(tw_shortterm_feb, 28, n);

snow_mass_prop_feb=reshape(snow_mass_prop_feb, 24*28, seasNr);

snow_mass_prop_feb1=nansum(snow_mass_prop_feb);

snow_mass_lan_feb=reshape(snow_mass_lan_feb, 24*28, seasNr);

snow_mass_lan_feb1=nansum(snow_mass_lan_feb);

%mar
[tw2_dailyamp_mar_2, tw2_dailyamp_prop_mar, tw2_dailyamp_lan_mar, tw2_dailyamp_rel_mar, snow_mass_prop_mar, snow_mass_lan_mar, snow_mass_prop_mar_sum, snow_mass_prop_mar_sum_std, snow_mass_prop_mar_sum_mean, snow_mass_lan_mar_sum, snow_mass_lan_mar_sum_std, snow_mass_lan_mar_sum_mean]=tagesgang(tw_shortterm_mar, 31, n);

snow_mass_prop_mar=reshape(snow_mass_prop_mar, 24*31, seasNr);

snow_mass_prop_mar1=nansum(snow_mass_prop_mar);

snow_mass_lan_mar=reshape(snow_mass_lan_mar, 24*31, seasNr);

snow_mass_lan_mar1=nansum(snow_mass_lan_mar);

%apr
[tw2_dailyamp_apr_2, tw2_dailyamp_prop_apr, tw2_dailyamp_lan_apr, tw2_dailyamp_rel_apr, snow_mass_prop_apr, snow_mass_lan_apr, snow_mass_prop_apr_sum, snow_mass_prop_apr_sum_std, snow_mass_prop_apr_sum_mean, snow_mass_lan_apr_sum, snow_mass_lan_apr_sum_std, snow_mass_lan_apr_sum_mean]=tagesgang(tw_shortterm_apr, 30, n);

snow_mass_prop_apr=reshape(snow_mass_prop_apr, 24*30, seasNr);

snow_mass_prop_apr1=nansum(snow_mass_prop_apr); 

snow_mass_lan_apr=reshape(snow_mass_lan_apr, 24*30, seasNr);

snow_mass_lan_apr1=nansum(snow_mass_lan_apr);          

teste(1,:)=tw2_dailyamp_rel_oct;
teste(2,:)=tw2_dailyamp_rel_nov;
teste(3,:)=tw2_dailyamp_rel_dec;
teste(4,:)=tw2_dailyamp_rel_jan;
teste(5,:)=tw2_dailyamp_rel_feb;
teste(6,:)=tw2_dailyamp_rel_mar;
teste(7,:)=tw2_dailyamp_rel_apr;
teste(isnan(teste)) = 0 ;

teste2(1,:)=tw2_dailyamp_prop_oct;
teste2(2,:)=tw2_dailyamp_prop_nov;
teste2(3,:)=tw2_dailyamp_prop_dec;
teste2(4,:)=tw2_dailyamp_prop_jan;
teste2(5,:)=tw2_dailyamp_prop_feb;
teste2(6,:)=tw2_dailyamp_prop_mar;
teste2(7,:)=tw2_dailyamp_prop_apr;

teste3(1,:)=tw2_dailyamp_lan_oct;
teste3(2,:)=tw2_dailyamp_lan_nov;
teste3(3,:)=tw2_dailyamp_lan_dec;
teste3(4,:)=tw2_dailyamp_lan_jan;
teste3(5,:)=tw2_dailyamp_lan_feb;
teste3(6,:)=tw2_dailyamp_lan_mar;
teste3(7,:)=tw2_dailyamp_lan_apr;

% meanvalue(1,:)=season_mean;
% meanvalue(2,:)=season_mean_rh;
% meanvalue(3,:)=season_mean_tw;


tab_snowmass_prop(:,1)=snow_mass_prop_oct1;
tab_snowmass_prop(:,2)=snow_mass_prop_nov1;
tab_snowmass_prop(:,3)=snow_mass_prop_dec1;
tab_snowmass_prop(:,4)=snow_mass_prop_jan1;
tab_snowmass_prop(:,5)=snow_mass_prop_feb1;
tab_snowmass_prop(:,6)=snow_mass_prop_mar1;
tab_snowmass_prop(:,7)=snow_mass_prop_apr1;


tab_snowmass_lan(:,1)=snow_mass_lan_oct1;
tab_snowmass_lan(:,2)=snow_mass_lan_nov1;
tab_snowmass_lan(:,3)=snow_mass_lan_dec1;
tab_snowmass_lan(:,4)=snow_mass_lan_jan1;
tab_snowmass_lan(:,5)=snow_mass_lan_feb1;
tab_snowmass_lan(:,6)=snow_mass_lan_mar1;
tab_snowmass_lan(:,7)=snow_mass_lan_apr1;


% 
% filename2=['schneistundenMonat_',location,'.txt'];
% dlmwrite(filename2,tab_tw_count,'delimiter', ' ' , 'newline', 'pc');
% 
% filename3=['propeller_',location,'.txt'];
% dlmwrite(filename3,tab_snowmass_prop,'delimiter', ' ' , 'newline', 'pc');
% 
% filename4=['lanze_',location,'.txt'];
% dlmwrite(filename4,tab_snowmass_lan,'delimiter', ' ' , 'newline', 'pc');
% 
% filename5=['snohours_seas_',location,'.txt'];
% dlmwrite(filename5,snow_hours_seas,'delimiter', ' ' , 'newline', 'pc');
% % %   PLOTS
% % subfunctions plot1, 2, 3 
% 
     
% % Distribution of TW
% 
%     figure('name','Distribution of Tw')
%     set(gca,'fontsize',12)
%     hist(tw2,-25:0.5:10);
%     grid on
%     axis tight   
% 
% 
% Only plot winter data:

% tl2m(summer)=NaN;
% rh(summer)=NaN;

tw_nosnow1=tl2m;
tw_nosnow1=tw_nosnow;

tw_snow1=tl2m;
tw_snow1=tw_snow;
lala=isnan(tw_snow1);
% 
%-----------------------------------------plot1-----------------------------------------------------------
 %subfunction plot1 makes three part figure of Temp, RH and wetbulb temp (snow and now
   %snow.
plot1(datum, tl2m, rh, tw_nosnow1, tw_snow1, jan_1, label);
% -----------------------------------------plot1-----------------------------------------------------------
% 
% 
% 
% -----------------------------------------plot2-----------------------------------------------------------
%  subfunction plot2 plots two part figure of seasonal snowmaking hours
%  (total and with different TW limits)
plot2(jan_15, snow_hours_seas, jan_1, location, alti, label);
% -----------------------------------------plot2-----------------------------------------------------------
% 
% -----------------------------------------plot3-----------------------------------------------------------
%  subfunction plot33 plots four part figure of temporal distribution of
%  snowmaking hours (2 bar graphs, 2 pie charts.

plot33(zz1, location, zz_day, zz_2);
%-----------------------------------------plot3-----------------------------------------------------------

%----------------------------------------------------plot4-----------------------------------------------------------
%Plot4 subfunction plots snow production potential in (m³/h) for propeller
%guns and lance gun, hourly and seasonal sum.


 plot4(location, datum, jan_1, jan_15, snow_mass_hour, snow_mass_seas, snow_mass_hour_lan, snow_mass_seas_lan, label)   ;
%----------------------------------------------------plot4-----------------------------------------------------------


%----------------------------------------------------plot5-----------------------------------------------------------
%subfunction plot5 plottet kumulative verteilungsfunktion
%plot5(location, tw2, tw_dist_2, tw_dist_4 );
%----------------------------------------------------plot5-----------------------------------------------------------


%-----------------------------------------plot6-----------------------------------------------------------
%subfunction Plot Beschneistunden (mittlere Anzahl, max min, täglich)

plot6( tw_shortterm_year, jahr_start);
%-----------------------------------------plot6-----------------------------------------------------------


%----------------------------------------------------plot7-----------------------------------------------------------
%subfunction plot7 plottet mittleren Tagesgang für Oktober-April
plot7( tw2_dailyamp_mean, tw2_dailyamp_rel, tw2_dailyamp_prop, tw2_dailyamp_lan);
%----------------------------------------------------plot7-----------------------------------------------------------


%---------------------------------------------------plot8------------------
%Plot mittlerer Tagesgang Monatsweise
m9=figure('name','Mittlerer Tagesgang saison')
subplot(2,2,1)
set(gca,'fontsize',14)
plot(0:23,tw2_dailyamp_oct_2,'b',0:23,tw2_dailyamp_nov_2,'-db',0:23,tw2_dailyamp_dec_2,'-sr',0:23,tw2_dailyamp_jan_2,'-or',0:23,tw2_dailyamp_feb_2,'-+r',0:23,tw2_dailyamp_mar_2,'-*k',0:23,tw2_dailyamp_apr_2,'-.k')
xlabel('Uhrzeit (UTC)')
ylabel('Feuchttemperatur (°C)')
grid off
axis tight
legend('oct','nov','dec','jan','feb','mar','apr')
title('Mittlerer Tagesgang der Feuchttemperatur')

subplot(2,2,2)
set(gca,'fontsize',14)
bar(teste,'stacked')
set(gca,'XTickLabel','Oct|Nov|Dec|Jan|Feb|Mar|Apr')
colormap hot
legend('0-5 Uhr','6-11 Uhr','12-17 Uhr','18-23 Uhr')
grid off
axis tight
title('Zeitliche Verteilung der künstlich erzeugbaren Schneemenge')

subplot(2,2,3)
set(gca,'fontsize',14)
plot(0:23,tw2_dailyamp_prop_oct,'b',0:23,tw2_dailyamp_prop_nov,'-db',0:23,tw2_dailyamp_prop_dec,'-sr',0:23,tw2_dailyamp_prop_jan,'-or',0:23,tw2_dailyamp_prop_feb,'-+r',0:23,tw2_dailyamp_prop_mar,'-*k',0:23,tw2_dailyamp_prop_apr,'-.k')
xlabel('Uhrzeit (UTC)')
ylabel('Schneileistung (m³/h)')
grid on
axis tight
legend('oct','nov','dec','jan','feb','mar','apr')
title('Mittlerer Tagesgang der Schneileistung (Propellererzeuger)')

subplot(2,2,4)
set(gca,'fontsize',14)
plot(0:23,tw2_dailyamp_lan_oct,'b',0:23,tw2_dailyamp_lan_nov,'-db',0:23,tw2_dailyamp_lan_dec,'-sr',0:23,tw2_dailyamp_lan_jan,'-or',0:23,tw2_dailyamp_lan_feb,'-+r',0:23,tw2_dailyamp_lan_mar,'-*k',0:23,tw2_dailyamp_lan_apr,'-.k')
xlabel('Uhrzeit (UTC)')
ylabel('Schneileistung (m³/h)')
grid on
axis tight
title('Mittlerer Tagesgang der Schneileistung (Lanzenerzeuger)')

saveas(m9, 'Mittlerer Tagesgang Saison.jpg');
saveas(m9, 'Mittlerer Tagesgang Saison.fig');