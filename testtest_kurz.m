% %Einlesen der MATLAB-Formatierten DWD Meteodaten aus files
% clear all
% close all
% 
% %***************
% %EINSTELLUNGEN
% %***************
% %==========================================================================
% location='Obergurgl';
% alti=1940; %Seehöhe
%  jahr_start=1953;
%  jahr_end=2014;
%  seasNr=jahr_end-jahr_start;
%  %--------------------------------------------------!!!!!!!-------------       
%      
% %==========================================================================
% % Daten einlesen
% file1=fopen('Daten_ok\ZAMG\Obergurgl1953_2014Tage.txt');
% C=textscan(file1, '%f %f %f %f %f', 'headerlines',0);
% fclose(file1);
% %**************************************
% %PROGRAM SECTION -> nicht verändern!
% %**************************************
% 
% %calculate mean air pressure using simple barometric formula with altitude (for psychrometric formula) 
% press=1013.25.*exp((alti.*(-1))./7290);
% 
% yyyy=C{1};
% mon=C{2};
% dd=C{3};
% %careful here
% rel_mean_d=C{4};
% t_mean_d=C{5};

function[tw_prob_year_last20, date_year] = testtest_kurz(jahr_start, jahr_end, seasNr, press, yyyy, mon, dd, rel_mean_d, t_mean_d); 

Q=find(mon==2 & dd==29);
yyyy(Q)=[];
mon(Q)=[];
dd(Q)=[];
rel_mean_d(Q)=[];
t_mean_d(Q)=[];

HH=zeros(size(yyyy));
mins=zeros(size(yyyy));
sec=zeros(size(yyyy));
datumvec=[yyyy, mon, dd,HH, mins, sec];
datestring=datestr(datumvec, 'yyyy-mm-dd');
date_winter=datenum(datestring, 'yyyy-mm-dd');
date_daily=date_winter;
  

AA=find(yyyy==1991 & mon==10 & dd==1);
AAA=find(yyyy==1992 & mon==4 & dd==30);
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
    


%Calculate longterm wet-bulb temperatur - select Temperature mean and rh mean to use when calculating Tw
 [tw_mean_d]=wetbulb(t_mean_d, rel_mean_d, press);  
 
%if wetblub takes too long, write output to file once and then read tw from file. 
% data_1(:,1)=tw_daily_longterm;
% filename1 =['tw_longterm.out'];
% %File schreiben
% dlmwrite(filename1,data_1,' ')


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
    
%count number of days per season where tw is less than -2, -3, -4 and greater than -2.
    for i=1:seasNr;
     count_tw_2(i)=length(find(tw_winter(:,i)<=-2));
     count_tw_3(i)=length(find(tw_winter(:,i)<=-3));
     count_tw_4(i)=length(find(tw_winter(:,i)<=-4));
     count_tw_5(i)=length(find(tw_winter(:,i)>-2));
    end
%calculate standard deviation of count variables.    
    count_tw_2_std=std(count_tw_2);
    count_tw_3_std=std(count_tw_3);
    count_tw_4_std=std(count_tw_4);
    count_tw_5_std=std(count_tw_5);
 
     
%Trend/Rauschverhältnis berechnen, linearen Trend berechnen
      
 %signifikanten Trend finden; Trend ist dann signifikant, wenn Trend/Rauschverhältnis >1.64 (90% Confidence-Level)
jan_15=jan_15' ;

[Tr_max2, trend2, xxxcount2, ycount2]=trend_rausch(count_tw_2, jan_15);


        
%CONTOUR Plot mit Entwicklung Anzahl Tage/Feuchtemperatur in Jahren
max(tw_mean_d)
tw_mean_d;
min(tw_mean_d)
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
    tw_longterm_feb_seas,tw_count_feb, tw_longterm_mar_seas,tw_count_mar, tw_longterm_apr_seas,tw_count_apr, tab_tw_count]=monthly_tw(mon, tw_mean_d, seasNr);

  
   
%      distributionplots_kurz(location, yyyy, mon, dd, tw_mean_d, t_mean_d,rel_mean_d, press, tw_longterm_oct, tw_longterm_nov, tw_longterm_dec, tw_longterm_jan, tw_longterm_feb, tw_longterm_mar, tw_longterm_apr, U)

%Fliri-Probability Plot; Probability to have artificial snow conditions on a specific day. Call fliri_klima subfunction. Input: seasNr, number of
   %days per month, tw_longterm_seas :tw_count_oct_prob
   
   %comment last loop in fliri subfct. so that zeros are given out for
   %period where no data exists.
   
%Shorter Period (1987/88-2006/07) last 20 seasons to detect changes in probability : tw_count_oct_prob_last20
   
[tw_count_oct_prob, tw_count_oct_prob_last20, tw_count_oct_prob_last20_2]=fliri_klima(seasNr, 31, tw_longterm_oct_seas); 
[tw_count_nov_prob, tw_count_nov_prob_last20, tw_count_nov_prob_last20_2]=fliri_klima(seasNr, 30, tw_longterm_nov_seas); 
[tw_count_dec_prob, tw_count_dec_prob_last20, tw_count_dec_prob_last20_2]=fliri_klima(seasNr, 31, tw_longterm_dec_seas); 
[tw_count_jan_prob, tw_count_jan_prob_last20, tw_count_jan_prob_last20_2]=fliri_klima(seasNr, 31, tw_longterm_jan_seas); 
[tw_count_feb_prob, tw_count_feb_prob_last20, tw_count_feb_prob_last20_2]=fliri_klima(seasNr, 28, tw_longterm_feb_seas); 
[tw_count_mar_prob, tw_count_mar_prob_last20, tw_count_mar_prob_last20_2]=fliri_klima(seasNr, 31, tw_longterm_mar_seas); 
[tw_count_apr_prob, tw_count_apr_prob_last20, tw_count_apr_prob_last20_2]=fliri_klima(seasNr, 30, tw_longterm_apr_seas); 

   
%   Merge monthly probabilities to yearly plot
  tw_prob_year(1:31,1:4)=tw_count_oct_prob;
  tw_prob_year(32:61,1:4)=tw_count_nov_prob;
  tw_prob_year(62:92,1:4)=tw_count_dec_prob;
  tw_prob_year(93:123,1:4)=tw_count_jan_prob;
  tw_prob_year(124:151,1:4)=tw_count_feb_prob;
  tw_prob_year(152:182,1:4)=tw_count_mar_prob;
  tw_prob_year(183:212,1:4)=tw_count_apr_prob; 
  
 date_year=date_winter(1:212); %dates for one season for plot (takes first season in time series but year is not shown in plot)

  
  %Merge monthly probabilities to yearly plot
  tw_prob_year_last20(1:31,1:4)=tw_count_oct_prob_last20;
  tw_prob_year_last20(32:61,1:4)=tw_count_nov_prob_last20;
  tw_prob_year_last20(62:92,1:4)=tw_count_dec_prob_last20;
  tw_prob_year_last20(93:123,1:4)=tw_count_jan_prob_last20;
  tw_prob_year_last20(124:151,1:4)=tw_count_feb_prob_last20;
  tw_prob_year_last20(152:182,1:4)=tw_count_mar_prob_last20;
  tw_prob_year_last20(183:212,1:4)=tw_count_apr_prob_last20; 
    
  %Merge monthly probabilities to yearly plot
  tw_prob_year_last20_2(1:31,1:4)=tw_count_oct_prob_last20_2;
  tw_prob_year_last20_2(32:61,1:4)=tw_count_nov_prob_last20_2;
  tw_prob_year_last20_2(62:92,1:4)=tw_count_dec_prob_last20_2;
  tw_prob_year_last20_2(93:123,1:4)=tw_count_jan_prob_last20_2;
  tw_prob_year_last20_2(124:151,1:4)=tw_count_feb_prob_last20_2;
  tw_prob_year_last20_2(152:182,1:4)=tw_count_mar_prob_last20_2;
  tw_prob_year_last20_2(183:212,1:4)=tw_count_apr_prob_last20_2; 

%     filename22=['twprobyear',location,'.out'];
%   filename222=['twprobyearlast20',location,'.out'];
%   
%    dlmwrite(filename22,tw_prob_year,'delimiter', ' ', 'newline', 'pc');
%    dlmwrite(filename222,tw_prob_year_last20,'delimiter', ' ', 'newline', 'pc');
% % Write files with number of snowdays per year and per year and month.
% % Insert desired file name. 
jahr=U;
jahr(1)=[];
jahr=jahr';
data=[jahr; count_tw_2];   

 
% filename1=['schneitage',location, '.txt'];
% dlmwrite(filename1,data,'delimiter', ' ', 'newline', 'pc');
% 
% data_all(:,1)=tw_count_oct;
% data_all(:,2)=tw_count_nov;
% data_all(:,3)=tw_count_dec;
% data_all(:,4)=tw_count_jan;
% data_all(:,5)=tw_count_feb;
% data_all(:,6)=tw_count_mar;
% data_all(:,7)=tw_count_apr;
%   filename2=['schneitage',location, 'all.txt'];
% dlmwrite(filename2,data_all,'delimiter', ' ', 'newline', 'pc');
% 
% last20_2=nanmean(tw_prob_year_last20_2)
% last20=nanmean(tw_prob_year_last20)
% 
% 
% meanvalue(:,1)=season_mean_t;
% meanvalue(:,2)=season_mean_rh;
% meanvalue(:,3)=season_mean_tw;
% 
% filename3=['means_',location,'.txt'];
% dlmwrite(filename3,meanvalue,'delimiter', ' ' , 'newline', 'pc');
%  
% %   PLOTS
%  
% 
% m= figure('name','Trend Klima');
%    
%     subplot(3,1,1)
%     set(gca,'fontsize',14)
%     plot(jan_15,season_mean_t,'r','linewidth',2)
%     grid on
%     axis tight
%     ylabel('Temperatur (°C)')
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     title(['Saisonale Temperatur ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold');
%     
%      subplot(3,1,2)
%     set(gca,'fontsize',14)
%     plot(jan_15,season_mean_tw,'k','linewidth',2)
%     grid on
%     axis tight
%     ylabel('Feuchttemperatur (°C)')
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     title(['Saisonale Feuchttemperatur ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold');
%     
%       subplot(3,1,3)
%     set(gca,'fontsize',14)
%     plot(jan_15,season_mean_rh,'b','linewidth',2)
%     grid on
%     axis tight
%     ylabel('Relative Feuchte (%)')
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     title(['Saisonale Feuchte ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold');
%     
% saveas(m, 'fig3.jpg');
% saveas(m, 'fig3.fig');
% 
%    %--------------------------------------------------------------------- 
%   mm=  figure('name','Number of days and Trend');
%   
%     set(gca,'fontsize',14)
%     plot(jan_15,count_tw_2,jan_15,count_tw_3,'-+r',jan_15,count_tw_4,'-*k',jan_15,count_tw_5,'-og')
%     grid on
%     axis tight
%     ylabel('Anzahl Tage')
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     title(['Anzahl der Tage pro Saison mit Tagesmittel d. Feuchttemperatur, ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
%     legend('\leq -2°C','\leq -3°C','\leq -4°C','> -2°C')    
%        
% saveas(m, 'fig4.jpg');
% saveas(m, 'fig4.fig');
% 
% 
% 
% 
% %          
%     %Contour Plot
%    m1= figure('name','Entwicklung Tf');
%     [C,h]=contourf(jan_15,a3,ZII,50);
%     colormap(jet)
%     set(h,'LineStyle','none');
%     set(gca,'Layer','bottom')
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     colorbar
%     title('Entwicklung der Feuchttemperatur')
%     set(gca,'fontsize',14)
%     axis tight
%     xlabel('Jahre')
%     ylabel('Feuchttemperatur (°C)')
%     hold on
%   %  bar(jan_15,diff_clim_tw,'k')
% saveas(m1, 'fig5.jpg');
% saveas(m1, 'fig5.fig');
% 
%     
% m2=    figure('name','Trend number of days monthly');
%     subplot(3,1,1)
%     set(gca,'fontsize',14)
%     plot(jan_15,tw_count_oct,'r',jan_15,tw_count_nov,'-.k')
%     grid off
%     axis tight
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     legend('Oct','Nov','location','best')
%     xlabel('Jahre')
%     %ylabel('Anomalie d. Tagesanzahl mit Tagesmittel Tf \leq -2°C')
%     title('Vorsaison (Oktober & November)')
%     
%     subplot(3,1,2)
%     set(gca,'fontsize',14)
%     plot(jan_15,tw_count_dec,'r',jan_15,tw_count_jan,'-.k',jan_15,tw_count_feb,':b')
%     grid off
%     axis tight
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     legend('Dez','Jan','Feb','location','best')
%     xlabel('Jahre')
%     ylabel('Tagesanzahl mit Tagesmittel Tf  \leq -2°C')
%     title('Hauptsaison (Dezember - Februar)')
%     
%     subplot(3,1,3)
%     set(gca,'fontsize',14)
%     plot(jan_15,tw_count_mar,'r',jan_15,tw_count_apr,'-.k')
%     grid off
%     axis tight
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     legend('Mar','Apr','location','best')
%     xlabel('Jahre')
%    % ylabel('Tagesanzahl mit Tagesmittel Tf  \leq -2°C')
%     title('Nachsaison (März - April)')
% saveas(m2, 'fig6.jpg');
% saveas(m2, 'fig6.fig');
% 
% 
% 
% dateplot=[date_daily(AA):date_daily(AAA-1)];    
% 
%  m6=   figure('name','Probability Year');
%     [haxes,hline1,hline2]=plotyy(date_year,tw_prob_year(:,1),date_year,tw_prob_year(:,2),'area','plot');
%     grid off
%     set(hline1,'FaceColor',[0.8 0.8 0.8])
%     set(hline2,'Color','b')
%     axes(haxes(1))
%     set(gca,'fontsize',14)
%     ylabel('Wahrscheinlichkeit (%)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
%     set(haxes(2), 'Box', 'Off');
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
%     plot(date_year,tw_prob_year(:,3),'r')
%     hold on
%     plot(date_year,tw_prob_year(:,4),'k')
%     title(['Wahrscheinlichkeit f. Tagesmittelwert Feuchttemp. \leq -2°C, Mittelwert und Extremwerte (basierend auf ',num2str(jahr_start),' - ',num2str(jahr_end),' ) ,',num2str(location),' ',num2str(alti),' m']);
%     axis tight
%     hline=refline(0,-2);
%     set(hline, 'Color', 'k', 'LineWidth', 2);
%     
%     saveas(m6, 'fig10.jpg');
% saveas(m6, 'fig10.fig');
%     
%     
%     
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
%     set(haxes(2), 'Box', 'Off');
%     set(haxes(1),'Box','Off') 
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
%     axis tight
%     axes(haxes(2))
%     set(gca,'fontsize',14,'YColor','k')
%     ylabel('Extremwerte und Mittelwert Feuchttemperatur (°C)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     set(gca,'YTick',[-25 -20 -15 -10 -5 0 5])
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
%     hold on
%     plot(date_year,tw_prob_year_last20(:,3),'r')
%     hold on
%     plot(date_year,tw_prob_year(:,4),'k')
%     title(['Wahrscheinlichkeit f. Tagesmittelwert Feuchttemp. \leq -2°C, Mittelwert und Extremwerte (basierend auf 1993-2014), ', num2str(location),' ',num2str(alti)])
%     axis tight
%     hline=refline(0,-2);
%     set(hline, 'Color', 'k', 'LineWidth', 2);
%     
%         saveas(m7, 'fig11.jpg');
% saveas(m7, 'fig11.fig');
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
