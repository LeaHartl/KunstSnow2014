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

location='Zugspitze';
alti=2964; %Seehöhe
 jahr_start=1950;
 jahr_end=2014;
 seasNr=jahr_end-jahr_start;
 %--------------------------------------------------!!!!!!!-------------       
     
%==========================================================================
% Daten einlesen
file1=fopen('Daten_ok\DWD\Zugspitze_Tage.out');
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


t_mean_d_2030=t_mean_d+1;
t_mean_d_2050=t_mean_d+1.8;


[tw_prob_year_last20_normal, date_year] = testtest_kurz(jahr_start, jahr_end, seasNr, press, yyyy, mon, dd, rel_mean_d, t_mean_d); 
[tw_prob_year_last20_2030, date_year] = testtest_kurz(jahr_start, jahr_end, seasNr, press, yyyy, mon, dd, rel_mean_d, t_mean_d_2030); 
[tw_prob_year_last20_2050, date_year] = testtest_kurz(jahr_start, jahr_end, seasNr, press, yyyy, mon, dd, rel_mean_d, t_mean_d_2050); 


 m6=   figure('name','Probability Year');
    h=area(date_year,tw_prob_year_last20_normal(:,1));
    grid off
    set(h,'FaceColor',[0.8 0.8 0.8])
    hold on
    h=plot(date_year,tw_prob_year_last20_2030(:,1), 'm:', 'Linewidth', 2);
    hold on
    h=plot(date_year,tw_prob_year_last20_2050(:,1), 'b:', 'Linewidth', 2);
% %     set(hline2,'Color','b')
%     axes(haxes(1))
     set(gca,'fontsize',14)
     ylabel('Wahrscheinlichkeit (%)')
     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
     set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
%     set(haxes(2), 'Box', 'Off');
%     set(haxes(1),'Box','Off') 
title([num2str(location),' ', num2str(alti), 'm']);
     datetick('x','mmm','keeplimits','keepticks')
     xlabel('Monate')
legend('1994-2014', '2030', '2050');
     
     %     axis tight
%     axes(haxes(2))
%     set(gca,'fontsize',14,'YColor','k')
%     ylabel('Extremwerte Feuchttemperatur (°C)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
% %     set(gca,'YTick',[-25 -20 -15 -10 -5 0 5])
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
% %     hold on
% %     plot(date_year,tw_prob_year(:,3),'r')
% %     hold on
% %     plot(date_year,tw_prob_year(:,4),'k')
    
%     axis tight
%     hline=refline(0,-2);
%     set(hline, 'Color', 'k', 'LineWidth', 2);
%     
%     saveas(m6, 'fig10.jpg');
% saveas(m6, 'fig10.fig');
