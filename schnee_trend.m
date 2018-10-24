%Einlesen der MATLAB-Formatierten DWD Meteodaten aus files
clear all
close all

%***************
%EINSTELLUNGEN
%***************
%==========================================================================
%Tageswerte Schneehöhe 

%ohne februar 29.

location='HPeissenberg';
alti=977; %Seehöhe
 jahr_start=1901;
 jahr_end=2014;
 seasNr=jahr_end-jahr_start;
%  dir='snow';
 %--------------------------------------------------!!!!!!!-------------       
     
%==========================================================================
% Daten einlesen
file1=fopen('SnowHPeissenberg\mehr als 30cm und mittelHPeissenberg.txt');
C=textscan(file1, '%f %f %f', 'headerlines',0);
fclose(file1);

%**************************************
%PROGRAM SECTION -> nicht verändern!
%**************************************

yyyy=C{1};
tage30plus=C{2};
mittel=C{3};

%ref periode 1960 bis 1990
pp=find(yyyy==1960);
pp2=find(yyyy==1990);

%----------------------------------------------------------------------
k=212; %Number of days per winter season (Oct-Apr)


%calculate standard deviation of count variables.    
    count_SH30_std=std(tage30plus);

%Calculate climatological reference period 60/61-90/91 out of seasonal means  
    clim_ref_count_SH30=round(nanmean(tage30plus(pp:pp2)));
    
%Abweichung der saisonalen Mittel vom klimatologischen Mittel berechnen
    diff_clim_count_SH30=tage30plus-clim_ref_count_SH30;
      
%Calculate 10-season running-mean
 
 [diff_clim_count_SH30_10]=running_10(diff_clim_count_SH30);   
     
%Trend/Rauschverhältnis berechnen, linearen Trend berechnen
      
 %signifikanten Trend finden; Trend ist dann signifikant, wenn Trend/Rauschverhältnis >1.64 (90% Confidence-Level)
 jan_15=yyyy ;

[Tr_max2, trend2, xxxcount2, ycount2]=trend_rausch_snow(diff_clim_count_SH30, jan_15);
Tr_max2
trend2
% 
% %---------------------------section seasonal means over-------------------
% %    Trendanalyse SH
   m= figure('name','Schneehöhe');
   
    subplot(2,1,1)
    set(gca,'fontsize',14)
    plot(jan_15,mittel,'r','linewidth',2)
    grid on
    axis tight
    ylabel('Schneehöhe (cm)')
    set(gca,'XTick',jan_15(1:10:end))
    set(gca,'XTickLabel',jan_15(1:10:end))
%     datetick('x','keeplimits','keepticks')
    title(['Saisonale Schneehöhe, ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    bar(jan_15,tage30plus,'b')
    grid on
    axis tight
    ylabel('Anzahl an Tagen mit SH  \geq 30 cm')
    set(gca,'XTick',jan_15(1:10:end))
    set(gca,'XTickLabel',jan_15(1:10:end))
    title(['Tage mit SH über 30 cm, ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
  
    
m2 = figure('name','Trend Tage');
    set(gca,'fontsize',14)
    hold on
    bar(jan_15,diff_clim_count_SH30,'b')
    plot(jan_15,diff_clim_count_SH30_10,'k','linewidth',2)
    if Tr_max2>1.64
%     size(ycount2)
%     size(xxxcount2)
% ycount2(:,xxxcount2)
    plot(jan_15,ycount2(:,xxxcount2),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Tage mit SH über 30 cm (Anomalie)')
    axis([min(jan_15) max(jan_15) fix(min(diff_clim_count_SH30)) round(max(diff_clim_count_SH30))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1),ymittl+2,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH30),-1)),' Tage'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+20,ymittl+1.5,['Klimatol. Mittel (1960-90): ',num2str(roundn(clim_ref_count_SH30,-1)),' Tage'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+40,ymittl+2,['T/R.: ',num2str(roundn((Tr_max2),-2)),' (',num2str(jan_15(xxxcount2)),'/',num2str(jan_15(xxxcount2+1)),' - ',num2str(jan_15(end-1)),'/',num2str(jan_15(end)),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_max2>1.64
    text(xpos(1)+40,ymittl+1.5,['Linearer Trend: ',num2str(roundn(trend2,-1)),'Tage/a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(1:10:end))
    set(gca,'XTickLabel',jan_15(1:10:end))
    set(gca,'fontsize',14)
    title(['Abweichung der Tage mit SH \geq 30cm ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.) vom langjährigen Mittel (1960-90), ',num2str(location)],'fontsize',12,'fontweight','bold')
    %line([jan_15(47)-175 jan_15(47)-175],[min(diff_clim_t) max(diff_clim_t)+2],'linewidth',2,'color',[0.7 0.7 0.7])
    



% 
%     
% m2=    figure('name','Trend number of days monthly');
%     subplot(3,1,1)
%     set(gca,'fontsize',14)
%     %plot(jan_15,SH_count_oct,'r',jan_15,SH_count_nov,'-.k')
%     plot(jan_15,SH_oct,'r',jan_15,SH_nov,'-.k')
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
%     %plot(jan_15,SH_count_dec,'r',jan_15,SH_count_jan,'-.k',jan_15,SH_count_feb,':b')
%     plot(jan_15,SH_dec,'r',jan_15,SH_jan,'-.k',jan_15,SH_feb,':b')
%     grid off
%     axis tight
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     legend('Dez','Jan','Feb','location','best')
%     xlabel('Jahre')
%     %ylabel('Tagesanzahl mit SH >= 30cm')
%     ylabel('Schneehöhe - Monatsmittel (cm)')
%     title('Hauptsaison (Dezember - Februar)')
%     
%     subplot(3,1,3)
%     set(gca,'fontsize',14)
% %     plot(jan_15,SH_count_mar,'r',jan_15,SH_count_apr,'-.k')
%     plot(jan_15,SH_mar,'r',jan_15,SH_apr,'-.k')
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
%     
% m3=    figure('name','Anomalie d. Anzahl an Tagen nach Monat: Vorsaison');
% %    Trend number of days
%     
%     subplot(2,1,1)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_count_SH_oct,'b')
%     plot(jan_15,diff_clim_count_SH_oct_10,'k','linewidth',2)
%     if Tr_count_SH_oct_max>1.64
%     plot(jan_15,y8(:,xxx8),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Oktober (Tage)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_oct)) round(max(diff_clim_count_SH_oct))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_oct),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_oct_max,-2)),' (',num2str(datestr(jan_15(xxx8),'yyyy')),'/',num2str(datestr(jan_15(xxx8+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_count_SH_oct_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_oct,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',14)
%      title('Vorsaison (Oktober & November) Anomalie d. Anzahl an Tagen mit SH \geq 30cm')
%     
%     subplot(2,1,2)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_count_SH_nov,'b')
%     plot(jan_15,diff_clim_count_SH_nov_10,'k','linewidth',2)
%     if Tr_count_SH_nov_max>1.64
%     plot(jan_15,y9(:,xxx9),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('November (Tage)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_nov)) round(max(diff_clim_count_SH_nov))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_nov),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_nov_max,-2)),' (',num2str(datestr(jan_15(xxx9),'yyyy')),'/',num2str(datestr(jan_15(xxx9+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_count_SH_nov_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_nov,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     
% saveas(m3, 'fig7.jpg');
% saveas(m3, 'fig7.fig');
%     
% %     
% 
% %für schneehöhenanomalie, nicht anzahl an tagen. 
% m3=    figure('name','Anomalie d. Anzahl an Tagen nach Monat: Vorsaison');
% %    Trend number of days
%     
%     subplot(2,1,1)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_SH_oct,'b')
%     plot(jan_15,diff_clim_SH_oct_10,'k','linewidth',2)
%     if Tr_SH_oct_max>1.64
%     plot(jan_15,y8(:,xxx8),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Oktober (cm)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_oct)) round(max(diff_clim_SH_oct))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Kilm. Mittel.: ',num2str(clim_ref_SH_oct)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_oct_max,-2)),' (',num2str(datestr(jan_15(xxx8),'yyyy')),'/',num2str(datestr(jan_15(xxx8+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_SH_oct_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_oct,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',14)
%      title('Vorsaison (Oktober & November) Anomalie d. Schneehöhe im Vergleich zum langjährigen Mittel.')
%     
%     subplot(2,1,2)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_SH_nov,'b')
%     plot(jan_15,diff_clim_SH_nov_10,'k','linewidth',2)
%     if Tr_SH_nov_max>1.64
%     plot(jan_15,y9(:,xxx9),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('November (cm)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_nov)) round(max(diff_clim_SH_nov))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Klim. Mittel: ',num2str(clim_ref_SH_nov)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_nov_max,-2)),' (',num2str(datestr(jan_15(xxx9),'yyyy')),'/',num2str(datestr(jan_15(xxx9+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_SH_nov_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_nov,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     
% saveas(m3, 'fig7.jpg');
% saveas(m3, 'fig7.fig');
%     
% 
% 
% 
% 
%     
%     
%  m4=   figure('name','Anomalie d. Anzahl an Tagen nach Monat: Hauptsaison');
%   %  Trend number of days
%     
%     subplot(3,1,1)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_count_SH_dec,'b')
%     plot(jan_15,diff_clim_count_SH_dec_10,'k','linewidth',2)
%     if Tr_count_SH_dec_max>1.64
%     plot(jan_15,y10(:,xxx10),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Dezember')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_dec)) round(max(diff_clim_count_SH_dec))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_dec),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_dec_max,-2)),' (',num2str(datestr(jan_15(xxx10),'yyyy')),'/',num2str(datestr(jan_15(xxx10+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_count_SH_dec_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_dec,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%    end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     title('Hauptsaison (Dezember, Januar, Februar) Anomalie d. Anzahl an Tagen mit SH \geq 30cm');
%     
%     
%     subplot(3,1,2)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_count_SH_jan,'b')
%     plot(jan_15,diff_clim_count_SH_jan_10,'k','linewidth',2)
%    if Tr_count_SH_jan_max>1.64
%     plot(jan_15,y11(:,xxx11),'-.r','linewidth',2)
%    end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Januar')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_jan)) round(max(diff_clim_count_SH_jan))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_jan),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_jan_max,-2)),' (',num2str(datestr(jan_15(xxx11),'yyyy')),'/',num2str(datestr(jan_15(xxx11+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_count_SH_jan_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_jan,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',14)
%         
%     subplot(3,1,3)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_count_SH_feb,'b')
%     plot(jan_15,diff_clim_count_SH_feb_10,'k','linewidth',2)
%     if Tr_count_SH_feb_max>1.64
%     plot(jan_15,y12(:,xxx12),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Februar')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_feb)) round(max(diff_clim_count_SH_feb))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_feb),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_feb_max,-2)),' (',num2str(datestr(jan_15(xxx12),'yyyy')),'/',num2str(datestr(jan_15(xxx12+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%    if Tr_count_SH_feb_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_feb,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%    end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     
%     saveas(m4, 'fig8.jpg');
% saveas(m4, 'fig8.fig');
%     
%     
% 
%  m4=   figure('name','Anomalie d. Anzahl an Tagen nach Monat: Hauptsaison');
%   %  Trend number of days
%     
%     subplot(3,1,1)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_SH_dec,'b')
%     plot(jan_15,diff_clim_SH_dec_10,'k','linewidth',2)
%     if Tr_SH_dec_max>1.64
%     plot(jan_15,y10(:,xxx10),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Dezember (cm)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_dec)) round(max(diff_clim_SH_dec))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Klim. Mittel.: ',num2str(clim_ref_SH_dec)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_dec_max,-2)),' (',num2str(datestr(jan_15(xxx10),'yyyy')),'/',num2str(datestr(jan_15(xxx10+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_SH_dec_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_dec,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%    end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     title('Hauptsaison (Dezember, Januar, Februar) Anomalie d. Schneehöhe im Vergleich zum langjährigen Mittel');
%    
%     
%     subplot(3,1,2)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_SH_jan,'b')
%     plot(jan_15,diff_clim_SH_jan_10,'k','linewidth',2)
%    if Tr_SH_jan_max>1.64
%     plot(jan_15,y11(:,xxx11),'-.r','linewidth',2)
%    end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Januar (cm)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_jan)) round(max(diff_clim_SH_jan))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Klim. Mittel.: ',num2str(clim_ref_SH_jan)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_jan_max,-2)),' (',num2str(datestr(jan_15(xxx11),'yyyy')),'/',num2str(datestr(jan_15(xxx11+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_SH_jan_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_jan,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',14)
%         
%     subplot(3,1,3)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_SH_feb,'b')
%     plot(jan_15,diff_clim_SH_feb_10,'k','linewidth',2)
%     if Tr_SH_feb_max>1.64
%     plot(jan_15,y12(:,xxx12),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Februar (cm)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_feb)) round(max(diff_clim_SH_feb))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Klim. Mittel.: ',num2str(clim_ref_SH_feb)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_feb_max,-2)),' (',num2str(datestr(jan_15(xxx12),'yyyy')),'/',num2str(datestr(jan_15(xxx12+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%    if Tr_SH_feb_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_feb,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%    end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     
%     saveas(m4, 'fig8.jpg');
% saveas(m4, 'fig8.fig');
% 
% 
% 
% 
% 
% 
%     
%   m5=  figure('name','Anomalie d. Anzahl an Tagen nach Monat: Nachsaison');
%    
%     subplot(2,1,1)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_count_SH_mar,'b')
%     plot(jan_15,diff_clim_count_SH_mar_10,'k','linewidth',2)
%     if Tr_count_SH_mar_max>1.64
%     plot(jan_15,y13(:,xxx13),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('März')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_mar)) round(max(diff_clim_count_SH_mar))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_mar),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_mar_max,-2)),' (',num2str(datestr(jan_15(xxx13),'yyyy')),'/',num2str(datestr(jan_15(xxx13+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_count_SH_mar_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_mar,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     title('Nachsaison (März & April)Anomalie d. Anzahl an Tagen mit SH \geq 30cm')
%     
%     subplot(2,1,2)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_count_SH_apr,'b')
%     plot(jan_15,diff_clim_count_SH_apr_10,'k','linewidth',2)
%     if Tr_count_SH_apr_max>1.64
%     plot(jan_15,y14(:,xxx14),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('April')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_SH_apr)) round(max(diff_clim_count_SH_apr))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_SH_apr),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_count_SH_apr_max,-2)),' (',num2str(datestr(jan_15(xxx14),'yyyy')),'/',num2str(datestr(jan_15(xxx14+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_count_SH_apr_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_count_SH_apr,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     
%     
%     saveas(m5, 'fig9.jpg');
% saveas(m5, 'fig9.fig');
% 
% 
%     
%   m5=  figure('name','Anomalie d. Anzahl an Tagen nach Monat: Nachsaison');
%    
%     subplot(2,1,1)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_SH_mar,'b')
%     plot(jan_15,diff_clim_SH_mar_10,'k','linewidth',2)
%     if Tr_SH_mar_max>1.64
%     plot(jan_15,y13(:,xxx13),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('März (cm)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_mar)) round(max(diff_clim_SH_mar))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Klim. Mittel: ',num2str(clim_ref_SH_mar)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_mar_max,-2)),' (',num2str(datestr(jan_15(xxx13),'yyyy')),'/',num2str(datestr(jan_15(xxx13+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_SH_mar_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_mar,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     title('Nachsaison (März & April) Anomalie d. Schneehöhe im Vergleich zum langjährigen Mittel')
%     
%     subplot(2,1,2)
%     set(gca,'fontsize',12,'fontweight','bold')
%     hold on
%     bar(jan_15,diff_clim_SH_apr,'b')
%     plot(jan_15,diff_clim_SH_apr_10,'k','linewidth',2)
%     if Tr_SH_apr_max>1.64
%     plot(jan_15,y14(:,xxx14),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('April (cm)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_SH_apr)) round(max(diff_clim_SH_apr))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl,['Klim. Mittel: ',num2str(clim_ref_SH_mar)],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_SH_apr_max,-2)),' (',num2str(datestr(jan_15(xxx14),'yyyy')),'/',num2str(datestr(jan_15(xxx14+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_SH_apr_max>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trend_SH_apr,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',12,'fontweight','bold')
%     
%     
%     saveas(m5, 'fig9.jpg');
% saveas(m5, 'fig9.fig');
% 
% 
% dateplot=[date_daily(AA):date_daily(AAA-1)];    
% 
%  m6=   figure('name','Probability Year');
%     [haxes,hline1,hline2]=plotyy(date_year,SH_prob_year(:,1),date_year,SH_prob_year(:,2),'area','plot');
%     grid off
%     set(hline1,'FaceColor',[0.8 0.8 0.8])
%     set(hline2,'Color','b')
%     axes(haxes(1))
%     set(gca,'fontsize',14)
%     ylabel('Wahrscheinlichkeit (%)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
%     axis tight
%     axes(haxes(2))
%     set(gca,'fontsize',14,'YColor','k')
%     ylabel('Extremwerte Schneehöhe (cm)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     set(gca,'YTick',[-25 -20 -15 -10 -5 0 5])
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
%     hold on
%     plot(date_year,SH_prob_year(:,3),'r')
%     title(['Wahrscheinlichkeit f. Schneehöhe \geq 30 cm, und Extremwerte (basierend auf ',num2str(jahr_start),' - ',num2str(jahr_end),' ) ,',num2str(location),' ',num2str(alti),' m']);
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
%     [haxes,hline1,hline2]=plotyy(date_year,SH_prob_year_last20(:,1),date_year,SH_prob_year_last20(:,2),'area','plot');
%     grid off
%     set(hline1,'FaceColor',[0.8 0.8 0.8])
%     set(hline2,'Color','b')
%     axes(haxes(1))
%     set(gca,'fontsize',14)
%     ylabel('Wahrscheinlichkeit (%)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     set(gca,'YTick',[10 20 30 40 50 60 70 80 90 100])
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
%     axis tight
%     axes(haxes(2))
%     set(gca,'fontsize',14,'YColor','k')
%     ylabel('Extremwerte Schneehöhe (cm)')
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     set(gca,'YTick',[-25 -20 -15 -10 -5 0 5])
%     datetick('x','mmm','keeplimits','keepticks')
%     xlabel('Monate')
%     hold on
%     plot(date_year,SH_prob_year_last20(:,3),'r')
%     title(['Wahrscheinlichkeit f. Tagesmittelwert Schneehöhe \geq 30 cm, und Extremwerte (basierend auf 1993-2014), ', num2str(location),' ',num2str(alti)])
%     axis tight
%     hline=refline(0,-2);
%     set(hline, 'Color', 'k', 'LineWidth', 2);
%     
%         saveas(m7, 'fig11.jpg');
% saveas(m7, 'fig11.fig');
%     
%     
%   m8=  figure('name','Comparison Probabilities last 20 years/20 years before');
%     area(date_year,SH_prob_year_last20(:,1),'facecolor','g')
%     hold on
%     area(date_year,SH_prob_year_last20_2(:,1))
%     grid on
%     ylabel('Wahrscheinlichkeit')
%     xlabel('Monate')
%     axis tight
%     set(gca,'XTick',[datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
%     datetick('x','mmm','keeplimits','keepticks')
%     jahr_20=[num2str(jahr_end-20),'/ ',num2str(jahr_end-1919),'-',num2str(jahr_end-1),'/ ',num2str(jahr_end-2000)];
%     jahr_40=[num2str(jahr_end-40),'/ ',num2str(jahr_end-1939),'-',num2str(jahr_end-21),'/ ',num2str(jahr_end-1940)];
%     %legend('1987/88-2006/07','1967/68-1986/87')
%     legend(jahr_20,jahr_40)
%     
%         saveas(m8, 'fig12.jpg');
% saveas(m8, 'fig12.fig');
% % % % % % % end
% % % % % 
% % % % % toc
% % % % % t=toc;
% % % % % 
% % % % %display(t)
% % % % 
% % % % 
% % % % 
