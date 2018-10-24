function[]=plot_trendTW(season_mean_tw, jan_15, jahr_start, jahr_end, location, alti, diff_clim_tw, diff_clim_tw_10, xxx3, Tr_tw_max, clim_ref_tw, trend_tw, y3, ...
    count_tw_2, count_tw_3, count_tw_4, count_tw_5,diff_clim_count_tw_2,diff_clim_count_tw_2_10, y4, xxx4, Tr_count_tw_2_max, trend_count_tw_2 ); 


   m= figure('name','Trend wet-bulb temp');
   
    subplot(2,1,1)
    set(gca,'fontsize',14)
    plot(jan_15,season_mean_tw,'r','linewidth',2)
    grid on
    axis tight
    ylabel('Feuchttemperatur (°C)')
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    title(['Saisonale Feuchttemperatur ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
%     if nofile~=1
%     line([datum(1) datum(1)],[min(season_mean_tw) max(season_mean_tw)],'linewidth',2,'color',[0.7 0.7 0.7])
%     else
%       line([datum_clim(1) datum_clim(1)],[min(season_mean_tw) max(season_mean_tw)],'linewidth',2,'color',[0.7 0.7 0.7])  
%     end
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    hold on
    bar(jan_15,diff_clim_tw,'b')
    plot(jan_15,diff_clim_tw_10,'k','linewidth',2)
    if Tr_tw_max>1.64
    plot(jan_15,y3(:,xxx3),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Feuchttemperaturanomalie (°C)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_tw)) round(max(diff_clim_tw))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl+1,['Stabw.: ',num2str(roundn(nanstd(diff_clim_tw),-1)),' °C'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl+1,['T/R.: ',num2str(roundn(Tr_tw_max,-2)),' (',num2str(datestr(jan_15(xxx3),'yyyy')),'/',num2str(datestr(jan_15(xxx3+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    text(xpos(1)+500,ymittl+0.5,['Klimatol. Mittel (1960-90): ',num2str(roundn(clim_ref_tw,-1)),' °C'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  

    if Tr_tw_max>1.64
    text(xpos(1)+15500,ymittl+1,['Linearer Trend: ',num2str(roundn(trend_tw,-2)),'°C/ a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',14)
    title(['Abweichung der saisonalen Feuchttemperatur ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.) vom langjährigen Mittel (1960-90), ',num2str(location)],'fontsize',12,'fontweight','bold')
%     line([jan_15(47)-175 jan_15(47)-175],[min(diff_clim_tw) max(diff_clim_tw)+2],'linewidth',2,'color',[0.7 0.7 0.7])

saveas(m, 'fig3.jpg');
saveas(m, 'fig3.fig');

%    %--------------------------------------------------------------------- 
%   mm=  figure('name','Number of days and Trend');
%     subplot(2,1,1)
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
%     subplot(2,1,2)
%     set(gca,'fontsize',14)
%     hold on
%     bar(jan_15,diff_clim_count_tw_2,'b')
%     plot(jan_15,diff_clim_count_tw_2_10,'k','linewidth',2)
%     if Tr_count_tw_2_max>1.64
%     plot(jan_15,y4(:,xxx4),'-.r','linewidth',2)
%     end
%     axis tight
%     grid off
%     xlabel('Jahre')
%     ylabel('Anomalie d. Anzahl an Tagen (Tagesmittel Tf   \leq -2°C)')
%     axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_2)) round(max(diff_clim_count_tw_2))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     text(xpos(1)+500,ymittl-25,['Stabw.: ',num2str(roundn(nanstd(diff_clim_count_tw_2),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl-25,['T/R.: ',num2str(roundn(Tr_count_tw_2_max,-2)),' (',num2str(datestr(jan_15(xxx4),'yyyy')),'/',num2str(datestr(jan_15(xxx4+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_count_tw_2_max>1.64
%     text(xpos(1)+15500,ymittl-25,['Linearer Trend: ',num2str(roundn(trend_count_tw_2,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
%     set(gca,'XTick',jan_15(2:10:end))
%     set(gca,'XTickLabel',jan_15(2:10:end))
%     datetick('x','keeplimits','keepticks')
%     set(gca,'fontsize',14)
%    %title(['Abweichung der Tagesanzahl mit Tagesmittel d. Feuchttemperatur pro Saison (Okt-Apr) \leq -2°C ',num2str(jahr_start_2),' - ',num2str(jahr_end),' vom langjährigen Mittel (1960-90), ',num2str(location)],'fontsize',12,'fontweight','bold')
%    
% saveas(mm, 'fig4.jpg');
% saveas(mm, 'fig4.fig');