function[]=plot_trendT(season_mean, jan_15, jahr_start, jahr_end, location, alti, diff_clim_t, diff_clim_t_10, xxx, Tr_max, clim_ref_t, trend, y); 

%Trendanalyse Lufttemperatur und RH
   m= figure('name','Trend air temperature');
   
    subplot(2,1,1)
    set(gca,'fontsize',14)
    plot(jan_15,season_mean,'r','linewidth',2)
    grid on
    axis tight
    ylabel('Temperatur (°C)')
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    title(['Saisonale Lufttemperatur ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
    %if nofile~=1
%    line([datum(1) datum(1)],[min(season_mean) max(season_mean)],'linewidth',2,'color',[0.7 0.7 0.7])
    %else
    %  line([datum_clim(1) datum_clim(1)],[min(season_mean) max(season_mean)],'linewidth',2,'color',[0.7 0.7 0.7])  
    %end
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    hold on
    bar(jan_15,diff_clim_t,'b')
    plot(jan_15,diff_clim_t_10,'k','linewidth',2)
    %if Tr_max>1.64
    size(y)
    size(xxx)
    plot(jan_15,y(:,xxx),'-.r','linewidth',2)
    %end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Temperaturanomalie (°C)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_t)) round(max(diff_clim_t))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+500,ymittl+2,['Stabw.: ',num2str(roundn(nanstd(diff_clim_t),-1)),' °C'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+500,ymittl+1.5,['Klimatol. Mittel (1960-90): ',num2str(roundn(clim_ref_t,-1)),' °C'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl+2,['T/R.: ',num2str(roundn((Tr_max),-2)),' (',num2str(datestr(jan_15(xxx),'yyyy')),'/',num2str(datestr(jan_15(xxx+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    if Tr_max>1.64
    text(xpos(1)+15500,ymittl+2,['Linearer Trend: ',num2str(roundn(trend,-1)),'°C/a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',14)
    title(['Abweichung der saisonalen Lufttemperatur ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.) vom langjährigen Mittel (1960-90), ',num2str(location)],'fontsize',12,'fontweight','bold')
    %line([jan_15(47)-175 jan_15(47)-175],[min(diff_clim_t) max(diff_clim_t)+2],'linewidth',2,'color',[0.7 0.7 0.7])
    
saveas(m, 'fig1.jpg');
saveas(m, 'fig1.fig');
