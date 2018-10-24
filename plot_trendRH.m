function[]=plot_trendRH(season_mean_rh, jan_15, jahr_start, jahr_end, location, alti, diff_clim_rh, diff_clim_rh_10, xxx2, Tr_rh_max, clim_ref_rh, trend_rh, y2); 


    
    m=figure('name','Trend relative humidity');
   
    subplot(2,1,1)
    set(gca,'fontsize',14)
    plot(jan_15,season_mean_rh,'r','linewidth',2)
    grid on
    axis tight
    ylabel('Relative Feuchte (%)')
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    title(['Saisonale Relative Feuchte ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.), ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
%     if nofile~=1
%     line([datum(1) datum(1)],[min(season_mean_rh) max(season_mean_rh)],'linewidth',2,'color',[0.7 0.7 0.7])
%     else
%       line([datum_clim(1) datum_clim(1)],[min(season_mean_rh) max(season_mean_rh)],'linewidth',2,'color',[0.7 0.7 0.7])  
%     end
    
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    hold on
    bar(jan_15,diff_clim_rh,'b')
    plot(jan_15,diff_clim_rh_10,'k','linewidth',2)
    if Tr_rh_max>1.64
    plot(jan_15,y2(:,xxx2),'-.r','linewidth',2)
    end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('Feuchteanomalie (%)')
    axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_rh)) round(max(diff_clim_rh))])
    xpos=get(gca,'xlim');
    xmittl=(xpos(1)+xpos(2))./2;
    ypos=get(gca,'ylim');
    ymittl=(ypos(1)+ypos(2))./2;
    text(xpos(1)+1100,ymittl+10,['Stabw.: ',num2str(roundn(nanstd(diff_clim_rh),-1)),' %'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    text(xpos(1)+8000,ymittl+10,['T/R.: ',num2str(roundn((Tr_rh_max),-2)),' (',num2str(datestr(jan_15(xxx2),'yyyy')),'/',num2str(datestr(jan_15(xxx2+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    text(xpos(1)+1100,ymittl+8,['Klimatol. Mittel (1960-90): ',num2str(roundn(clim_ref_rh,0)),' %'],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
    if Tr_rh_max>1.64
    text(xpos(1)+15500,ymittl+10,['Linearer Trend: ',num2str(roundn(trend_rh,-1)),'%/ a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
    end
    set(gca,'XTick',jan_15(2:10:end))
    set(gca,'XTickLabel',jan_15(2:10:end))
    datetick('x','keeplimits','keepticks')
    set(gca,'fontsize',14)
    title(['Abweichung der saisonalen relativen Feuchte ',num2str(jahr_start),' - ',num2str(jahr_end),' (Okt.-Apr.) vom langjährigen Mittel (1960-90), ',num2str(location)],'fontsize',12,'fontweight','bold')
%     line([jan_15(47)-175 jan_15(47)-175],[min(diff_clim_rh) max(diff_clim_rh)+10],'linewidth',2,'color',[0.7 0.7 0.7])


saveas(m, 'fig2.jpg');
saveas(m, 'fig2.fig');