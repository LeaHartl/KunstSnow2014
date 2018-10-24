function[]=plot4(location, datum2, jan_1, jan_15, snow_mass_hour, snow_mass_seas, snow_mass_hour_lan, snow_mass_seas_lan, label)   ;

 [cc, dd]=size(snow_mass_hour);
 snow_mass_hour=reshape(snow_mass_hour,1,  cc*dd );
 
  [cc1, dd1]=size(snow_mass_hour_lan);
 snow_mass_hour_lan=reshape(snow_mass_hour_lan,1,  cc1*dd1 );
 


%Plot snow production potential in (m³/h) for propeller guns
 m=   figure('name','snow production potential')
    subplot(2,1,1)
    set(gca,'fontsize',14)
    hold on
    plot(datum2,snow_mass_hour,'b')
    grid off 
    axis tight
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    %datetick('x','keeplimits','keepticks')
    ylabel('Produktionspotential (m³/h)')
    title(['Produktionspotential mit 1 Propellererzeuger (stündliche Summe), ',num2str(location)])
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    hold on
    plot(jan_15,snow_mass_seas,'.r','MarkerSize',25)
    
%     if Tr_maxProp>1.64
%     plot(jan_15,yProp(:,xxxProp),'-.k','linewidth',2)
%     end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('saisonale Summe (m³)')
    %axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_oct)) round(max(diff_clim_count_tw_oct))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     %text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_prop),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_maxProp,-2)),' (',num2str(datestr(jan_15(xxxProp),'yyyy')),'/',num2str(datestr(jan_15(xxxProp+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_maxProp>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trendProp,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
    
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    axis 'tight'
    %datetick('x','keeplimits','keepticks')
    
    title(['Produktionspotential mit 1 Propellererzeuger (saisonale Summe), ',num2str(location)])
    
    saveas(m, 'snow prod pot.jpg');
saveas(m, 'snow prod pot.fig');
    
    
     %Plot snow production potential in (m³/h) for lance gun
  m2=  figure('name','snow production potential')
    subplot(2,1,1)
    set(gca,'fontsize',14)
    hold on
    plot(datum2,snow_mass_hour_lan,'b')
    grid off 
    axis tight
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    %datetick('x','keeplimits','keepticks')
    ylabel('Produktionspotential (m³/h)')
    title(['Produktionspotential mit 1 Lanzenerzeuger (stündliche Summe), ',num2str(location)])
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    hold on
    plot(jan_15,snow_mass_seas_lan,'.r','MarkerSize',25)
    
%      if Tr_maxLan>1.64
%     plot(jan_15,yLan(:,xxxLan),'-.k','linewidth',2)
%     end
    axis tight
    grid off
    xlabel('Jahre')
    ylabel('saisonale Summe (m³)')
%     %axis([min(jan_15)-150 max(jan_15)+150 fix(min(diff_clim_count_tw_oct)) round(max(diff_clim_count_tw_oct))])
%     xpos=get(gca,'xlim');
%     xmittl=(xpos(1)+xpos(2))./2;
%     ypos=get(gca,'ylim');
%     ymittl=(ypos(1)+ypos(2))./2;
%     %text(xpos(1)+500,ymittl,['Stabw.: ',num2str(roundn(std(diff_clim_count_tw_oct),-1))],'fontweight','bold','BackgroundColor','none','Fontsize',12);  
%     text(xpos(1)+8000,ymittl,['T/R.: ',num2str(roundn(Tr_maxLan,-2)),' (',num2str(datestr(jan_15(xxxLan),'yyyy')),'/',num2str(datestr(jan_15(xxxLan+1),'yy')),' - ',num2str(datestr(jan_15(end-1),'yyyy')),'/',num2str(datestr(jan_15(end),'yy')),')'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     if Tr_maxLan>1.64
%     text(xpos(1)+15500,ymittl,['Linearer Trend: ',num2str(roundn(trendLan,-1)),' / a'],'fontweight','bold','BackgroundColor','none','Fontsize',12); 
%     end
    
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    axis 'tight'
    %datetick('x','keeplimits','keepticks')
    ylabel('saisonale Summe (m³)')
    title(['Produktionspotential mit 1 Lanzenerzeuger (saisonale Summe), ',num2str(location)])
    
    saveas(m2, 'snow prod pot erzeuger.jpg');
saveas(m2, 'snow prod pot erzeuger.fig');