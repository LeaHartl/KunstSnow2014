function[]=plot2(jan_15, snow_hours_seas, jan_1, location, alti, label);   

%Plot seasonal snowmaking hours with different TW limits
m=    figure('name','snowmaking hours diff Tw limits')
    subplot(2,1,1)
    set(gca,'fontsize',14)
    bar(jan_15,snow_hours_seas(2:4,:)','group')
    grid on
    axis tight
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    %datetick('x','keeplimits','keepticks')
    ylabel('Beschneistunden')
    colormap summer
    legend('-2°C','-3°C','-4°C')
   title([num2str(location),', ',num2str(alti),' m'])
    
    subplot(2,1,2)
    set(gca,'fontsize',14)
    plot(jan_15,snow_hours_seas(2,:),'-s',jan_15,snow_hours_seas(3,:),'-*',jan_15,snow_hours_seas(4,:),'-o');
    lsline;
    grid on
    axis tight
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    %datetick('x','keeplimits','keepticks')
    ylabel('Beschneistunden')
    colormap summer
    legend('-2°C','-3°C','-4°C')
    title([num2str(location),', ',num2str(alti),' m'])
    
 saveas(m, 'snowmakinghours.jpg');
saveas(m, 'snowmakinghours.fig');   
% %     figure(2)
% %        plot(jan_15,season_mean_tw,'k','linewidth',2)
% %     legend('Keine Beschneiung','Beschneiung möglich','Mittelwert Saison (Okt-Apr)','location','best','orientation','horizontal')
% %     %legend('boxoff')
% %     h3 = axes('Position',get(h1,'Position'));
% %     set(gca,'fontsize',12)
% %     bar(jan_15,snow_hours_seas(1,:),0.05,'b'); 
% %     set(h3,'YAxisLocation','right','Color','none','XTickLabel',[])
% %     set(h3,'XLim',get(h1,'XLim'))
% %     set(h3,'XTick',jan_1)
% %     ylabel('Beschneistunden')
% %     legend('Beschneistunden','location','best')
% %     legend('boxoff')