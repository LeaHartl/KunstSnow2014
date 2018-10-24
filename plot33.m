function[]=plot33(zz, location, zz_day, zz_2);    

%Time-dependance of snowmaking
 m=   figure('name','time dependance')
    
   %tägliche Variabilität
    subplot(2,2,1)
    set(gca,'fontsize',14)
    bar(0:23,zz)
    grid on
    axis tight
    xlabel('Uhrzeit (UTC)')
    ylabel('Anzahl Beschneistunden')
    title(num2str(location))
    axis([0 23 fix(min(zz))-100 round(max(zz))])     
    subplot(2,2,2)
    pie3(zz_day);
    title('Relative zeitliche Verteilung der möglichen Beschneistunden (9-17/17-8 Uhr)')
    
   %monatliche Variabilität 
    subplot(2,2,3)
    set(gca,'fontsize',14)
    bar(1:7,zz_2,0.3,'k')
    grid on
    axis tight
    xlabel('Monate')
    ylabel('Anzahl Beschneistunden')
    title(num2str(location))
    set(gca,'XTicklabel','10|11|12|1|2|3|4')
    
    subplot(2,2,4)
    set(gca,'fontsize',14)
    explode=[0 0 1 1 1 1 0];
    pie3(zz_2,explode);
    title('Relative zeitliche Verteilung der möglichen Beschneistunden (Monate)')
    colorbar('location','southoutside','Xticklabel',{'Okt','Nov','Dez','Jan','Feb','Mar','Apr'})
    
saveas(m, 'time dependance.jpg');
saveas(m, 'time dependance.fig');