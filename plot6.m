function[]=plot6( tw_shortterm_year, jahr_start);   
%Plot Beschneistunden (mittlere Anzahl, max min, täglich)


AA= datenum('2006-10-01', 'yyyy-mm-dd');
AAA=datenum('2007-04-30', 'yyyy-mm-dd');
aaa=[AA:1:AAA]

dateplot=aaa;
size(dateplot);
size(tw_shortterm_year(:,1))


   m= figure('name','mittlere Anzahl Beschneistunden')
    set(gca,'fontsize',14)
    h=area(dateplot,tw_shortterm_year(:,1));
    set(h,'FaceColor',[.5 .5 .5])
    axis tight
    grid on
%     set(gca,'XTick', dateplot)
%     set(gca,'XTickLabel', dateplot)
    set(gca,'XTick',[datenum(jahr_start,10,1,0,0,0) datenum(jahr_start,11,1,0,0,0) datenum(jahr_start,12,1,0,0,0) datenum(jahr_start+1,1,1,0,0,0) datenum(jahr_start+1,2,1,0,0,0) datenum(jahr_start+1,3,1,0,0,0) datenum(jahr_start+1,4,1,0,0,0)]);
    datetick('x','mmm','keeplimits')%,'keepticks')
    hold on
    plot(dateplot,tw_shortterm_year(:,2),'r',dateplot,tw_shortterm_year(:,3),'b')
    legend('Mittel','Max','Min')
    ylabel('Potentielle Beschneistunden pro Tag')
    
        saveas(m, 'Mittlere Anzahl Beschneistunden.jpg');
saveas(m, 'Mittlere Anzahl Beschneistunden.fig');