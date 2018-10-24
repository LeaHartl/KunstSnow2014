function[]=plot1(datum2, tl2m, rh, tw_nosnow, tw_snow, jan_1, label);  


 

% label=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 1960 ...
%      NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 1970, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 1980, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 1990, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 2010 ...
%      NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 2010, NaN, NaN, NaN, NaN];

m=figure('name','Temperature Plot');
    subplot(3,1,1)
    set(gca,'fontsize',14)
    plot(datum2,tl2m,'r')
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
%    datetick('x','keeplimits','keepticks')
    grid on
    axis tight
    ylabel('2 m Lufttemperatur(°C)')
    
     
    subplot(3,1,2)
    set(gca,'fontsize',14)
    title('RH Plot')
    plot(datum2,rh,'g')
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    %datetick('x','keeplimits','keepticks')
    grid on
    axis tight
    ylabel('Relative Feuchte (%)')
    
    
     %Plot color-coded Tw Plot
    
    subplot(3,1,3)
    set(gca,'fontsize',14)
    hold on
    h=stem(datum2,tw_nosnow,'r');
    grid on 
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    %datetick('x','keeplimits','keepticks')
    xlabel('Time (CET)')
    ylabel('Feuchttemperatur (°C)')
    set(h,'BaseValue',-1.5,'Marker','none')
    h1=gca;
    h2=stem(datum2,tw_snow,'g');
    set(h2,'BaseValue',-1.5,'Marker','none')
    set(gca,'XTick',jan_1)
    set(gca,'XTickLabel',label)
    %datetick('x','keeplimits','keepticks')
    axis tight
    
    
saveas(m, 'Temperatur Plot.jpg');
saveas(m, 'Temperatur Plot.fig');