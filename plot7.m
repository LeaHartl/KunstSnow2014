function[]=plot7(tw2_dailyamp_mean, tw2_dailyamp_rel, tw2_dailyamp_prop, tw2_dailyamp_lan);
   
%    tw2_dailyamp_mean=tw2_dailyamp_mean(winter);
%    tw2_dailyamp_rel=tw2_dailyamp_rel(winter);
%    tw2_dailyamp_prop=tw2_dailyamp_prop(winter);
%    tw2_dailyamp_lan=tw2_dailyamp_lan(winter);
    
%Plot mittlerer Tagesgang Okt-Apr
 m =  figure('name','Mittlerer Tagesgang Saison')
    subplot(2,2,1)
    set(gca,'fontsize',14)
    plot(0:23,tw2_dailyamp_mean,'r')
    xlabel('Uhrzeit (UTC)')
    ylabel('Feuchttemperatur (°C)')
    grid on
    axis tight
    title('Mittlerer Tagesgang der Feuchttemperatur')
        
    subplot(2,2,2)
    set(gca,'fontsize',14)
    h=pie3(tw2_dailyamp_rel);
    colormap jet
    textObjs = findobj(h,'Type','text');
    oldStr = get(textObjs,{'String'})
    val = get(textObjs,{'Extent'});
    oldExt = cat(1,val{:});
 %  Names = {'0-5 Uhr: ';'6-11 Uhr: ';'18-23 Uhr: '};
 %  Names = {'0-5 Uhr: ';'18-23 Uhr: '};
    Names = {'0-5 Uhr: ';'6-11 Uhr: ';'12-17 Uhr: ';'18-23 Uhr: '};
    newStr = strcat(Names,oldStr);
    set(textObjs,{'String'},newStr);
    grid on 
    title('Zeitliche Verteilung der künstlich erzeugbaren Schneemenge')
    
    subplot(2,2,3)
    set(gca,'fontsize',14)
    bar(0:23,tw2_dailyamp_prop,'b')
    xlabel('Uhrzeit (UTC)')
    ylabel('Schneileistung (m³/h)')
    grid on
    axis tight
    title('Mittlerer Tagesgang der Schneileistung (Propellererzeuger)')
    
    subplot(2,2,4)
    set(gca,'fontsize',14)
    bar(0:23,tw2_dailyamp_lan,'k')
    xlabel('Uhrzeit (UTC)')
    ylabel('Schneileistung (m³/h)')
    grid on
    axis tight
    title('Mittlerer Tagesgang der Schneileistung (Lanzenerzeuger)')
%     
        saveas(m, 'zeitlVert.jpg');
saveas(m, 'zeitlVert.fig');