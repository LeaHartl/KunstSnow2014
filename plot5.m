function[]=plot5(location, tw2, tw_dist_2, tw_dist_4);    
%PLOT kumulative Verteilungsfunktion
  m=  figure('name','CDFPLOT')
    set(gca,'fontsize',12)
    h1=cdfplot(tw2);
    hold on
    h2=cdfplot(tw_dist_2);
    h3=cdfplot(tw_dist_4);
    legend('Alle Werte','Stundenwerte von Tagen mit Tagesmittel <=-2°C','Stundenwerte mit Tagesmittel <= -2°C und > -3°C','location','best')
    set(h2,'color','r','linestyle','--')
    set(h3,'color','g','linestyle','-.')
   
    title(' ')
    
        saveas(m, 'cdf.jpg');
saveas(m, 'cdf.fig');