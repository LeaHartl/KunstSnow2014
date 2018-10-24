    
    file1=fopen('zug_zukunft.txt');

C=textscan(file1, '%f %f', 'headerlines',0);
fclose(file1)

future_time=C{1};
tw_fut=C{2};
a=  size(tw_fut) 

[yy_f, mon_f, dd_f]=datevec(future_time);


BB=find(yy_f==2029 & mon_f==10 & dd_f==1)
BBB=find(yy_f==2049 & mon_f==4 & dd_f==30)

yy_f=yy_f(BB:BBB);
dd_f=dd_f(BB:BBB);
mon_f=mon_f(BB:BBB);
tw_fut=tw_fut(BB:BBB);



oct_f= find(mon_f==10);
nov_f= find(mon_f==11);
dec_f= find(mon_f==12);
jan_f= find(mon_f==1);
feb_f= find(mon_f==2);
mar_f= find(mon_f==3);
apr_f= find(mon_f==4);


% % Distribution of TW
    
    [foct_f, xioctf]=ksdensity(tw_fut(oct_f));
    tw_count_oct_f=size(find((tw_fut(oct_f))<=-2),1);
%-----------------nov-----------    
    [fnov_f, xinovf]=ksdensity(tw_fut(nov_f));
    tw_count_nov_f=size(find((tw_fut(nov_f))<=-2),1);  
%------------------------dec-----------
    [fdec_f, xidecf]=ksdensity(tw_fut(dec_f));
    tw_count_dec_f=size(find((tw_fut(dec_f))<=-2),1);
%---------------jan-----------------------    
    [fjan_f, xijanf]=ksdensity(tw_fut(jan_f));
    tw_count_jan_f=size(find((tw_fut(jan_f))<=-2),1); 
%---------------feb-----------------------        
    [ffeb_f, xifebf]=ksdensity(tw_fut(feb_f));
    tw_count_feb_f=size(find((tw_fut(feb_f))<=-2),1);
%-------------------mar----------------------
    [fmar_f, ximarf]=ksdensity(tw_fut(mar_f));
    tw_count_mar_f=size(find((tw_fut(mar_f))<=-2),1);
%---------------------apr-------------------
    [fapr_f, xiaprf]=ksdensity(tw_fut(apr_f));
    tw_count_apr_f=size(find((tw_fut(apr_f))<=-2),1);
    
    
figure('name','Distribution of Tw IPCC')
    set(gca,'fontsize',12)
plot( xioctf, foct_f, 'k-.')%  
    
%     subplot(4,2,1)
%     plot(xioct1,foct1, 'r', xioct2, foct2, 'b', xioct_f, foct_f, 'k-.')%  
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Okt.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,2)
%     plot(xinov1,fnov1, 'r', xinov2, fnov2, 'b', xinov_f, fnov_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Nov.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,3)
%     plot(xidec1,fdec1, 'r', xidec2, fdec2, 'b', xidec_f, fdec_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Dez.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,4)
%     plot(xijan1,fjan1, 'r', xijan2, fjan2, 'b', xijan_f, fjan_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Jan.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,5)
%     plot(xifeb1,ffeb1, 'r', xifeb2, ffeb2, 'b', xifeb3, ffeb_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Feb.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,6)
%     plot(ximar1,fmar1, 'r', ximar2, fmar2, 'b', ximar3, fmar_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Mar.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     subplot(4,2,7)
%     plot(xiapr1,fapr1, 'r', xiapr2, fapr2, 'b', xiapr3, fapr_f, 'k-.')%
%     axis([-30 20 0 0.2])
%     hline =line([-2, -2], ylim, 'Color', 'k');
%     title('Apr.');
%     %xlabel('Feuchttemperatur (°C)');
%     %ylabel('relative Häufigkeit');
%     
%     legend('1974-94', '1994-2014', '2030-2050');  
