function[]=trend_version2(location, alti, diff_clim_t, diff_clim_rh, diff_clim_tw,diff_clim_count_tw_2 , jan_15, diff_clim_count_tw_oct,diff_clim_count_tw_nov,...
    diff_clim_count_tw_dec,diff_clim_count_tw_jan,diff_clim_count_tw_feb,diff_clim_count_tw_mar,diff_clim_count_tw_apr );  


% central_year_sorted, window_with, window_with_sorted, sen_mk_sorted, sen_mk_sorted_2, sen_mk_sorted_3, sen_mk_sorted_4,...
%     sen_mk_sorted_oct, sen_mk_sorted_nov, sen_mk_sorted_dec , sen_mk_sorted_jan, sen_mk_sorted_feb, sen_mk_sorted_mar, sen_mk_sorted_apr,  min_trend_abs,  max_trend_abs, sig_mk_4_sorted




   %sen_mk_sorted_oct, sen_mk_sorted_nov, sen_mk_sorted_dec , sen_mk_sorted_jan, sen_mk_sorted_feb, sen_mk_sorted_mar, sen_mk_sorted_apr]=trend_rausch(diff_clim_t, diff_clim_rh, diff_clim_tw,diff_clim_count_tw_2 , jan_15);  
%, window_with_sorted, sen_mk_sorted, sen_mk_sorted_2,...
  %  sen_mk_sorted_3, sen_mk_sorted_4, sen_mk_sorted_oct, sen_mk_sorted_nov, sen_mk_sorted_dec, sen_mk_sorted_jan, sen_mk_sorted_feb, sen_mk_sorted_mar, sen_mk_sorted_apr]

%NEU: Mann-Kendall Trendtest (verteilungsfrei) für Publikation
   %Beschneiungsklimatol.
   
    %1-D Version with fixed Window end 2007
%     for i=1:size(diff_clim_t,1)-5 
%     x1(:,1)=jan_15(i:end);
%     x1(:,2)=diff_clim_t(i:end);
%     [taub_mk(i,:) h_mk(i,:) sig_mk(i,:) Z_mk(i,:) S_mk(i,:) sigma_mk(i,:) sen_mk(i,:) plotofslope_mk(i,:)] = ktaub(x1,0.05);
%     clear x1
%     end     
   
%Irrtumswahrscheinlichkeit p_stat für statistische Trendanalyse nach
%Mann-Kendall definieren.
p_stat=0.05;

    %2-D Version with variabel window with (Running trend analysis)
    k=0;
   % cde=1
size(diff_clim_t);
diff_clim_t;
    for i=1:length(diff_clim_t)-5;
        
      for j=10:length(diff_clim_t)-k;
         % cde=2

        %Air temp
        x1(:,1)=jan_15(i:j+k);
        x1(:,2)=diff_clim_t(i:j+k);
       diff_clim_t(i:10);
       [taub_mk(j-9,i) tau h_mk(j-9,i) sig_mk(j-9,i) Z_mk(j-9,i) S_mk(j-9,i) sigma_mk(j-9,i) sen_mk(j-9,i) n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(x1,p_stat, 0);
        %[taub_mk(j-9,i) tau_mk(j-9,i) h_mk(j-9,i) sig_mk(j-9,i) Z_mk(j-9,i) S_mk(j-9,i) sigma_mk(j-9,i) sen_mk(j-9,i) n_mk(j-9,i) plotofslope_mk(j-9,i)...
         %   CIlower(j-9,i) CIupper(j-9,i) D(j-9,i) Dall(j-9,i) C3(j-9,i) nsigma(j-9,i)] = ktaub(x1,p_stat, 1);
        
      %taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma
        [p33,S33] = polyfit(x1(:,1),x1(:,2),1);
        y33(1:size(x1,1),1) = polyval(p33,x1(:,1));
        trend_tl_lin(j-9,i)=(max(y33+abs(min(y33))))./(size(x1,1));
        if y33(end)<y33(1);
            trend_tl_lin(j-9,i)=trend_tl_lin(j-9,i)*(-1);
        end
        clear p33 S33 y33
        


         %RH
        x2(:,1)=jan_15(i:j+k);
        x2(:,2)=diff_clim_rh(i:j+k);
        [taub_mk_2(j-9,i) tau h_mk_2(j-9,i) sig_mk_2(j-9,i) Z_mk_2(j-9,i) S_mk_2(j-9,i) sigma_mk_2(j-9,i) sen_mk_2(j-9,i) n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(x2,p_stat);
        
        [p33,S33] = polyfit(x2(:,1),x2(:,2),1);
        y33(1:size(x2,1),1) = polyval(p33,x2(:,1));
        trend_rh_lin(j-9,i)=(max(y33+abs(min(y33))))./(size(x2,1));
        if y33(end)<y33(1)
            trend_rh_lin(j-9,i)=trend_rh_lin(j-9,i)*(-1);
        end
        clear p33 S33 y33
        
        %Wet-bulb temp
        x3(:,1)=jan_15(i:j+k);
        x3(:,2)=diff_clim_tw(i:j+k);
        [taub_mk_3(j-9,i) tau h_mk_3(j-9,i) sig_mk_3(j-9,i) Z_mk_3(j-9,i) S_mk_3(j-9,i) sigma_mk_3(j-9,i) sen_mk_3(j-9,i) n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(x3,p_stat);

        [p33,S33] = polyfit(x3(:,1),x3(:,2),1);
        y33(1:size(x3,1),1) = polyval(p33,x3(:,1));
        trend_tw_lin(j-9,i)=(max(y33+abs(min(y33))))./(size(x3,1));
        if y33(end)<y33(1);
            trend_tw_lin(j-9,i)=trend_tw_lin(j-9,i)*(-1);
        end
        clear p33 S33 y33
        %Beschneitage -2°C Tw
        x4(:,1)=jan_15(i:j+k);
        x4(:,2)=diff_clim_count_tw_2(i:j+k);
        [taub_mk_4(j-9,i) tau h_mk_4(j-9,i) sig_mk_4(j-9,i) Z_mk_4(j-9,i) S_mk_4(j-9,i) sigma_mk_4(j-9,i) sen_mk_4(j-9,i) n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(x4,p_stat);

        [p33,S33] = polyfit(x4(:,1),x4(:,2),1);
        y33(1:size(x4,1),1) = polyval(p33,x4(:,1));
        trend_tw_count_lin(j-9,i)=(max(y33+abs(min(y33))))./(size(x4,1));
        if y33(end)<y33(1);
            trend_tw_count_lin(j-9,i)=trend_tw_count_lin(j-9,i)*(-1);
        
        end
        clear p33 S33 y33
        

      central_year(j-9,i)=str2double(datestr(mean(x1(:,1)),'yyyy')); %(str2num(datestr(jan_15(j+k),'yyyy'))+str2num(datestr(jan_15(i),'yyyy')))./2;
      window_with(j-9,i)=j+k-i+1;
      
     
      clear x1
      clear x2
      clear x3
      clear x4
      end
         k=k+1;
    end

    abc=1
    central_year(central_year==0)=NaN;
    clear min_year
    min_year=min(min(central_year));
    max_year=max(max(central_year));
    
    central_year_sorted=NaN(size(central_year));
    window_with_sorted=NaN(size(window_with));
    h_mk_sorted=NaN(size(h_mk));
    sen_mk_sorted=NaN(size(sen_mk));
    
    
    h_mk_sorted_2=NaN(size(h_mk_2));
    sen_mk_sorted_2=NaN(size(sen_mk_2));
    h_mk_sorted_3=NaN(size(h_mk_3));
    sen_mk_sorted_3=NaN(size(sen_mk_3));
    trend_tw_lin_sorted=NaN(size(trend_tw_lin));
    trend_tl_lin_sorted=NaN(size(trend_tl_lin));
    trend_rh_lin_sorted=NaN(size(trend_rh_lin));
    h_mk_sorted_4=NaN(size(h_mk_4));
    sen_mk_sorted_4=NaN(size(sen_mk_4));
    trend_tw_count_lin_sorted=NaN(size(trend_tw_count_lin));
    taub_mk_3_sorted=NaN(size(taub_mk_3));
    sig_mk_3_sorted=NaN(size(sig_mk_3));
    sig_mk_4_sorted=NaN(size(sig_mk_4));
    
       for j=1:size(central_year,2)
           [m,n]=find(central_year==min_year);
               for k=1:size(m,1)
            central_year_sorted(m(k),j)=central_year(m(k),n(k));
            window_with_sorted(m(k),j)=window_with(m(k),n(k));
            h_mk_sorted(m(k),j)=h_mk(m(k),n(k));
            sen_mk_sorted(m(k),j)=sen_mk(m(k),n(k));
            
            h_mk_sorted_2(m(k),j)=h_mk_2(m(k),n(k));
            sen_mk_sorted_2(m(k),j)=sen_mk_2(m(k),n(k));
            h_mk_sorted_3(m(k),j)=h_mk_3(m(k),n(k));
            sen_mk_sorted_3(m(k),j)=sen_mk_3(m(k),n(k));
            trend_tw_lin_sorted(m(k),j)=trend_tw_lin(m(k),n(k));
            trend_tl_lin_sorted(m(k),j)=trend_tl_lin(m(k),n(k));
            trend_rh_lin_sorted(m(k),j)=trend_rh_lin(m(k),n(k));
            h_mk_sorted_4(m(k),j)=h_mk_4(m(k),n(k));
            sen_mk_sorted_4(m(k),j)=sen_mk_4(m(k),n(k));
            trend_tw_count_lin_sorted(m(k),j)=trend_tw_count_lin(m(k),n(k));
            taub_mk_3_sorted(m(k),j)=taub_mk_3(m(k),n(k));
            sig_mk_3_sorted(m(k),j)=sig_mk_3(m(k),n(k));
            sig_mk_4_sorted(m(k),j)=sig_mk_4(m(k),n(k));
               end
                    if min_year<max_year
                        min_year=min_year+1;
                    end
       end
    
 sen_mk_sorted(h_mk_sorted==0)=NaN;
 sen_mk_sorted(sen_mk_sorted==0)=NaN;
 sen_mk_sorted_2(h_mk_sorted_2==0)=NaN;
 sen_mk_sorted_2(sen_mk_sorted_2==0)=NaN;
 sen_mk_sorted_3(h_mk_sorted_3==0)=NaN; 
 sen_mk_sorted_3(sen_mk_sorted_3==0)=NaN;
 trend_tw_lin_sorted(h_mk_sorted_3==0)=NaN;
 trend_tl_lin_sorted(h_mk_sorted==0)=NaN;
 trend_rh_lin_sorted(h_mk_sorted_2==0)=NaN;
 sen_mk_sorted_4(h_mk_sorted_4==0)=NaN;
 sen_mk_sorted_4(sen_mk_sorted_4==0)=NaN;
 trend_tw_count_lin_sorted(h_mk_sorted_4==0)=NaN;
 trend_tw_lin_sorted(trend_tw_lin_sorted==0)=NaN;
 trend_tl_lin_sorted(trend_tl_lin_sorted==0)=NaN;
 trend_rh_lin_sorted(trend_rh_lin_sorted==0)=NaN;
 trend_tw_count_lin_sorted(trend_tw_count_lin_sorted==0)=NaN;
 sig_mk_3_sorted(h_mk_sorted_3==0)=NaN;
 sig_mk_4_sorted(h_mk_sorted_4==0)=NaN;

 

 


 %--------------------------------
 
    k=0;
    clear p33 S33 y33
    for i=1:length(diff_clim_t)-5
      for j=10:length(diff_clim_t)-k
        
        %Monthly snowmaking days
        
        %oct
        x_oct(:,1)=jan_15(i:j+k);
        x_oct(:,2)=diff_clim_count_tw_oct(i:j+k);
        [taub_mk_oct(j-9,i) tau h_mk_oct(j-9,i) sig_mk_oct(j-9,i) Z_mk_oct(j-9,i) S_mk_oct(j-9,i) sigma_mk_oct(j-9,i) sen_mk_oct(j-9,i) n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(x_oct,p_stat);
        
        [p33,S33] = polyfit(x_oct(:,1),x_oct(:,2),1);
        y33(1:size(x_oct,1),1) = polyval(p33,x_oct(:,1));
        trend_tw_oct_lin(j-9,i)=(max(y33+abs(min(y33))))./(size(x_oct,1));
        if y33(end)<y33(1)
            trend_tw_oct_lin(j-9,i)=trend_tw_oct_lin(j-9,i)*(-1);
        end
        clear p33 S33 y33
        
      clear x_oct
      end
         k=k+1;
    end
   
    h_mk_sorted_oct=NaN(size(h_mk_oct));
    sen_mk_sorted_oct=NaN(size(sen_mk_oct));
    trend_tw_oct_lin_sorted=NaN(size(trend_tw_oct_lin));
    sig_mk_sorted_oct=NaN(size(sig_mk_oct));
    
     min_year=min(min(central_year));
    
     for j=1:size(central_year,2)
           [m,n]=find(central_year==min_year);
               for k=1:size(m,1)
            h_mk_sorted_oct(m(k),j)=h_mk_oct(m(k),n(k));
            sen_mk_sorted_oct(m(k),j)=sen_mk_oct(m(k),n(k));
            trend_tw_oct_lin_sorted(m(k),j)=trend_tw_oct_lin(m(k),n(k));
            sig_mk_sorted_oct(m(k),j)=sig_mk_oct(m(k),n(k));
               end
                    if min_year<max_year
                        min_year=min_year+1;
                    end
       end
    
 trend_tw_oct_lin_sorted(trend_tw_oct_lin_sorted==0)=NaN;
 sen_mk_sorted_oct(h_mk_sorted_oct==0)=NaN;
 sen_mk_sorted_oct(sen_mk_sorted_oct==0)=NaN;
 trend_tw_oct_lin_sorted(h_mk_sorted_oct==0)=NaN;
 sig_mk_sorted_oct(h_mk_sorted_oct==0)=NaN;
 
 % sen_mk_sorted_oct
    %nov
    
     k=0;
    for i=1:length(diff_clim_t)-5
      for j=10:length(diff_clim_t)-k
        x_nov(:,1)=jan_15(i:j+k);
        x_nov(:,2)=diff_clim_count_tw_nov(i:j+k);
        [taub_mk_nov(j-9,i) tau h_mk_nov(j-9,i) sig_mk_nov(j-9,i) Z_mk_nov(j-9,i) S_mk_nov(j-9,i) sigma_mk_nov(j-9,i) sen_mk_nov(j-9,i) n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(x_nov,p_stat);
        
        [p44,S44] = polyfit(x_nov(:,1),x_nov(:,2),1);
        y44(1:size(x_nov,1),1) = polyval(p44,x_nov(:,1));
        trend_tw_nov_lin(j-9,i)=(max(y44+abs(min(y44))))./(size(x_nov,1));
        if y44(end)<y44(1)
            trend_tw_nov_lin(j-9,i)=trend_tw_nov_lin(j-9,i)*(-1);
        end
        clear p44 S44 y44
        
      clear x_nov
      end
         k=k+1;
    end
   
    h_mk_sorted_nov=NaN(size(h_mk_nov));
    sen_mk_sorted_nov=NaN(size(sen_mk_nov));
    trend_tw_nov_lin_sorted=NaN(size(trend_tw_nov_lin));
    sig_mk_sorted_nov=NaN(size(sig_mk_nov));
    
     min_year=min(min(central_year));
    
     for j=1:size(central_year,2)
           [m,n]=find(central_year==min_year);
               for k=1:size(m,1)
            h_mk_sorted_nov(m(k),j)=h_mk_nov(m(k),n(k));
            sen_mk_sorted_nov(m(k),j)=sen_mk_nov(m(k),n(k));
            trend_tw_nov_lin_sorted(m(k),j)=trend_tw_nov_lin(m(k),n(k));
            sig_mk_sorted_nov(m(k),j)=sig_mk_nov(m(k),n(k));
               end
                    if min_year<max_year
                        min_year=min_year+1;
                    end
       end
    
 trend_tw_nov_lin_sorted(trend_tw_nov_lin_sorted==0)=NaN;
 sen_mk_sorted_nov(h_mk_sorted_nov==0)=NaN;
 sen_mk_sorted_nov(sen_mk_sorted_nov==0)=NaN;
 trend_tw_nov_lin_sorted(h_mk_sorted_nov==0)=NaN;
 sig_mk_sorted_nov(h_mk_sorted_nov==0)=NaN;
 
 %dec
    
     k=0;
    for i=1:length(diff_clim_t)-5
      for j=10:length(diff_clim_t)-k
        x_dec(:,1)=jan_15(i:j+k);
        x_dec(:,2)=diff_clim_count_tw_dec(i:j+k);
         [taub_mk_dec(j-9,i) tau h_mk_dec(j-9,i) sig_mk_dec(j-9,i) Z_mk_dec(j-9,i) S_mk_dec(j-9,i) sigma_mk_dec(j-9,i) sen_mk_dec(j-9,i) n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(x_dec,p_stat);
        
        [p44,S44] = polyfit(x_dec(:,1),x_dec(:,2),1);
        y44(1:size(x_dec,1),1) = polyval(p44,x_dec(:,1));
        trend_tw_dec_lin(j-9,i)=(max(y44+abs(min(y44))))./(size(x_dec,1));
        if y44(end)<y44(1)
            trend_tw_dec_lin(j-9,i)=trend_tw_dec_lin(j-9,i)*(-1);
        end
        clear p44 S44 y44
        
      clear x_dec
      end
         k=k+1;
    end
   
    h_mk_sorted_dec=NaN(size(h_mk_dec));
    sen_mk_sorted_dec=NaN(size(sen_mk_dec));
    trend_tw_dec_lin_sorted=NaN(size(trend_tw_dec_lin));
    sig_mk_sorted_dec=NaN(size(sig_mk_dec));
    
     min_year=min(min(central_year));
    
     for j=1:size(central_year,2)
           [m,n]=find(central_year==min_year);
               for k=1:size(m,1)
            h_mk_sorted_dec(m(k),j)=h_mk_dec(m(k),n(k));
            sen_mk_sorted_dec(m(k),j)=sen_mk_dec(m(k),n(k));
            trend_tw_dec_lin_sorted(m(k),j)=trend_tw_dec_lin(m(k),n(k));
            sig_mk_sorted_dec(m(k),j)=sig_mk_dec(m(k),n(k));
               end
                    if min_year<max_year
                        min_year=min_year+1;
                    end
       end
    
 trend_tw_dec_lin_sorted(trend_tw_dec_lin_sorted==0)=NaN;
 sen_mk_sorted_dec(h_mk_sorted_dec==0)=NaN;
 sen_mk_sorted_dec(sen_mk_sorted_dec==0)=NaN;
 trend_tw_dec_lin_sorted(h_mk_sorted_dec==0)=NaN;
 sig_mk_sorted_dec(h_mk_sorted_dec==0)=NaN;
 

 
 %Find max and min values of trend nb. of snowmaking days for all months to
 %unify colormaps
 min_trend(1)=min(min(sen_mk_sorted_oct));
 min_trend(2)=min(min(sen_mk_sorted_nov));
 min_trend(3)=min(min(sen_mk_sorted_dec));

 
 min_trend_abs=min(min_trend);
 
 [max_trend(1), I(1)]=max(max(sen_mk_sorted_oct));
 [max_trend(2), I(2)]=max(max(sen_mk_sorted_nov));
 [max_trend(3), I(3)]=max(max(sen_mk_sorted_dec));

 max_trend_abs=max(max_trend);
 
% k=cumsum(days_feb);  

% 
% [max_tr(1),ind(1)] = find(sen_mk_sorted_oct==max(sen_mk_sorted_oct(:)));
% [max_tr(2),ind(2)] = find(sen_mk_sorted_nov==max(sen_mk_sorted_nov(:)));
% [max_tr(3),ind(3)] = find(sen_mk_sorted_dec==max(sen_mk_sorted_dec(:)));



[xx_oct, yy_oct]= max(abs(sen_mk_sorted_oct(:)));

data(3,1)=xx_oct;

data(3,2)=window_with_sorted(yy_oct);

data(3,3)=central_year_sorted(yy_oct);

data(3,4)=sig_mk_sorted_oct(yy_oct);


clear BB central_year_sorted_1 window_with1

BB=find(~isnan(sen_mk_sorted_oct));

window_with1=window_with(BB);

sen_mk_sorted_oct_1=sen_mk_sorted_oct(BB);

central_year_sorted_1=central_year_sorted(BB);

sig_mk_oct_sorted_1=sig_mk_sorted_oct(BB);



[xx_oct2, yy_oct2]=max(window_with1(:));



if isempty(BB)
data(4,2)=NaN;
data(4,1)=NaN;
data(4,3)=NaN;
data(4,4)=NaN;
else
data(4,2)=xx_oct2;
data(4,1)=sen_mk_sorted_oct_1(yy_oct2);
data(4,3)=central_year_sorted_1(yy_oct2);
data(4,4)=sig_mk_oct_sorted_1(yy_oct2);
end

%----------------------November

[xx_nov, yy_nov]= max(abs(sen_mk_sorted_nov(:)));

data(5,1)=xx_nov;
data(5,2)=window_with_sorted(yy_nov);
data(5,3)=central_year_sorted(yy_nov);
data(5,4)=sig_mk_sorted_nov(yy_nov);


clear BB central_year_sorted_1 window_with1

BB=find(~isnan(sen_mk_sorted_nov));

window_with1=window_with(BB);
sen_mk_sorted_nov_1=sen_mk_sorted_nov(BB);
central_year_sorted_1=central_year_sorted(BB);
sig_mk_nov_sorted_1=sig_mk_sorted_nov(BB);

[xx_nov2, yy_nov2]=max(window_with1(:));


if isempty(BB)
data(6,2)=NaN;
data(6,1)=NaN;
data(6,3)=NaN;
data(6,4)=NaN;
else
data(6,2)=xx_nov2;
data(6,1)=sen_mk_sorted_nov_1(yy_nov2);
data(6,3)=central_year_sorted_1(yy_nov2);
data(6,4)=sig_mk_nov_sorted_1(yy_nov2);
end
%----------------------Dezember

[xx_dec, yy_dec]= max(abs(sen_mk_sorted_dec(:)));

data(7,1)=xx_dec;
data(7,2)=window_with_sorted(yy_dec);
data(7,3)=central_year_sorted(yy_dec);
data(7,4)=sig_mk_sorted_dec(yy_dec);


clear BB central_year_sorted_1 window_with1

BB=find(~isnan(sen_mk_sorted_dec));

window_with1=window_with(BB);
sen_mk_sorted_dec_1=sen_mk_sorted_dec(BB);
central_year_sorted_1=central_year_sorted(BB);
sig_mk_dec_sorted_1=sig_mk_sorted_dec(BB);

[xx_dec2, yy_dec2]=max(window_with1(:));

if isempty(BB)
data(8,2)=NaN;
data(8,1)=NaN;
data(8,3)=NaN;
data(8,4)=NaN;    
else
data(8,2)=xx_dec2;
data(8,1)=sen_mk_sorted_dec_1(yy_dec2);
data(8,3)=central_year_sorted_1(yy_dec2);
data(8,4)=sig_mk_dec_sorted_1(yy_dec2);
end





filename1=['table_',location,'.txt'];
dlmwrite(['mk_new\',filename1],data,'delimiter', ' ' , 'newline', 'pc');
% 
% %*-------------------------
% 
% 
%     
mk_oct=    figure(1)
   % title(['Running trend analysis with linear trends, significant after Mann-Kendall @ 95%, ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
    [h2,h]=contourf(central_year_sorted(1,:),window_with(:,1),sen_mk_sorted_oct.*365,50);
    set(h,'LineStyle','none');
    
    xlabel('Central years of window width','fontsize',24)
    ylabel('Window (years)','fontsize',24)
    %caxis([min_trend_abs.*365 max_trend_abs.*365])
    caxis([-2 2])
    colorbar('location','northoutside');%,'Outerposition',[0.071 0.96 0.433 0.045])
    grid on
    title([num2str(location),', ',num2str(alti),'m'],'fontsize',24)
    set(gca,'fontsize',24);
    %text(central_year_sorted(1,2),window_with(end-5,1),'Nb. snowmaking days Oct','fontweight','bold','backgroundcolor',[.9 .9 .9],'fontsize',9)
    hold on
    line([min(min(central_year_sorted)) nanmean(nanmean(central_year_sorted)) max(max(central_year_sorted))], [min(min(window_with_sorted))  max(max(window_with_sorted)) min(min(window_with_sorted))],'color','k')
    
load('MyColormaps','mycmap')
colormap(gca,mycmap)
saveas(mk_oct, ['mk_new\', location 'oct.jpg']);
saveas(mk_oct, ['mk_new\', location 'oct.fig']);
saveas(mk_oct, ['mk_new\', location 'oct.eps']);
    
mk_nov=    figure(2)
    [h2,h]=contourf(central_year_sorted(1,:),window_with(:,1),sen_mk_sorted_nov.*365,50);
    set(h,'LineStyle','none');

    xlabel('Central years of window width','fontsize',24)
    ylabel('Window (years)','fontsize',24)
    caxis([-2 2])
    colorbar('location','northoutside');%,'Outerposition',[0.514 0.96 0.433 0.045])
    grid on
    title([num2str(location),', ',num2str(alti),'m'],'fontsize',24)
        set(gca,'fontsize',24);
   % text(central_year_sorted(1,2),window_with(end-5,1),'Nb. snowmaking days Nov','fontweight','bold','backgroundcolor',[.9 .9 .9],'fontsize',9)
    hold on
    line([min(min(central_year_sorted)) nanmean(nanmean(central_year_sorted)) max(max(central_year_sorted))], [min(min(window_with_sorted))  max(max(window_with_sorted)) min(min(window_with_sorted))],'color','k')
  
load('MyColormaps','mycmap')
colormap(gca,mycmap)    
saveas(mk_nov, ['mk_new\', location 'nov.jpg']);
saveas(mk_nov, ['mk_new\', location 'nov.fig']);
saveas(mk_nov, ['mk_new\', location 'nov.eps']);
    
mk_dec=    figure(3)
    [h2,h]=contourf(central_year_sorted(1,:),window_with(:,1),sen_mk_sorted_dec.*365,50);
    set(h,'LineStyle','none');

    xlabel('Central years of window width','fontsize',24)
    ylabel('Window (years)','fontsize',24)
    caxis([-2 2])
    colorbar('location','northoutside');%,'Outerposition')%,[0.071 0.488 0.433 0.045])
    grid on
    title([num2str(location),', ',num2str(alti),'m'],'fontsize',24)
        set(gca,'fontsize',24);
    %text(central_year_sorted(1,2),window_with(end-5,1),'Nb. snowmaking days Dec','fontweight','bold','backgroundcolor',[.9 .9 .9],'fontsize',9)
    hold on
    line([min(min(central_year_sorted)) nanmean(nanmean(central_year_sorted)) max(max(central_year_sorted))], [min(min(window_with_sorted))  max(max(window_with_sorted)) min(min(window_with_sorted))],'color','k')

load('MyColormaps','mycmap')
colormap(gca,mycmap)
saveas(mk_dec, ['mk_new\', location 'dec.jpg']);
saveas(mk_dec, ['mk_new\', location 'dec.fig']);
saveas(mk_dec, ['mk_new\', location 'dec.eps']);
% 
% mk_4=    figure('name','p-Value Mann-Kendall Snowmaking days')
%     title(['p-Value Mann-Kendall snowmaking days, ',num2str(location),', ',num2str(alti),' m'],'fontsize',12,'fontweight','bold')
%     subplot(2,2,1)
%     [h2,h]=contourf(central_year_sorted(1,:),window_with(:,1),sig_mk_4_sorted,50);
%     set(h,'LineStyle','none');
%     xlabel('Central years of window with')
%     ylabel('Window (years)')
%     %caxis([min_trend_abs.*365 max_trend_abs.*365])
%     colorbar('location','northoutside');%,'Outerposition',[0.071 0.96 0.433 0.045])
%     grid on
%     title('p-Value Mann-Kendall snowmaking days (days/season)','fontsize',8)
%     text(central_year_sorted(1,2),window_with(end-5,1),'p-Value Mann-Kendall snowmaking days','fontweight','bold','backgroundcolor',[.9 .9 .9],'fontsize',9)
%     hold on
%     line([min(min(central_year_sorted)) nanmean(nanmean(central_year_sorted)) max(max(central_year_sorted))], [min(min(window_with_sorted))  max(max(window_with_sorted)) min(min(window_with_sorted))],'color','k')
% 
% saveas(mk_4, ['mk\', location 'figmk4.jpg']);
% saveas(mk_4, ['mk\', location 'figmk4.fig']);
% saveas(mk_4, ['mk\', location 'figmk4.eps']);
%   
%   
%   
%   
%   