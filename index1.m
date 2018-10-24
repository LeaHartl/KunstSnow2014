function [] = index1(location, yyyy, mon, dd, tw_mean_d, t_mean_d,rel_mean_d, press, tw_longterm_oct, tw_longterm_nov, tw_longterm_dec, tw_longterm_jan, tw_longterm_feb, tw_longterm_mar, tw_longterm_apr, U)

AA=find(yyyy==1974 & mon==10 & dd==1);
AAA=find(yyyy==1994 & mon==4 & dd==30);

yyyy1=yyyy(AA:AAA);
dd1=dd(AA:AAA);
mon1=mon(AA:AAA);
tw_mean_d1=tw_mean_d(AA:AAA);
t_mean_d1=t_mean_d(AA:AAA);
rel_mean_d1=rel_mean_d(AA:AAA);

oct1= find(mon1==10);
nov1= find(mon1==11);
dec1= find(mon1==12);
jan1= find(mon1==1);
feb1= find(mon1==2);
mar1= find(mon1==3);
apr1= find(mon1==4);

AA2=find(yyyy==1994 & mon==10 & dd==1);
AAA2=find(yyyy==2014 & mon==4 & dd==30);

yyyy2=yyyy(AA2:AAA2);
dd2=dd(AA2:AAA2);
mon2=mon(AA2:AAA2);
tw_mean_d2=tw_mean_d(AA2:AAA2);
t_mean_d2=t_mean_d(AA2:AAA2);
rel_mean_d2=rel_mean_d(AA2:AAA2);

oct2= find(mon2==10);
nov2= find(mon2==11);
dec2= find(mon2==12);
jan2= find(mon2==1);
feb2= find(mon2==2);
mar2= find(mon2==3);
apr2= find(mon2==4);

ano_oct=nanmean(tw_mean_d2(oct2))-tw_mean_d2(oct2);
ano_nov=nanmean(tw_mean_d2(nov2))-tw_mean_d2(nov2);
ano_dec=nanmean(tw_mean_d2(dec2))-tw_mean_d2(dec2);
ano_jan=nanmean(tw_mean_d2(jan2))-tw_mean_d2(jan2);
ano_feb=nanmean(tw_mean_d2(feb2))-tw_mean_d2(feb2);
ano_mar=nanmean(tw_mean_d2(mar2))-tw_mean_d2(mar2);
ano_apr=nanmean(tw_mean_d2(apr2))-tw_mean_d2(apr2);


% % Distribution of TW

    [foct1, xioct1]=ksdensity(tw_mean_d1(oct1));
    tw_count_oct1=size(find(tw_mean_d1(oct1)<=-2),1);
    
    mxoct1=max(foct1);
    
    XX2=find(xioct1<=-2);
    XX2oct1=foct1(XX2(end));
    
    XX3=find(xioct1<=-3);

    XX3oct1=foct1(XX3(end));

    ind_oct1=(mxoct1-XX2oct1)/(XX2oct1-XX3oct1); 
    
    [foct2, xioct2]=ksdensity(tw_mean_d2(oct2));
    tw_count_oct2=size(find(tw_mean_d2(oct2)<=-2),1);
    
    mxoct2=max(foct2);
    
    XX2=find(xioct2<=-2);
    XX2oct2=foct2(XX2(end));
    
    XX3=find(xioct2<=-3);
    XX3oct2=foct2(XX3(end));
    
    ind_oct2=(mxoct2-XX2oct2)/(XX2oct2-XX3oct2);

%-----------------nov-----------    
    [fnov1, xinov1]=ksdensity(tw_mean_d1(nov1));
    tw_count_nov1=size(find(tw_mean_d1(nov1)<=-2),1);
    
        mxnov1=max(fnov1);
    
    XX2=find(xinov1<=-2);
    XX2nov1=fnov1(XX2(end));
    
    XX3=find(xinov1<=-3);
    XX3nov1=fnov1(XX3(end));
    
    ind_nov1=(mxnov1-XX2nov1)/(XX2nov1-XX3nov1); 
    
    [fnov2, xinov2]=ksdensity(tw_mean_d2(nov2));
    tw_count_nov2=size(find(tw_mean_d2(nov2)<=-2),1);
    
            mxnov2=max(fnov2);
    
    XX2=find(xinov2<=-2);
    XX2nov2=fnov2(XX2(end));
    
    XX3=find(xinov2<=-3);
    XX3nov2=fnov2(XX3(end));
    
    ind_nov2=(mxnov2-XX2nov2)/(XX2nov2-XX3nov2); 


%------------------------dec-----------
    
    [fdec1, xidec1]=ksdensity(tw_mean_d1(dec1));
    tw_count_dec1=size(find(tw_mean_d1(dec1)<=-2),1);
    
        mxdec1=max(fdec1);
    
    XX2=find(xidec1<=-2);
    XX2dec1=fdec1(XX2(end));
    
    XX3=find(xidec1<=-3);
    XX3dec1=fdec1(XX3(end));
    
    ind_dec1=(mxdec1-XX2dec1)/(XX2dec1-XX3dec1); 
    
    
    [fdec2, xidec2]=ksdensity(tw_mean_d2(dec2));
    tw_count_dec2=size(find(tw_mean_d2(dec2)<=-2),1);
    

            mxdec2=max(fdec2);
    
    XX2=find(xidec2<=-2);
    XX2dec2=fdec2(XX2(end));
    
    XX3=find(xidec2<=-3);
    XX3dec2=fdec2(XX3(end));
    
    ind_dec2=(mxdec2-XX2dec2)/(XX2dec2-XX3dec2); 
            
%---------------jan-----------------------    
    
    
    [fjan1, xijan1]=ksdensity(tw_mean_d1(jan1));
    tw_count_jan1=size(find(tw_mean_d1(jan1)<=-2),1);
    
        mxjan1=max(fjan1);
    
    XX2=find(xijan1<=-2);
    XX2jan1=fjan1(XX2(end));
    
    XX3=find(xijan1<=-3);
    XX3jan1=fjan1(XX3(end));
    
    ind_jan1=(mxjan1-XX2jan1)/(XX2jan1-XX3jan1); 
    
    
    
    [fjan2, xijan2]=ksdensity(tw_mean_d2(jan2));
    tw_count_jan2=size(find(tw_mean_d2(jan2)<=-2),1);
    
    
            mxjan2=max(fjan2);
    
    XX2=find(xijan2<=-2);
    XX2jan2=fjan2(XX2(end));
    
    XX3=find(xijan2<=-3);
    XX3jan2=fjan2(XX3(end));
    
    ind_jan2=(mxjan2-XX2jan2)/(XX2jan2-XX3jan2); 

%---------------feb-----------------------        
    
    [ffeb1, xifeb1]=ksdensity(tw_mean_d1(feb1));
    tw_count_feb1=size(find(tw_mean_d1(feb1)<=-2),1);
    
        mxfeb1=max(ffeb1);
    
    XX2=find(xifeb1<=-2);
    XX2feb1=ffeb1(XX2(end));
    
    XX3=find(xifeb1<=-3);
    XX3feb1=ffeb1(XX3(end));
    
    ind_feb1=(mxfeb1-XX2feb1)/(XX2feb1-XX3feb1); 
    
    

    [ffeb2, xifeb2]=ksdensity(tw_mean_d2(feb2));
    tw_count_feb2=size(find(tw_mean_d2(feb2)<=-2),1);
    
            mxfeb2=max(ffeb2);
    
    XX2=find(xifeb2<=-2);
    XX2feb2=ffeb2(XX2(end));
    
    XX3=find(xifeb2<=-3);
    XX3feb2=ffeb2(XX3(end));
    
    ind_feb2=(mxfeb2-XX2feb2)/(XX2feb2-XX3feb2); 
    
%-------------------mar----------------------
    
    [fmar1, ximar1]=ksdensity(tw_mean_d1(mar1));
    tw_count_mar1=size(find(tw_mean_d1(mar1)<=-2),1);
    
        mxmar1=max(fmar1);
    
    XX2=find(ximar1<=-2);
    XX2mar1=fmar1(XX2(end));
    
    XX3=find(ximar1<=-3);
    XX3mar1=fmar1(XX3(end));
    
    ind_mar1=(mxmar1-XX2mar1)/(XX2mar1-XX3mar1); 
    

    [fmar2, ximar2]=ksdensity(tw_mean_d2(mar2));
    tw_count_mar2=size(find(tw_mean_d2(mar2)<=-2),1);
    
        mxmar2=max(fmar2);
    
    XX2=find(ximar2<=-2);
    XX2mar2=fmar2(XX2(end));
    
    XX3=find(ximar2<=-3);
    XX3mar2=fmar2(XX3(end));
    
    ind_mar2=(mxmar2-XX2mar2)/(XX2mar2-XX3mar2);

%---------------------apr-------------------
    
    [fapr1, xiapr1]=ksdensity(tw_mean_d1(apr1));
    tw_count_apr1=size(find(tw_mean_d1(apr1)<=-2),1);
    
        mxapr1=max(fapr1);
    
    XX2=find(xiapr1<=-2);
    XX2apr1=fapr1(XX2(end));
    
    XX3=find(xiapr1<=-3);
    XX3apr1=fapr1(XX3(end));
    
    ind_apr1=(mxapr1-XX2apr1)/(XX2apr1-XX3apr1); 
    
    
    [fapr2, xiapr2]=ksdensity(tw_mean_d2(apr2));
    tw_count_apr2=size(find(tw_mean_d2(apr2)<=-2),1);
    
            mxapr2=max(fapr2);
    
    XX2=find(xiapr2<=-2);
    XX2apr2=fapr2(XX2(end));
    
    XX3=find(xiapr2<=-3);
    XX3apr2=fapr2(XX3(end));
    
    ind_apr2=(mxapr2-XX2apr1)/(XX2apr2-XX3apr2); 
    
%------------------------------------------------------------------

data_all(1,1)=ind_oct1;
data_all(2,1)=ind_oct2;

data_all(1,2)=ind_nov1;
data_all(2,2)=ind_nov2;

data_all(1,3)=ind_dec1;
data_all(2,3)=ind_dec2;

data_all(1,4)=ind_jan1;
data_all(2,4)=ind_jan2;

data_all(1,5)=ind_feb1;
data_all(2,5)=ind_feb2;

data_all(1,6)=ind_mar1;
data_all(2,6)=ind_mar2;

data_all(1,7)=ind_apr1;
data_all(2,7)=ind_apr2;


filename2=['tabelle.txt'];
dlmwrite(filename2,data_all,'delimiter', ' ' , 'newline', 'pc');
    
    
  