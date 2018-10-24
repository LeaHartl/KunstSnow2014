function[tw_snow, tw_nosnow, tw_snow_2, tw_nosnow_2, tw_snow_3, tw_nosnow_3,tw_snow_4, tw_nosnow_4, snow_hours, snow_hours_seas]=prep_tw_plot(tw2, seasNr, tl2m);
    hoursperseason=212*24;
%Prepare color coded Tw plot

        %create 2 tw series with values >-1.5 and <-1.5

        tw2=tw2';
    tw_snow=tw2;
    tw_nosnow=tw2;

    tw_snow(tw_snow>=-1.5)=NaN;
    tw_nosnow(tw_nosnow<-1.5)=NaN;
    
        %Create additional tw-series with lower limits
        tw_snow_2=tw2;
        tw_nosnow_2=tw2;

        tw_snow_2(tw_snow_2>=-2)=NaN;
        tw_nosnow_2(tw_nosnow_2<-2)=NaN;
    
        %Create additional tw-series with lower limits
        tw_snow_3=tw2;
        tw_nosnow_3=tw2;

        tw_snow_3(tw_snow_3>=-3)=NaN;
        tw_nosnow_3(tw_nosnow_3<-3)=NaN;
        
        %Create additional tw-series with lower limits
        tw_snow_4=tw2;
        tw_nosnow_4=tw2;

        tw_snow_4(tw_snow_4>=-4)=NaN;
        tw_nosnow_4(tw_nosnow_4<-4)=NaN;
%---------------------------------------------------------------------------------     
         
    %Calculate total snowmaking hours
    qqq=isnan(tw_snow);
    qqq(qqq==1)=[];
    snow_hours=size(qqq,2);
    
    %Calculate seasonal snowmaking hours for Tw limit=-1.5°C
    
    tw_snow1=reshape(tw_snow, hoursperseason, seasNr);
    for i=1:seasNr;  
    Y=isnan(tw_snow1(:,i));
    snow_hours_seas(1,i)=length(find(Y==0));
    end
    
    clear Y
    
      tw_snow_21=reshape(tw_snow_2, hoursperseason, seasNr);
    for i=1:seasNr;  
    Y=isnan(tw_snow_21(:,i));
    snow_hours_seas(2,i)=length(find(Y==0));
    end
    
    clear Y
    
        tw_snow_31=reshape(tw_snow_3, hoursperseason, seasNr);
    for i=1:seasNr;  
       Y=isnan(tw_snow_31(:,i));
    snow_hours_seas(3,i)=length(find(Y==0));
    end
    
    clear Y
        tw_snow_41=reshape(tw_snow_4, hoursperseason, seasNr);
    for i=1:seasNr;  
      Y=isnan(tw_snow_41(:,i));
    snow_hours_seas(4,i)=length(find(Y==0));
    end
    
    clear Y
%     qqq2=isnan(tw_snow(1:test2(1)-1));
%     qqq2(qqq2==1)=[];
%     snow_hours_seas(1,1)=size(qqq2,2);
%     clear qqq2    
%     for i=2:size(test2,1)
%         clear qqq2
%         qqq2=isnan(tw_snow(test2(i-1):test2(i)-1));
%         qqq2(qqq2==1)=[];
%         snow_hours_seas(1,i)=size(qqq2,2);
%     end
%     clear qqq2
%     qqq2=isnan(tw_snow(test2(i):size(tl2m,1)));
%     qqq2(qqq2==1)=[];
%     snow_hours_seas(1,i+1)=size(qqq2,2);
    
    
    clear qqq2
    %Calculate seasonal snowmaking hours for Tw limit=-2°C
%     qqq2=isnan(tw_snow_2(1:test2(1)-1));
%     qqq2(qqq2==1)=[];
%     snow_hours_seas(2,1)=size(qqq2,2);
%     clear qqq2    
%     for i=2:size(test2,1)
%         clear qqq2
%         qqq2=isnan(tw_snow_2(test2(i-1):test2(i)-1));
%         qqq2(qqq2==1)=[];
%         snow_hours_seas(2,i)=size(qqq2,2);
%     end
%     clear qqq2
%     qqq2=isnan(tw_snow_2(test2(i):size(tl2m,1)));
%     qqq2(qqq2==1)=[];
%     snow_hours_seas(2,i+1)=size(qqq2,2);
%     clear qqq2
%     
%     %Calculate seasonal snowmaking hours for Tw limit=-3°C
%     qqq2=isnan(tw_snow_3(1:test2(1)-1));
%     qqq2(qqq2==1)=[];
%     snow_hours_seas(3,1)=size(qqq2,2);
%     clear qqq2    
%     for i=2:size(test2,1)
%         clear qqq2
%         qqq2=isnan(tw_snow_3(test2(i-1):test2(i)-1));
%         qqq2(qqq2==1)=[];
%         snow_hours_seas(3,i)=size(qqq2,2);
%     end
%     clear qqq2
%     qqq2=isnan(tw_snow_3(test2(i):size(tl2m,1)));
%     qqq2(qqq2==1)=[];
%     snow_hours_seas(3,i+1)=size(qqq2,2);
%     clear qqq2
%     
%     %Calculate seasonal snowmaking hours for Tw limit=-4°C
%     qqq2=isnan(tw_snow_4(1:test2(1)-1));
%     qqq2(qqq2==1)=[];
%     snow_hours_seas(4,1)=size(qqq2,2);
%     clear qqq2    
%     for i=2:size(test2,1)
%         clear qqq2
%         qqq2=isnan(tw_snow_4(test2(i-1):test2(i)-1));
%         qqq2(qqq2==1)=[];
%         snow_hours_seas(4,i)=size(qqq2,2);
%     end
%     clear qqq2
%     qqq2=isnan(tw_snow_4(test2(i):size(tl2m,1)));
%     qqq2(qqq2==1)=[];
%     snow_hours_seas(4,i+1)=size(qqq2,2);
%     clear qqq2
%     
%     
%     %delete summer seasons from snow_hours_seas  %ADJUST ACCORDINGLY!!
%     %size(snow_hours_seas)
%     snow_hours_seas=snow_hours_seas(:,1:2:end);
%     %size(snow_hours_seas)
%     snow_hours_seas(:,end)=[];
%     
%     