function[season_mean]=seas_mean(seasNr, var);   

%Calculate seasonal mean (1st of Oct-30th of April). var is the placeholder
%for the variable given to function in main program (i.e. temp, rh...)

%     
%     season_mean(1)=nanmean(var(1:test2(1)-1));
%     for i=2:size(test2,1)
%         season_mean(i)=nanmean(var(test2(i-1):test2(i)-1));
%     end
%     
%     season_mean(i+1)=nanmean(var(test2(i):size(var,1)));
%     season_mean=season_mean';
%     
%     
%     %delete summer seasons   %ADJUST ACCORDINGLY!!
%     
%     season_mean=season_mean(1:2:end);
%     season_mean(end)=[];
% size(var)
% seasNr

size(var)

hoursperseason=212*24;

season_mean=reshape(var, hoursperseason, seasNr);

season_mean=nanmean(season_mean);