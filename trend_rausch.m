function[Tr_max, trend, xxx, y]=trend_rausch(data, jan_15); 
%Trend/Rauschverhältnis berechnen, linearen Trend berechnen
     
%signifikanten Trend finden; Trend ist dann signifikant, wenn Trend/Rauschverhältnis >1.64 (90% Confidence-Level)


%   -------------------------------------------------------   
%      falls kein Trend/Abweichung vom saisonalen Mittel besteht data aus
%      Nullen und restliche subfunction geht nicht. Output in diesem Fall
%      auf Null setzen. 
         if sum(data)==0;
         Tr_max=0;
         trend=0;
         xxx=0;
         y=zeros(64, 59);
         y(y==0)=NaN;
         else
   %------------------------------------------------------
%     size(data)
%     size(jan_15)
   Tr=zeros(length(data)-5,1);

   for i=1:length(data)-5  
       good=~(isnan(data(i:end)));
       
       % [p(i,:),S(i,:)] = polyfit(jan_15(i:end),data(i:end),1);
        abc=jan_15(i:end);
        cde=data(i:end);
         [p(i,:),S(i,:)] = polyfit(abc(good),cde(good),1);
        y(i:length(jan_15),i) = polyval(p(i,:),jan_15(i:end));
        Tr(i,:)=(nanmax(y(:,i))+abs(nanmin(y(:,i))))./nanstd(data(i:end));
   end
   y(y==0)=NaN;
  
   xxx=find(Tr==nanmax(Tr));
   Tr_max=nanmax(Tr);
   

   mx=nanmax(y(:,xxx));
   ab= abs(min(y(:,xxx)));
    %added mean function to fix weird problem with april data kitzbühel.
    %should not affect correct data... remove both means for original
    %script. 
   % trend=(mean(max(y(:,xxx))+mean(abs(min(y(:,xxx))))))./(length(data(xxx:end)));
    
      trend=(nanmax(y(:,xxx))+abs(nanmin(y(:,xxx))))./(length(data(xxx:end)));
         end
         
         
