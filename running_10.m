function[data_10]=running_10(data);  
%Calculate 10-season running-mean
   
    windowSize = 10;
    data_10=filter(ones(1,windowSize)/windowSize,1,data);