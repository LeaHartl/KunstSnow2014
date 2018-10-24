function[tw2_dailyamp_2, tw2_dailyamp_prop, tw2_dailyamp_lan, tw2_dailyamp_rel, snow_mass_prop, snow_mass_lan, snow_mass_prop_sum, snow_mass_prop_sum_std, snow_mass_prop_sum_mean, snow_mass_lan_sum, snow_mass_lan_sum_std, snow_mass_lan_sum_mean]=tagesgang(tw_shortterm,q, n);

%q = number of days per month.

Q=reshape(tw_shortterm, length(tw_shortterm)/n, n);
tw2_dailyamp=nanmean(Q,2);
%tw2_dailyamp_oct=nanmean(tw_shortterm_oct,2);  %this line replaces the
%three above lines in orig. code. 



tw2_dailyamp_2=reshape(tw2_dailyamp,24,q);
tw2_dailyamp_2=nanmean(tw2_dailyamp_2,2);
    
%Umrechnung in Schneileistung
tw2_dailyamp_prop=-4.83.*tw2_dailyamp_2+3.94;
tw2_dailyamp_lan=-3.94.*tw2_dailyamp_2+4.23;
tw2_dailyamp_lan(tw2_dailyamp_lan<12.11)=0;
tw2_dailyamp_prop(tw2_dailyamp_prop<13.6)=0;
% tw2_dailyamp_lan(tw2_dailyamp_lan<0)=0;
% tw2_dailyamp_prop(tw2_dailyamp_prop<0)=0;
tw2_dailyamp_lan(tw2_dailyamp_lan>51)=51;
tw2_dailyamp_prop(tw2_dailyamp_prop>72)=72;


%Prozentueller Anteil Schneileistung in 6 Stunden Teilen
tw2_dailyamp_rel(1)=roundn((sum(tw2_dailyamp_prop(1:6))./sum(tw2_dailyamp_prop)).*100,0);
tw2_dailyamp_rel(2)=roundn((sum(tw2_dailyamp_prop(7:12))./sum(tw2_dailyamp_prop)).*100,0);
tw2_dailyamp_rel(3)=roundn((sum(tw2_dailyamp_prop(13:18))./sum(tw2_dailyamp_prop)).*100,0);
tw2_dailyamp_rel(4)=roundn((sum(tw2_dailyamp_prop(19:24))./sum(tw2_dailyamp_prop)).*100,0);


    %absolut oktober (kein mittlerer Tagesgang)
    snow_mass_prop=-4.83.*tw_shortterm+3.94;
%     snow_mass_prop(snow_mass_prop<13.6)=0;
    snow_mass_prop(snow_mass_prop>72)=72;
    snow_mass_prop(snow_mass_prop<0)=0;    
    snow_mass_prop_sum=nansum(snow_mass_prop);
    snow_mass_prop_sum_std=std(snow_mass_prop_sum);
    snow_mass_prop_sum_mean=mean(snow_mass_prop_sum);
    
    snow_mass_lan=-3.94.*tw_shortterm+4.23;
    %snow_mass_lan(snow_mass_lan<12.11)=0;
    snow_mass_lan(snow_mass_lan<0)=0;
    snow_mass_lan(snow_mass_lan>51)=51;
    snow_mass_lan_sum=nansum(snow_mass_lan);
    snow_mass_lan_sum_std=std(snow_mass_lan_sum);
    snow_mass_lan_sum_mean=mean(snow_mass_lan_sum);