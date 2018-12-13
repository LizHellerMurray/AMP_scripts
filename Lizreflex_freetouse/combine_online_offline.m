%combine online (analysis) and offline (analysis_off)
%%%%% prior to this script load in both analysis files, by saving the
%%%%% offline one as "analysis_off"

 close all 

%avg
field4 =analysis.neg.ST.Full ;
field4a  = analysis_off.neg.ST.Full;
combine.fullneg = [field4, field4a];

combine.fullneg(combine.fullneg<-5) = nan;
combine.fullneg(combine.fullneg>5) = nan;

if ~isempty(combine.fullneg)
combine.negAvg_excludeOver5ST = nanmean(combine.fullneg,2);
end
%%%%%%

field12 =analysis.pos.ST.Full ;
field12a  = analysis_off.pos.ST.Full;
combine.fullpos = [field12, field12a];

combine.fullpos(combine.fullpos<-5) = nan;
combine.fullpos(combine.fullpos>5) = nan;

if ~isempty(combine.fullpos)
combine.posAvg_excludeOver5ST = nanmean(combine.fullpos,2);
end




field1 = analysis.graphNegKeepALL;
field1a = analysis_off.graphNegKeepALL;
combine.graphNegKeepALL = [field1, field1a];

field2 =analysis.graphNegexcludeALL;
field2a =analysis_off.graphNegexcludeALL;
combine.graphNegexcludeALL = [field2, field2a];

field3 = analysis.graphNegKeep_excludeOver5ST ;
field3a = analysis_off.graphNegKeep_excludeOver5ST;
combine.graphNegKeep_excludeOver5ST = [field3, field3a];


field5 = analysis.graphNegKeepFollow_excludeOver5ST;
field5a =analysis_off.graphNegKeepFollow_excludeOver5ST;
combine.graphNegKeepFollow_excludeOver5ST = [field5, field5a];

field6 = analysis.graphNegexcludeFollow_excludeOver5ST;
field6a =  analysis_off.graphNegexcludeFollow_excludeOver5ST;
combine.graphNegexcludeFollow_excludeOver5ST = [field6, field6a];

field7 = analysis.graphNegKeepOpp_excludeOver5ST;
field7a = analysis_off.graphNegKeepOpp_excludeOver5ST;
combine.graphNegKeepOpp_excludeOver5ST = [field7, field7a];

field8  = analysis.graphNegexcludeOpp_excludeOver5ST; 
field8a = analysis_off.graphNegexcludeOpp_excludeOver5ST ;
combine.graphNegexcludeOpp_excludeOver5ST = [field8, field8a];

field17 = analysis.graphSTposmean ;
field17a = analysis_off.graphSTposmean;
combine.graphSTposmean = [field17, field17a];

field19 = analysis.graphNegexclude_excludeOver5ST ;
field19a = analysis_off.graphNegexclude_excludeOver5ST;
combine.graphNegexclude_excludeOver5ST = [field19, field19a];

%pos
field9 = analysis.graphPosKeepALL;
field9a = analysis_off.graphPosKeepALL;
combine.graphPosKeepALL = [field9, field9a];

field10 =analysis.graphPosexcludeALL;
field10a =analysis_off.graphPosexcludeALL;
combine.graphPosexcludeALL = [field10, field10a];

field11 = analysis.graphPosKeep_excludeOver5ST ;
field11a = analysis_off.graphPosKeep_excludeOver5ST;
combine.graphPosKeep_excludeOver5ST = [field11, field11a];



field13 = analysis.graphPosKeepFollow_excludeOver5ST;
field13a =analysis_off.graphPosKeepFollow_excludeOver5ST;
combine.graphPosKeepFollow_excludeOver5ST = [field13, field13a];

field14 = analysis.graphPosexcludeFollow_excludeOver5ST;
field14a =  analysis_off.graphPosexcludeFollow_excludeOver5ST;
combine.graphPosexcludeFollow_excludeOver5ST = [field14, field14a];

field15 = analysis.graphPosKeepOpp_excludeOver5ST;
field15a = analysis_off.graphPosKeepOpp_excludeOver5ST;
combine.graphPosKeepOpp_excludeOver5ST = [field15, field15a];

field16  = analysis.graphPosexcludeOpp_excludeOver5ST; 
field16a = analysis_off.graphPosexcludeOpp_excludeOver5ST ;
combine.graphPosexcludeOpp_excludeOver5ST = [field16, field16a];

field18 = analysis.graphSTnegmean ;
field18a = analysis_off.graphSTnegmean;
combine.graphSTnegmean = [field18, field18a];

field20 = analysis.graphPosexclude_excludeOver5ST ;
field20a = analysis_off.graphPosexclude_excludeOver5ST;
combine.graphPosexclude_excludeOver5ST = [field20, field20a];




%calculate relevant information

%make graphs to save
% % title ('Average all responses')
% % xlabel ('time (s)')
% % ylabel ('ST: compared to trial baseline')
% % legend('Shift up', 'Shift down')
% % line([0 0], get(gca, 'ylim'), 'color', 'black')
% % hold on
% % plot (t, combine.graphSTposmean, 'r')
% % plot(t, combine.graphSTposmean, 'b')
% % %savefig(gcf, 'averages')
% % hold off
% % close all

t = [-.2: 1/44101: 1]; %for graphing looking at -200 to 1000, with 0 being onset of perturbation

title ('Average all responses Excluded Over 5 ST')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

combine.baseline.Mean.neg = mean(combine.negAvg_excludeOver5ST(1:(length(combine.negAvg_excludeOver5ST)*.2)));
combine.baseline.SD.neg = std(combine.negAvg_excludeOver5ST(1:(length(combine.negAvg_excludeOver5ST)*.2)));

%start looking for peak 60 ms after onset of perturbation until 600 ms
%after pertubation onset
startindex = 0.06;
endindex = 0.6;

hold on
plot (t, combine.posAvg_excludeOver5ST, 'r','LineWidth', 2)
plot(t, combine.negAvg_excludeOver5ST, 'b','LineWidth', 2)
legend('Shift up', 'Shift down')
line([0 0], get(gca, 'ylim'), 'color', 'black')
line([0 0], get(gca, 'ylim'), 'color', 'black') %dumb way of making the line span the entire y axis limit

%mark start and end of analysis period
line([startindex startindex], get(gca, 'ylim'), 'Color', 'black', 'LineStyle', '--')
line([endindex endindex], get(gca, 'ylim'), 'Color', 'black', 'LineStyle', '--')


thresholdneg = 2*combine.baseline.SD.neg;
%thresholdpos = 2*combine.baseline.SD.pos;

analysisSection = combine.negAvg_excludeOver5ST((length(combine.negAvg_excludeOver5ST)*.26):(length(combine.negAvg_excludeOver5ST)*.8), :);

save('combine.mat', 'combine')







savefig(gcf, 'averages_ExcludeOver5ST')
hold off
close all


title ('shift down (greater then abs 5 ST removed')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

hold on
if ~isempty(combine.graphNegexclude_excludeOver5ST)
plot(t,combine.graphNegexclude_excludeOver5ST,'k-.')
end
if ~isempty(combine.graphNegKeep_excludeOver5ST)
plot(t,combine.graphNegKeep_excludeOver5ST, '-','LineWidth', 2)
end
line([0 0], get(gca, 'ylim'), 'color', 'black')
if ~isempty(combine.graphNegexclude_excludeOver5ST) || ~isempty(combine.graphNegKeep_excludeOver5ST)
savefig(gcf, 'ShiftDown_over5ST_excluded.fig')
else 
end
hold off
close all

title ('shift down (follow): over 5 ST excluded')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

hold on

if ~isempty(combine.graphNegexcludeFollow_excludeOver5ST)
plot(t,combine.graphNegexcludeFollow_excludeOver5ST,'k-.')
end

if ~isempty(combine.graphNegKeepFollow_excludeOver5ST)
plot(t,combine.graphNegKeepFollow_excludeOver5ST, '-','LineWidth', 2)
end
line([0 0], get(gca, 'ylim'), 'color', 'black')
if ~isempty(combine.graphNegexcludeFollow_excludeOver5ST) || ~isempty(combine.graphNegKeepFollow_excludeOver5ST)
savefig(gcf, 'ShiftDown_follow.fig')
else
end
hold off
close all

title ('shift down (Oppose) over 5 ST excluded')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

hold on

if ~isempty(combine.graphNegexcludeOpp_excludeOver5ST)
plot(t,combine.graphNegexcludeOpp_excludeOver5ST,'k-.')
end

if ~isempty(combine.graphNegKeepOpp_excludeOver5ST)
plot(t,combine.graphNegKeepOpp_excludeOver5ST, '-','LineWidth', 2)
end
line([0 0], get(gca, 'ylim'), 'color', 'black')
if ~isempty(combine.graphNegKeepOpp_excludeOver5ST) || ~isempty(combine.graphNegexcludeOpp_excludeOver5ST)
savefig(gcf, 'ShiftDown_oppose.fig')
else
end
hold off
close all



title ('shift up (greater then abs 5 ST removed')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

hold on

if ~isempty(combine.graphPosexclude_excludeOver5ST)
plot(t,combine.graphPosexclude_excludeOver5ST ,'k-.')
end

if ~isempty(combine.graphPosKeep_excludeOver5ST)
plot(t,combine.graphPosKeep_excludeOver5ST, '-','LineWidth', 2)
end
line([0 0], get(gca, 'ylim'), 'color', 'black')
if ~isempty(combine.graphPosexclude_excludeOver5ST) || ~isempty(combine.graphPosKeep_excludeOver5ST)
savefig(gcf, 'ShiftUp_over5ST_excluded.fig')
else
end
hold off
close all


title ('shift up (follow) over 5 ST excluded')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

hold on

if ~isempty(combine.graphPosexcludeFollow_excludeOver5ST)
plot(t,combine.graphPosexcludeFollow_excludeOver5ST,'k-.')
end

if ~isempty(combine.graphPosKeepFollow_excludeOver5ST)
plot(t,combine.graphPosKeepFollow_excludeOver5ST, '-','LineWidth', 2)
end
line([0 0], get(gca, 'ylim'), 'color', 'black')
if ~isempty(combine.graphPosexcludeFollow_excludeOver5ST) || ~isempty(combine.graphPosKeepFollow_excludeOver5ST)
savefig(gcf, 'ShiftUp_follow.fig')
else
end
hold off
close all

title ('shift up (oppose)over 5 ST excluded')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

hold on

if ~isempty(combine.graphPosexcludeOpp_excludeOver5ST)
plot(t,combine.graphPosexcludeOpp_excludeOver5ST,'k-.')
end

if ~isempty(combine.graphPosKeepOpp_excludeOver5ST)
plot(t,combine.graphPosKeepOpp_excludeOver5ST, '-','LineWidth', 2)
end

line([0 0], get(gca, 'ylim'), 'color', 'black')
if ~isempty(combine.graphPosexcludeOpp_excludeOver5ST) || ~isempty(combine.graphPosKeepOpp_excludeOver5ST)
savefig(gcf, 'ShiftUp_oppose.fig')
else
end

hold off
close all
