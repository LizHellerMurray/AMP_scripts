%checking follow/oppose by subtracting 200 ms - 0 ms value and seeing
%direction
%2 analysis types
%#1 --> looking for peak between 60 ms to 600 ms, latency is time of that peak
%#2 --> after 60 ms, find a part that is 2 SD away from mean, that is the
%latency. Then from that point to 200 ms after take average

clear all
close all


combinedalready =  input('is this combined(1) or analysis(2) mat file: ');
if combinedalready == 1
    load combine.mat
    
    
    combinepos = nanmean(combine.graphPosKeep_excludeOver5ST,2);
 combineneg = nanmean(combine.graphNegKeep_excludeOver5ST,2);
  combinestdpos = nanstd(combine.graphPosKeep_excludeOver5ST,0,2);
 combeinestdneg = nanstd(combine.graphNegKeep_excludeOver5ST,0,2);
 allNeg = combine.graphNegKeep_excludeOver5ST;
allPos = combine.graphPosKeep_excludeOver5ST;
 
 save negAvg_excludeOver5ST combineneg
 save posAvg_excludeOver5ST combinepos
  save posSTD combinestdpos
 save negSTD combeinestdneg
else
load analysis.mat
combinepos = nanmean(analysis.graphPosKeep_excludeOver5ST,2);
 combineneg = nanmean(analysis.graphNegKeep_excludeOver5ST,2);
 combinestdpos = nanstd(analysis.graphNegKeep_excludeOver5ST,0,2);
 combeinestdneg = nanstd(analysis.graphNegKeep_excludeOver5ST,0,2);
 allNeg = analysis.graphNegKeep_excludeOver5ST;
allPos = analysis.graphPosKeep_excludeOver5ST;

 save negAvg_excludeOver5ST combineneg
 save posAvg_excludeOver5ST combinepos
 save posSTD combinestdpos
 save negSTD combeinestdneg
end
graphNeg = nanmean(allNeg);
graphPos = nanmean(allPos);
t = [-.2: 1/44101: 1]; %for graphing looking at -200 to 1000, with 0 being onset of perturbation

title ('USE: Average all responses')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

load('negAvg_excludeOver5ST.mat')
load('posAvg_excludeOver5ST.mat')

%baseline is first 200 ms (time -200 to 0)
baselinenegmean = mean(combineneg(1:(length(combineneg)*.2)));
baselinenegSD = mean(combeinestdneg(1:(length(combeinestdneg)*.2)));
baselineposmean = mean(combinepos(1:(length(combinepos)*.2)));
baselineposSD = mean(combinestdpos(1:(length(combinestdpos)*.2)));

threshnegSD = std(combineneg(1:(length(combineneg)*.2)));
threshposSD = std(combinepos(1:(length(combinepos)*.2)));

%looking for direction 
startindex = 0.0; %0 ms after start of perturbation 
endindex = 0.20; %check 200 ms period
hold on
plot (t, combinepos, 'r','LineWidth', 2)
plot(t,combineneg, 'b','LineWidth', 2)
legend('Shift up', 'Shift down')
line([0 0], get(gca, 'ylim'), 'color', 'black')
line([0 0], get(gca, 'ylim'), 'color', 'black') %dumb way of making the line span the entire y axis limit
%line([startindex startindex], get(gca, 'ylim'), 'Color', 'black', 'LineStyle', '--')
%line([endindex endindex], get(gca, 'ylim'), 'Color', 'black', 'LineStyle', '--')

%looking for time 0 to 200 ms (index start at time - 200 so have to move
%up)
DIRstart_index = (startindex+.2) * 44100;
DIRend_index = (endindex+.2) * 44100;


%find average of direction in 60-260 ms
%posDir = mean(combinepos(DIRstart_index:DIRend_index));
%negDir = mean(combineneg(DIRstart_index:DIRend_index));

posDir = combinepos(DIRend_index) - combinepos(DIRstart_index); %shift up (Expect be going neg)
negDir = combineneg(DIRend_index) - combineneg(DIRstart_index); %shift down (Expect going pos)



if posDir < 0
    SUdirection =1; %negative direction, oppose
else
    SUdirection =2; %positive direction, follow
end

if negDir > 0
    SDdirection =1; %positive direction, oppose
else
    SDdirection = 2; %negative direction, follow
end


%SUdirection = input('Does the shiftup 1st peak oppose/neg (1) or follow/pos (2): ');
%SDdirection = input('Does the shiftdown 1st peak oppose/pos (1) or follow/neg (2): ');

%find the first point that deviates 2 SD away from the baseline (based on
%direction in first 60-260ms)



%mark start and end of analysis period
%line([startindex startindex], get(gca, 'ylim'), 'Color', 'black', 'LineStyle', '--')
%line([endindex endindex], get(gca, 'ylim'), 'Color', 'black', 'LineStyle', '--')


%mark 2 SD for pos and neg
if SDdirection == 1 %opposing shift down response
    thresholdneg = baselinenegmean + (2*threshnegSD);

else %following shift down response
    thresholdneg = baselinenegmean - (2*threshnegSD);

end

if SUdirection ==1 %opposing shiftup response
thresholdpos = baselineposmean - (2*threshposSD);
else %followingshift up response
thresholdpos = baselineposmean + (2*threshposSD);
end




line( get(gca, 'xlim'),  [thresholdneg thresholdneg],'Color', 'blue', 'LineStyle', '--')
line( get(gca, 'xlim'),  [thresholdneg thresholdneg],'Color', 'blue', 'LineStyle', '--') %dumb way to make it span the whole thing
line( get(gca, 'xlim'),  [thresholdpos thresholdpos],'Color', 'red', 'LineStyle', '--')

%analysisstart = input('start of analysis section: ');
%analysisend = input('end of analysis section: ');

%% peak magnitude and latency at that time 
analysisstart = 0.06; %60 ms
analysisend = 0.6; %600 ms

analysisstart_index = (analysisstart+.2) * 44100; % make it start at 0, so 60ms after pert onset
analysisend_index = (analysisend+.2) * 44100; %600 ms after pert onset


posanalysis = [];    
neganalysis = [];
posanalysis = [combinepos(analysisstart_index:analysisend_index)];
neganalysis = [combineneg(analysisstart_index:analysisend_index)];

if SUdirection ==1
Magpos = min(posanalysis);
Magposindex = find(posanalysis == Magpos,1);
latencypos = (Magposindex/44100);
else
Magpos = max(posanalysis);
Magposindex = (find(posanalysis == Magpos,1));
latencypos = (Magposindex/44100);

end


if SDdirection ==1
Magneg = max(neganalysis);
Magnegindex = find(neganalysis == Magneg,1);
latencyneg = (Magnegindex/44100);
else
Magneg = min(neganalysis);
Magnegindex = find(neganalysis == Magneg,1);
latencyneg = (Magnegindex/44100);
end

%change index for plottting because analysis starts at 60 ms
plot((latencypos+.06), Magpos, '.', 'Color', 'black', 'MarkerSize', 25)
plot((latencyneg+.06), Magneg, '.', 'Color', 'black', 'MarkerSize', 25)



datafile.shiftupMagMax = Magpos;
datafile.shiftupLatencyMax = latencypos+0.06; %latency relative to perturbation start
datafile.shiftupVariability = baselineposSD;

datafile.shiftdownMagMax = Magneg;
datafile.shiftdownLatencyMax = latencyneg+0.06; %latency relative to perturbation start
datafile.shiftdownVariability = baselinenegSD;

%% find latency as 1st point over 2 SD, magnitude is 200 ms average after latency

%look for latency between 60 and 600 ms
fineLatstart = (0.06+.2) *44100; %start or pert period line
findLatend = (0.6+.2) *44100;

posLat = [];    
negLat = [];
posLat = [combinepos(fineLatstart:findLatend)];
negLat = [combineneg(fineLatstart:findLatend)];


if SUdirection ==1
FirstLatencyPos = (((find(posLat < thresholdpos,1))/44100)+0.06); %latency relative to perturbation onset
else
FirstLatencyPos = (((find(posLat > thresholdpos,1))/44100)+0.06);
end


if SDdirection ==1
FirstLatencyNeg = (((find(negLat > thresholdneg,1))/44100)+0.06);
else
FirstLatencyNeg = (((find(negLat < thresholdneg,1))/44100)+0.06);
end

%after find latency, fine 200-400 ms afterward to new analysis period
analysisLatStartPos = (FirstLatencyPos+.2)*44100;
analysisLatEndPos = (FirstLatencyPos + .4)*44100;
analysisLatStartNeg = (FirstLatencyNeg+.2)*44100;
analysisLatEndNeg = (FirstLatencyNeg + .4)*44100;

posLatMag = [];
negLatMag = [];
posLatMag = mean(combinepos(analysisLatStartPos:analysisLatEndPos));
negLatMag = mean(combineneg(analysisLatStartNeg:analysisLatEndNeg));



 plot(FirstLatencyPos, posLatMag, '*', 'Color', 'red', 'MarkerSize', 25)
 plot(FirstLatencyNeg, negLatMag, '*', 'Color', 'blue', 'MarkerSize', 25)


datafile.shiftupMagavg = posLatMag;
datafile.shiftupLatency = FirstLatencyPos;


datafile.shiftdownMagavg = negLatMag;
datafile.shiftdownLatency = FirstLatencyNeg;

savefig(gcf, 'Average_USE.fig')
hold off 
close all

figure
hold on
title ({'shift down : n=' num2str(size(allNeg,2))})
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on
plot(t,allNeg, '-','LineWidth', 2)
line([0 0], get(gca, 'ylim'), 'color', 'black')
line([0 0], get(gca, 'ylim'), 'color', 'black')
savefig(gcf, 'ShiftDown')
hold off
close all

figure
hold on
title ({'shift Up : n=' num2str(size(allPos,2))})
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on
plot(t,allPos, '-','LineWidth', 2)
line([0 0], get(gca, 'ylim'), 'color', 'black')
line([0 0], get(gca, 'ylim'), 'color', 'black')
savefig(gcf, 'ShiftUp')
hold off
close all

% figure
% hold on
% title ('Average all responses')
% xlabel ('time (s)')
% ylabel ('ST: compared to trial baseline')
% legend('Shift up', 'Shift down')
% plot(t,graphPos, 'r', '-','LineWidth', 2)
% plot(t,graphNeg, 'b', '-','LineWidth', 2)

%define start and end points looking for following/opposing
startDir = 0; %0 ms, start of perturbation
endDir = .2; %200 ms after start

startindex = (startDir +0.2)*44100;
endindex = (endDir + 0.2)*44100;

negOpp = [];
negFollow = [];

%run through all neg (shiftdown) trials
for i = 1:size(allNeg,2)
negDir = allNeg(endindex, i) - allNeg(startindex,i);

if negDir > 0 %oppose
    negOpp = [negOpp allNeg(:,i)];
else
    negFollow = [negFollow allNeg(:,i)];
end

end

figure
hold on
title ({'shift down follow: n=' num2str(size(negFollow,2))})
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on
plot(t,negFollow, '-','LineWidth', 2)
line([0 0], get(gca, 'ylim'), 'color', 'black')
line([0 0], get(gca, 'ylim'), 'color', 'black')
savefig(gcf, 'ShiftDown_follow')
hold off
close all

figure
hold on2) 
title ({'shift down oppose: n=' num2str(size(negOpp,2))})
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on
plot(t,negOpp, '-','LineWidth', 2)
line([0 0], get(gca, 'ylim'), 'color', 'black')
line([0 0], get(gca, 'ylim'), 'color', 'black')
savefig(gcf, 'ShiftDown_oppose')
hold off
close all


posOpp = [];
posFollow = [];
posNull = [];



for i = 1:size(allPos,2)
    posDir = allPos(endindex, i) - allPos(startindex,i);
    
    if posDir < 0 %oppose
        
        posOpp = [posOpp allPos(:,i)];
    else
        posFollow = [posFollow allPos(:,i)];
    end
    
end

figure
hold on
title ({'shift up follow: n=' num2str(size(posFollow,2))})
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
plot(t,posFollow, '-','LineWidth', 2)
line([0 0], get(gca, 'ylim'), 'color', 'black')
line([0 0], get(gca, 'ylim'), 'color', 'black')
savefig(gcf, 'ShiftUp_follow')
hold off
close all

figure
hold on
title ({'shift up oppose: n=' num2str(size(posOpp,2))})
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
plot(t,posOpp, '-','LineWidth', 2)
line([0 0], get(gca, 'ylim'), 'color', 'black')
line([0 0], get(gca, 'ylim'), 'color', 'black')
savefig(gcf, 'ShiftUp_oppose')
hold off
close all

save datafile datafile



%analysisSectionneg = combineneg((length(combineneg)*.26):(length(combineneg)*.8), :);
%analysisSectionpos = combinepos((length(combinepos)*.26):(length(combinepos)*.8), :);