clear all
close all

files = dir('*.mat') ;
load(files.name);

f = fieldnames(keepdata.ch1data);

trialinvalid = find(cell2mat(keepdata.average(2:end,3))); %find nonzero elements that should be removed

%%initialize matrices
concatmatrix_baseline_neg = [];
concatmatrix_reflex_full_neg = [];
concatmatrix_reflex_neg = [];

concatmatrix_baseline_pos = [];
concatmatrix_reflex_full_pos = [];
concatmatrix_reflex_pos = [];

avgneg = [];
avgpos = [];

avgnegbase = [];
avgposbase = [];

STneg = [];
STpos = [];

STavgnegUSE = [];
STavgposUSE = [];

STnegbase = [];
STposbase = [];

%STnegfull = [];
%STposfull = [];

STnegUSE = [];
STposUSE = [];

oppNegValue = [];
followNegValue = [];

oppPosValue = [];
followPosValue = [];


oppNegUSE = [];
followNegUSE = [];

oppPosUSE = [];
followPosUSE = [];

STnegfull = [];
STposfull= [];



start = 44100*.15;
end2 = 44100*.3;

%% split trials into shift up and shift down
fulllist = [2:1:length(keepdata.average(:,1))]'; %number of trials based on their index - index 1 = trial header
negative = find(contains(keepdata.average(:,1), 'NEG')); %negative trials, anything that has "NEG" in the name
positive = setdiff(fulllist, negative); %all the other trials in the list that do not have NEG



%remove trials that are not valid
negative = setdiff(negative, trialinvalid);
positive = setdiff(positive, trialinvalid);

%change the indexing back to match  trial numbers
negative = negative-1;
positive = positive -1;
%trialinvalid = trialinvalid - 1;

%start with the assumption there are no usable pos or neg trials
negTrue = 0; 
posTrue = 0;
for i = 1:length([negative;positive;trialinvalid])
   
    if any(trialinvalid==i)
        %skip
    else
    extranan = [];
    fullneg = [];
    fullpos = [];
    if any(negative==i) %for all the negative trials
        negTrue = 1; %if entered this loop, there is a valid shiftdown (neg trial)
        concatmatrix_baseline_neg = [concatmatrix_baseline_neg keepdata.baselinech1data.(f{i})(:,3)]; %concatenate the F0 from -200 to 0 (start pert)
      
       %in case trial ended early
        extranan = [zeros(1, 44101 - length(keepdata.ch1data.(f{i})(:,3)))'];
        extranan(extranan==0) = NaN;
        fullneg = [keepdata.ch1data.(f{i})(:,3); extranan];
        
        concatmatrix_reflex_full_neg = [concatmatrix_reflex_full_neg fullneg]; %concat. the F0 form 0 to 1000 ms
        concatmatrix_reflex_neg = [concatmatrix_reflex_full_neg(start:end2, :)]; %concat the F0 from 150 to 300 ms
        
       
    elseif any(positive==i)
        posTrue = 1; %if entered this loop, there is a valid shiftup (pos trial)
        concatmatrix_baseline_pos = [concatmatrix_baseline_pos keepdata.baselinech1data.(f{i})(:,3)];
              
              
              
          %in case trial ended early
        extranan = [zeros(1, 44101 - length(keepdata.ch1data.(f{i})(:,3)))'];
        extranan(extranan==0) = NaN;
        fullpos = [keepdata.ch1data.(f{i})(:,3); extranan];
        
        
        concatmatrix_reflex_full_pos = [concatmatrix_reflex_full_pos fullpos];
        concatmatrix_reflex_pos = [concatmatrix_reflex_full_pos(start:end2, :)];
    end
    end
end

%% See "negTrue" loop for all notes: Same in pos in negative
if negTrue ==1 %if there is a negative trial, enter here
    
    %change zeros to nan's so don't get averaged in
    concatmatrix_reflex_full_neg(concatmatrix_reflex_full_neg==0)=nan;
    concatmatrix_baseline_neg(concatmatrix_baseline_neg==0)=nan;
    concatmatrix_reflex_pos(concatmatrix_reflex_pos==0)=nan;
    
    %calculate the average of the baseline
  avgnegbase = nanmean(concatmatrix_baseline_neg);
  
   %concat the baseline and the full reflex (-200 to 1000 ms)
    fullnegative = [concatmatrix_baseline_neg; concatmatrix_reflex_full_neg];

 for aa = 1:length(avgnegbase)
       STneg = [STneg 12*log2(concatmatrix_reflex_neg(:,aa)/avgnegbase(:,aa))]; %ST of 150 to 300 ms to determine direction
       STnegbase = [STnegbase 12*log2(concatmatrix_baseline_neg(:,aa)/avgnegbase(:,aa))]; %ST of baseline to determine variability
       STnegfull = [STnegfull 12*log2(fullnegative(:,aa)/avgnegbase(:,aa))]; 
 
 end
 
 STnegmean = nanmean(STnegfull,2); %calculate the mean across all trials, ignoring "nans"
 
 % determine direction of change for each trial

%directionneg = [mean(STneg) - mean(STnegbase)]';
idx_oppNeg = find(mean(STneg)>0);
idx_followNeg = find(mean(STneg)<0);

% find indicies where the baseline ST (in reference to the mean baseline) is
%greater than an 1 ST (*******currently a placeholder, not elimianting
%anything, but code is here so can add in a variable when want. Larson
%need to do sensitivity analysis for .1, .15. and .2 (maybe more?)

%calculate absolute value of average ST
maxabsSTnegbase = max(abs(STnegbase))';

%find indicies where max value is greater than 1 st (need to decide what this actual number is). If greater than need to
%eliminate because too variable
idx_maxabsSTnegbase = find(maxabsSTnegbase > 1); %%%placeholder number, change this when decide on criteria

oppNeg = setdiff(idx_oppNeg,idx_maxabsSTnegbase);
followNeg = setdiff(idx_followNeg,idx_maxabsSTnegbase);

% calculate ST in relation to average baseline
 for ee  = 1: length(avgnegbase)
     STnegUSE = [STnegUSE 12*log2(fullnegative(:,ee)/avgnegbase(:,ee))];
     
 end
 

 %calculate average ST at each data point across all trials within a
 %subject that are included 
  for hh = 1: length(STnegUSE)
     
     if any(idx_oppNeg==hh) %for all the negative trials

             oppNegValue = [oppNegValue STnegUSE(:, hh)];
         
     elseif any (idx_followNeg == hh)
         followNegValue = [followNegValue STnegUSE(:, hh)];
     else
     end
  end
 
   for cd = 1: size(oppNegValue,1)
     if ~isempty(oppNegValue)
    oppNegUSE = [oppNegUSE; mean(oppNegValue(cd,:))];
     else
         oppNegUSE = [];
     end
 end

 for cf = 1:size(followNegValue,1)
       if ~isempty(followNegValue)
    followNegUSE = [followNegUSE; mean(followNegValue(cf,:))];
       else
           followNegUSE = [];
       end
       
 end
 
oppNegUSE_downsamp = downsample(oppNegUSE, 10);
followNegUSE_downsamp= downsample(followNegUSE, 10);

% average regardless of direction
STnegmean = nanmean(STnegfull,2);


%% saving into structure

%shift down (negative trials)
analysis.neg.baseline.full = concatmatrix_baseline_neg;
analysis.neg.baseline.avg = avgnegbase;

analysis.neg.ST.baseline = STnegbase;
analysis.neg.ST.reflexwindow  = STneg;
analysis.neg.ST.average = STnegmean;
analysis.neg.ST.Full = STnegfull;

analysis.neg.reflex.full = concatmatrix_reflex_full_neg;
analysis.neg.reflex.window = concatmatrix_reflex_neg;

analysis.neg.baselinereflex = fullnegative;
analysis.neg.full.opp.avg =  oppNegUSE; 
analysis.neg.full.follow.avg =followNegUSE;
analysis.neg.full.opp.values =  oppNegValue; 
analysis.neg.full.follow.values =followNegValue;
analysis.neg.downsamp.opp = oppNegUSE_downsamp;
analysis.neg.downsamp.follow = followNegUSE_downsamp;

analysis.neg.idx.removebase = idx_maxabsSTnegbase;
%analysis.neg.idx.direction = directionneg;
analysis.neg.idx.opp = oppNeg;
analysis.neg.idx.follow = followNeg;
else
  analysis.neg.baseline.full = [];
analysis.neg.baseline.avg = [];

analysis.neg.ST.baseline = [];
analysis.neg.ST.reflexwindow  = [];
analysis.neg.ST.average = [];
analysis.neg.ST.Full = [];

analysis.neg.reflex.full = [];
analysis.neg.reflex.window = [];

analysis.neg.baselinereflex = [];
analysis.neg.full.opp.avg =  []; 
analysis.neg.full.follow.avg =[];
analysis.neg.full.opp.values =  []; 
analysis.neg.full.follow.values =[];
analysis.neg.downsamp.opp = [];
analysis.neg.downsamp.follow = [];

analysis.neg.idx.removebase = [];
%analysis.neg.idx.direction = [];
analysis.neg.idx.opp = [];
analysis.neg.idx.follow = [];

end

%%
if posTrue ==1 
    concatmatrix_reflex_full_pos(concatmatrix_reflex_full_pos==0)=nan;
    concatmatrix_baseline_pos(concatmatrix_baseline_pos==0)=nan;
    concatmatrix_reflex_neg(concatmatrix_reflex_neg==0)=nan;
    
    avgposbase = nanmean(concatmatrix_baseline_pos);
    
    fullpositive = [concatmatrix_baseline_pos; concatmatrix_reflex_full_pos];

    
    for   bb = 1:length(avgposbase)
       STpos = [STpos 12*log2(concatmatrix_reflex_pos(:,bb)/avgposbase(:,bb))]; %ST of 150 to 300 ms to determine direction
       STposbase = [STposbase 12*log2(concatmatrix_baseline_pos(:,bb)/avgposbase(:,bb))]; %ST of baseline to determine variability
       STposfull = [STposfull 12*log2(fullpositive(:,bb)/avgposbase(:,bb))];
 
    end
    
    STposmean = nanmean(STposfull,2);
    %%
    %directionpos = [mean(STpos) - mean(STposbase)]';
 idx_oppPos = find(mean(STpos)<0);
idx_followPos  = find(mean(STpos)>0);

maxabsSTposbase = max(abs(STposbase))';
idx_maxabsSTposbase = find(maxabsSTposbase > 1); %%%placeholder number, change this when decide on criteria

oppPos = setdiff(idx_oppPos,idx_maxabsSTposbase);
followPos = setdiff(idx_followPos,idx_maxabsSTposbase);

 for ff = 1: length(avgposbase)
        STposUSE = [STposUSE 12*log2(fullpositive(:,ff)/avgposbase(:,ff))]; 
 end
 
 
   for jj = 1: length(STposUSE)
     
     if any(idx_oppPos==jj) %for all the negative trials

             oppPosValue = [oppPosValue STposUSE(:, jj)];
         
     elseif any (idx_followPos == jj)
         followPosValue = [followPosValue STposUSE(:, jj)];
     else
     end
 end
 
 for cg = 1: size(oppPosValue,1)
    if ~isempty(oppPosValue)
     oppPosUSE = [oppPosUSE; mean(oppPosValue(cg,:))];
    else
        oppPosUSE = [];
    
    end
 end
 
  for ce = 1:size(followPosValue,1) 
      if ~isempty(followPosValue)
    followPosUSE = [followPosUSE; mean(followPosValue(ce,:))];
      else
          followPosUSE = [];
      end
  end
 
  oppPosUSE_downsamp = downsample(oppPosUSE, 10);
followPosUSE_downsamp = downsample(followPosUSE, 10);

STposmean =  nanmean(STposfull,2);




%shift up (positive trials)
analysis.pos.baseline.full = concatmatrix_baseline_pos;
analysis.pos.baseline.avg = avgposbase;

analysis.pos.ST.baseline = STposbase;
analysis.pos.ST.reflexwindow  = STpos;
analysis.pos.ST.average = STposmean;
analysis.pos.ST.Full = STposfull;

analysis.pos.reflex.full = concatmatrix_reflex_full_pos;
analysis.pos.reflex.window = concatmatrix_reflex_pos;

analysis.pos.baselinereflex = fullpositive; %%use this for looking at HZ over the -200 to 1000, 0 is onset of perturbation
analysis.pos.full.opp.avg =  oppPosUSE; 
analysis.pos.full.follow.avg =followPosUSE;
analysis.pos.full.opp.values =  oppPosValue; 
analysis.pos.full.follow.values =followPosValue;
analysis.pos.downsamp.opp.avg = oppPosUSE_downsamp;
analysis.pos.downsamp.follow.avg = followPosUSE_downsamp;

analysis.pos.idx.removebase = idx_maxabsSTposbase;
%analysis.pos.idx.direction = directionpos;
analysis.pos.idx.opp = oppPos;
analysis.pos.idx.follow = followPos;


else
    %do nothing
    %shift up (positive trials)
analysis.pos.baseline.full = [];
analysis.pos.baseline.avg = [];

analysis.pos.ST.baseline = [];
analysis.pos.ST.reflexwindow  = [];
analysis.pos.ST.average = [];
analysis.pos.ST.Full = [];

analysis.pos.reflex.full = [];
analysis.pos.reflex.window = [];

analysis.pos.baselinereflex = []; %%use this for looking at HZ over the -200 to 1000, 0 is onset of perturbation
analysis.pos.full.opp.avg =  []; 
analysis.pos.full.follow.avg =[];
analysis.pos.full.opp.values =  []; 
analysis.pos.full.follow.values =[];
analysis.pos.downsamp.opp.avg = [];
analysis.pos.downsamp.follow.avg = [];

analysis.pos.idx.removebase = [];
%analysis.pos.idx.direction = [];
analysis.pos.idx.opp = [];
analysis.pos.idx.follow = [];

end


%% graphing data



t = [-.2: 1/44101: 1]; %for graphing

if negTrue == 1 %shiftdown
    
%% all shift down
title ('shift down')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on
excludedgroup = [];
keepgroup = [];

for ga = 1:size(followNegValue,2)
 if any(idx_maxabsSTnegbase==ga)
     excludedgroup = [excludedgroup followNegValue(:,ga)];
     
 else
     keepgroup = [keepgroup followNegValue(:,ga)];
 end
end

for gb = 1:size(oppNegValue,2)
 if any(idx_maxabsSTnegbase==gb) 
     excludedgroup = [excludedgroup oppNegValue(:,gb)];
 else
     keepgroup = [keepgroup oppNegValue(:,gb)];
 end
end

if ~isempty(excludedgroup)
plot(t,excludedgroup,'k-.')
analysis.graphNegexcludeALL = excludedgroup;
else 
    analysis.graphNegexcludeALL = [];
end
hold on

if ~isempty(keepgroup)
plot(t,keepgroup, '-','LineWidth', 2)
analysis.graphNegKeepALL = keepgroup;

else
    analysis.graphNegKeepALL = [];
end

line([0 0], get(gca, 'ylim'), 'color', 'black')
if ~isempty(excludedgroup) || ~isempty(keepgroup)
savefig(gcf, 'ShiftDown.fig')
else
end

hold off
close all

%% shift down with trials removed greater than 5 ST, assuming pitch tracking error
figure
hold on
title ('shift down (greater then abs 5 ST removed')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

if ~isempty(excludedgroup)
plot(t,excludedgroup,'k-.')
analysis.graphNegexclude_excludeOver5ST = excludedgroup;
else
  analysis.graphNegexclude_excludeOver5ST = [];  
end
hold on

if ~isempty(keepgroup)
keepgroup(keepgroup<-5)=nan;
keepgroup(keepgroup>5)=nan;
analysis.graphNegKeep_excludeOver5ST = keepgroup;
plot(t,keepgroup, '-','LineWidth', 2)
else
    analysis.graphNegKeep_excludeOver5ST = [];
end
line([0 0], get(gca, 'ylim'), 'color', 'black')

if ~isempty(excludedgroup) || ~isempty(keepgroup)
savefig(gcf, 'ShiftDown_over5ST_excluded.fig')
else 
end


hold off
close all
else
     analysis.graphNegKeepALL = [];
    analysis.graphNegexcludeALL = [];
    analysis.graphNegKeep_excludeOver5ST = [];

    analysis.graphNegKeepFollow_excludeOver5ST  = [];
    analysis.graphNegexcludeFollow_excludeOver5ST = [];
    analysis.graphNegKeepOpp_excludeOver5ST  = [];
    analysis.graphNegexcludeOpp_excludeOver5ST = [];  
      analysis.graphNegexclude_excludeOver5ST = [];
end

if posTrue ==1 %shift up

%shift up, all trials    
title ('shift up')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

hold on
excludedgroup = [];
keepgroup = [];
for gc = 1:size(followPosValue,2)
 if any(idx_maxabsSTposbase==gc) %if there are any trials with variable baselines > .15 semitone
     excludedgroup = [excludedgroup followPosValue(:,gc)];
 else
     keepgroup = [keepgroup followPosValue(:,gc)];
 end
end
   

for gd = 1:size(oppPosValue,2)
 if any(idx_maxabsSTposbase==gd) %if there are any trials with variable baselines > .15 semitone
     excludedgroup = [excludedgroup oppPosValue(:,gd)];
 else
     keepgroup = [keepgroup oppPosValue(:,gd)];
 end
end
   
if ~isempty(excludedgroup)
plot(t,excludedgroup,'k-.')
analysis.graphPosexcludeALL = excludedgroup;
else
  analysis.graphPosexcludeALL = [];  
end
hold on

if ~isempty(keepgroup)
plot(t,keepgroup, '-','LineWidth', 2)
analysis.graphPosKeepALL = keepgroup;
else
    analysis.graphPosKeepALL = [];
end

line([0 0], get(gca, 'ylim'), 'color', 'black')

if ~isempty(excludedgroup) || ~isempty(keepgroup)
savefig(gcf, 'ShiftUp.fig')
else
end

hold off
close all

%% excluding values over 5 ST, assuming pitch tracking error
figure
hold on
title ('shift up (greater then abs 5 ST removed')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')

if ~isempty(excludedgroup)
plot(t,excludedgroup,'k-.')
analysis.graphPosexclude_excludeOver5ST = excludedgroup;
else
  analysis.graphPosexclude_excludeOver5ST = [];  
end
hold on


if ~isempty(keepgroup)
keepgroup(keepgroup<-5)=nan;
keepgroup(keepgroup>5)=nan;
plot(t,keepgroup, '-','LineWidth', 2)
analysis.graphPosKeep_excludeOver5ST = keepgroup;
else
    analysis.graphPosKeep_excludeOver5ST = [];
end
line([0 0], get(gca, 'ylim'), 'color', 'black')

if ~isempty(excludedgroup) || ~isempty(keepgroup)
savefig(gcf, 'ShiftUp_over5ST_excluded.fig')
else
end

hold off
close all
else
     analysis.graphPosKeepALL = [];
    analysis.graphPosexcludeALL = [];
    analysis.graphPosKeep_excludeOver5ST = [];
    analysis.graphPosKeepFollow_excludeOver5ST  = [];
    analysis.graphPosexcludeFollow_excludeOver5ST = [];
    analysis.graphPosKeepOpp_excludeOver5ST  = [];
    analysis.graphPosexcludeOpp_excludeOver5ST = [];
    analysis.graphPosexclude_excludeOver5ST = []; 
end

%save averages of both if available
%save ST means if available
figure
if exist('STposmean') ~=0
plot (t, STposmean, 'r')
analysis.graphSTposmean = STposmean;
else
    analysis.graphSTposmean = [];
    STposmean = [];
end
hold on
if exist('STnegmean') ~= 0
plot(t, STnegmean, 'b')
analysis.graphSTnegmean = STnegmean;
else
    analysis.graphSTnegmean = [];
    STnegmean = [];
end
line([0 0], get(gca, 'ylim'), 'color', 'black')
title ('Average all responses')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
legend('Shift up', 'Shift down')

if ~isempty(STnegmean) || ~isempty(STposmean)
savefig(gcf, 'averages')
else
end

hold off
close all

if posTrue == 1
if ~isempty(analysis.graphPosKeep_excludeOver5ST) 
analysis.posAvg_excludeOver5ST = nanmean(analysis.graphPosKeep_excludeOver5ST,2);
hold on
plot (t, analysis.posAvg_excludeOver5ST, 'r')
else
    analysis.posAvg_excludeOver5ST = [];
end
else
    analysis.posAvg_excludeOver5ST = [];
end

if negTrue ==1
if ~isempty(analysis.graphNegKeep_excludeOver5ST) 
analysis.negAvg_excludeOver5ST = nanmean(analysis.graphNegKeep_excludeOver5ST,2);
plot(t, analysis.negAvg_excludeOver5ST, 'b')
else
    analysis.negAvg_excludeOver5ST = [];
end
else
     analysis.negAvg_excludeOver5ST = [];
end


line([0 0], get(gca, 'ylim'), 'color', 'black')
title ('Average all responses Excluded Over 5 ST')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
legend('Shift up', 'Shift down')

if posTrue ==1 || negTrue ==1
savefig(gcf, 'averages_ExcludeOver5ST')
else 
end

hold off 
close all


%%%% seperate follow and opp graphs
if negTrue == 1
title ('shift down (follow): over 5 ST excluded')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on

%all following shift down (split up here in case want to differentiate them
%later
excludedgroup = [];
keepgroup = [];

for ga = 1:size(followNegValue,2)
 if any(idx_maxabsSTnegbase==ga)
     excludedgroup = [excludedgroup followNegValue(:,ga)];
     
 else
     keepgroup = [keepgroup followNegValue(:,ga)];
 end
end


if ~isempty(excludedgroup)
plot(t,excludedgroup,'k-.')
analysis.graphNegexcludeFollow_excludeOver5ST = excludedgroup;
else
   analysis.graphNegexcludeFollow_excludeOver5ST = []; 
end

hold on
if ~isempty(keepgroup)
keepgroup(keepgroup<-5)=nan;
keepgroup(keepgroup>5)=nan;
plot(t,keepgroup, '-','LineWidth', 2)
analysis.graphNegKeepFollow_excludeOver5ST = keepgroup;
else
analysis.graphNegKeepFollow_excludeOver5ST = [];
end


line([0 0], get(gca, 'ylim'), 'color', 'black')

if ~isempty(keepgroup) || ~isempty(excludedgroup)
savefig(gcf, 'ShiftDown_follow.fig')
else
end

hold off
close all

title ('shift down (Oppose) over 5 ST excluded')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on

excludedgroup = [];
keepgroup = [];

for gb = 1:size(oppNegValue,2)
 if any(idx_maxabsSTnegbase==gb) %if there are any trials with variable baselines > .15 semitone
     excludedgroup = [excludedgroup oppNegValue(:,gb)];
 else
     keepgroup = [keepgroup oppNegValue(:,gb)];
 end
end


if ~isempty(excludedgroup)
plot(t,excludedgroup,'k-.')
analysis.graphNegexcludeOpp_excludeOver5ST = excludedgroup;
else
    analysis.graphNegexcludeOpp_excludeOver5ST = [];
end
hold on

if ~isempty(keepgroup)
keepgroup(keepgroup<-5)=nan;
keepgroup(keepgroup>5)=nan;
plot(t,keepgroup, '-','LineWidth', 2)
analysis.graphNegKeepOpp_excludeOver5ST = keepgroup;
else
    analysis.graphNegKeepOpp_excludeOver5ST = [];

end
line([0 0], get(gca, 'ylim'), 'color', 'black')

if ~isempty(keepgroup) || ~isempty(excludedgroup)
savefig(gcf, 'ShiftDown_oppose.fig')
else
end

hold off
close all

else
end

if posTrue ==1
title ('shift up (follow) over 5 ST excluded')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on
excludedgroup = [];
keepgroup = [];

for gc = 1:size(followPosValue,2)
 if any(idx_maxabsSTposbase==gc) %if there are any trials with variable baselines > .15 semitone
     excludedgroup = [excludedgroup followPosValue(:,gc)];
 else
     keepgroup = [keepgroup followPosValue(:,gc)];
 end
end

if ~isempty(excludedgroup)
plot(t,excludedgroup,'k-.')
analysis.graphPosexcludeFollow_excludeOver5ST = excludedgroup;
else
    analysis.graphPosexcludeFollow_excludeOver5ST = [];
end
hold on

if ~isempty(keepgroup)
keepgroup(keepgroup<-5)=nan;
keepgroup(keepgroup>5)=nan;
plot(t,keepgroup, '-','LineWidth', 2)
analysis.graphPosKeepFollow_excludeOver5ST = keepgroup;
else
    analysis.graphPosKeepFollow_excludeOver5ST = [];

end
line([0 0], get(gca, 'ylim'), 'color', 'black')

if ~isempty(excludedgroup) || ~isempty(keepgroup)
savefig(gcf, 'ShiftUp_follow.fig')
else
end



hold off
close all


title ('shift up (oppose)over 5 ST excluded')
xlabel ('time (s)')
ylabel ('ST: compared to trial baseline')
hold on
excludedgroup = [];
keepgroup = [];
for gd = 1:size(oppPosValue,2)
 if any(idx_maxabsSTposbase==gd) %if there are any trials with variable baselines > .15 semitone
     excludedgroup = [excludedgroup oppPosValue(:,gd)];
 else
     keepgroup = [keepgroup oppPosValue(:,gd)];
 end
end



if ~isempty(excludedgroup)
plot(t,excludedgroup,'k-.')
analysis.graphPosexcludeOpp_excludeOver5ST = excludedgroup;
else
 analysis.graphPosexcludeOpp_excludeOver5ST = [];
end
hold on

if ~isempty(keepgroup)
keepgroup(keepgroup<-5)=nan;
keepgroup(keepgroup>5)=nan;
plot(t,keepgroup, '-','LineWidth', 2)
analysis.graphPosKeepOpp_excludeOver5ST = keepgroup;
else
analysis.graphPosKeepOpp_excludeOver5ST = [];
end

line([0 0], get(gca, 'ylim'), 'color', 'black')

if ~isempty(keepgroup) || ~isempty(excludedgroup)
savefig(gcf, 'ShiftUp_oppose.fig')
else
end

hold off
close all

else
        
end

save('analysis.mat', 'analysis')
hold off
close all



