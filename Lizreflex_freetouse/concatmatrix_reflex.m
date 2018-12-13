clear all
close all

files = dir('*.mat') ;
load(files.name);

f = fieldnames(keepdata.ch1data);

trialinvalid = find(cell2mat(keepdata.average(2:end,3))); %find nonzero elements that should be removed

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

start = 44100*.15;
end2 = 44100*.3;


fulllist = [2:1:length(keepdata.average(:,1))]'; %number of trials based on their index - index 1 = trial header
negative = find(contains(keepdata.average(:,1), 'NEG')); %negative trials, anything that has "NEG" in the name
positive = setdiff(fulllist, negative); %all the other trials in the list that do not have NEG



%remove trials that are not valid
negative = setdiff(negative, trialinvalid);
positive = setdiff(positive, trialinvalid);

for i = 1:length(f)
    if any(negative==i) %for all the negative trials
        concatmatrix_baseline_neg = [concatmatrix_baseline_neg keepdata.baselinech1data.(f{i})(:,3)]; %concatenate the F0 from -200 to 0 (start pert)
      
        concatmatrix_reflex_full_neg = [concatmatrix_reflex_full_neg keepdata.ch1data.(f{i})(:,3)]; %concat. the F0 form 0 to 1000 ms
        concatmatrix_reflex_neg = [concatmatrix_reflex_full_neg(start:end2, :)]; %concat the F0 from 150 to 300 ms
       
    elseif any(positive==i)
        concatmatrix_baseline_pos = [concatmatrix_baseline_pos keepdata.baselinech1data.(f{i})(:,3)];
        
        concatmatrix_reflex_full_pos = [concatmatrix_reflex_full_pos keepdata.ch1data.(f{i})(:,3)];
        concatmatrix_reflex_pos = [concatmatrix_reflex_full_pos(start:end2, :)];
    end
    
end

%calculate the average of the baseline
  avgnegbase = mean(concatmatrix_baseline_neg);
  avgposbase = mean(concatmatrix_baseline_pos);
  
 %concat the baseline and the full reflex (-200 to 1000 ms)
fullpositive = [concatmatrix_baseline_pos; concatmatrix_reflex_full_pos];
fullnegative = [concatmatrix_baseline_neg; concatmatrix_reflex_full_neg];
  
 for aa = 1:length(avgnegbase)
       STneg = [STneg 12*log2(concatmatrix_reflex_neg(:,aa)/avgnegbase(:,aa))]; %ST of 150 to 300 ms to determine direction
       STnegbase = [STnegbase 12*log2(concatmatrix_baseline_neg(:,aa)/avgnegbase(:,aa))]; %ST of baseline to determine variability
      % STnegfull = [STnegfull 12*log2(concatmatrix_reflex_full_neg(:,aa)/avgnegbase(:,aa))]; 
 
 end
       
for   bb = 1:length(avgposbase)
       STpos = [STpos 12*log2(concatmatrix_reflex_pos(:,bb)/avgposbase(:,bb))]; %ST of 150 to 300 ms to determine direction
       STposbase = [STposbase 12*log2(concatmatrix_baseline_pos(:,bb)/avgposbase(:,bb))]; %ST of baseline to determine variability
      % STposfull = [STposfull 12*log2(concatmatrix_reflex_full_pos(:,bb)/avgposbase(:,bb))];
 
end



%% determine direction of change for each trial

directionneg = [mean(STneg) - mean(STnegbase)]';
directionpos = [mean(STpos) - mean(STposbase)]';

idx_oppNeg = find(directionneg>0);
idx_followNeg = find(directionneg<0);

idx_oppPos = find(directionpos<0);
idx_followPos  = find(directionpos>0);

%% find indicies where the baseline ST (in reference to the mean baseline) is
%greater than an .15 ST

%calculate absolute value of average ST
maxabsSTnegbase = max(abs(STnegbase))';
maxabsSTposbase = max(abs(STposbase))';

%find indicies where max value is less than .15. If greater than need to
%eliminate because too variable
idx_maxabsSTnegbase = find(maxabsSTnegbase > 0.15);
idx_maxabsSTposbase = find(maxabsSTposbase > 0.15);


oppNeg = setdiff(idx_oppNeg,idx_maxabsSTnegbase);
followNeg = setdiff(idx_followNeg,idx_maxabsSTnegbase);

oppPos = setdiff(idx_oppPos,idx_maxabsSTposbase);
followPos = setdiff(idx_followPos,idx_maxabsSTposbase);



% for cc = 1: size(STpos,1)
%     STavgposUSE = [STavgposUSE; mean(STpos(cc,:))];
% end
% 
%  for dd = 1: size(STneg,1)
%     STavgnegUSE = [STavgnegUSE; mean(STneg(dd,:))];
%  end 

%calculate 0 to 1000 ms in ST in reference to the average baseline
 for ee  = 1: length(avgnegbase)
     STnegUSE = [STnegUSE 12*log2(fullnegative(:,ee)/avgnegbase(:,ee))];
     
 end
 
 for ff = 1: length(avgposbase)
        STposUSE = [STposUSE 12*log2(fullpositive(:,ff)/avgposbase(:,ff))]; 
 end
 
 %calculate average ST at each data point across all trials within a
 %subject that are included

 for hh = 1: length(STnegUSE)
     
     if any(oppNeg==hh) %for all the negative trials

             oppNegValue = [oppNegValue STnegUSE(:, hh)];
         
     elseif any (followNeg == hh)
         followNegValue = [followNegValue STnegUSE(:, hh)];
     else
     end
 end
 
  for jj = 1: length(STposUSE)
     
     if any(oppPos==jj) %for all the negative trials

             oppPosValue = [oppPosValue STposUSE(:, jj)];
         
     elseif any (followPos == jj)
         followPosValue = [followPosValue STposUSE(:, jj)];
     else
     end
 end
 
 for cc = 1: size(oppPosValue,1)
    if ~isempty(oppPosValue)
     oppPosUSE = [oppPosUSE; mean(oppPosValue(cc,:))];
    else
        oppPosUSE = [];
    
    end
    
     if ~isempty(oppNegValue)
    oppNegUSE = [oppNegUSE; mean(oppNegValue(cc,:))];
     else
         oppNegUSE = [];
     end
     
      if ~isempty(followPosValue)
    followPosUSE = [followPosUSE; mean(followPosValue(cc,:))];
      else
          followPosUSE = [];
      end
      
       if ~isempty(followNegValue)
    followNegUSE = [followNegUSE; mean(followNegValue(cc,:))];
       else
           followNegUSE = [];
       end
       
 end


oppPosUSE_downsamp = downsample(oppPosUSE, 10);
oppNegUSE_downsamp = downsample(oppNegUSE, 10);

followPosUSE_downsamp = downsample(followPosUSE, 10);
followNegUSE_downsamp= downsample(followNegUSE, 10);


%% saving into structure
%shift up (positive trials)
analysis.pos.baseline.full = concatmatrix_baseline_pos;
analysis.pos.baseline.avg = avgposbase;

analysis.pos.ST.baseline = STposbase;
analysis.pos.ST.reflexwindow  = STpos;

analysis.pos.reflex.full = concatmatrix_reflex_full_pos;
analysis.pos.reflex.window = concatmatrix_reflex_pos;

analysis.pos.baselinereflex = fullpositive;
analysis.pos.full.opp.avg =  oppPosUSE; 
analysis.pos.full.follow.avg =followPosUSE;
analysis.pos.full.opp.values =  oppPosValue; 
analysis.pos.full.follow.values =followPosValue;
analysis.pos.downsamp.opp.avg = oppPosUSE_downsamp;
analysis.pos.downsamp.follow.avg = followPosUSE_downsamp;

analysis.pos.idx.removebase = idx_maxabsSTposbase;
analysis.pos.idx.direction = directionpos;
analysis.pos.idx.opp = oppPos;
analysis.pos.idx.follow = followPos;

%shift down (negative trials)
analysis.neg.baseline.full = concatmatrix_baseline_neg;
analysis.neg.baseline.avg = avgnegbase;

analysis.neg.ST.baseline = STnegbase;
analysis.neg.ST.reflexwindow  = STneg;

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
analysis.neg.idx.direction = directionneg;
analysis.neg.idx.opp = oppNeg;
analysis.neg.idx.follow = followNeg;

save('analysis.mat', 'analysis')
