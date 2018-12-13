function LDPAR_Eclipse_Reflex_Analysis_nonaligned()
%for sustained reflex perturbations
    % Calculates delay times for each trial for mic and headphone stimuli
    % Calculates perturbation onset times for perturbed trials
    % Calculates average of 150-300 ms after perturbation begins

%% Script adapted from sustained reflex analysis script for pediatric perturbation (LHM) for PDGAP by Defne Abur 12/6/2017
    % 100 ms baseline period and 850 ms following onset
    % shift up and shift down 
    % for reflex, selection window is edited to be same as timing window so you can
    % see all possible times the reflex might begin 
    % Added the following to end of script, they can be repeated/re-run offline using script output of 'keepdata.mat'
        % takes 100 ms baseline period (as a vector)
        % converts 100 ms baseline period vector to cents relative to mean (single value) of baseline period
        % saves this as a normalized baseline vector EXCEPT:
               %trials with > 15 cent variability will be excluded
               %trials where there is pitch mistracking will be excluded
               %excluded trial nums will be saved and output at end 
        % takes 850 ms after onset of perturbation (as a vector) 
        % converts 850 ms perturbation period (as a vector) to cents relative to mean (single value) of baseline period
        % saves this as a normalized pert. vector
        % saves an average of the norm. 850ms vectors matrix as a single 850 ms that is a mean of all the perturbed trials
         
%% Script adapted from PDGAP for DFPAR by Hasini Weerathunge 05/10/2018
   % all edits done with comment (edit hasini 05/10/2018)
   % Edit 1: Time alignment of headphone signal with mic signal by
   % subtracting the delay between the mic and headphone signals, is
   % removed in identifying the reflex perturbation onset time in mic
   % signal ( common for all Pead reflex scripts as well as DFPAR/LDPAR )
   % Edit 2 : Time alignment of headphone signal with mic signal by
   % subtracting the delay between mic and headphone signals, is also
   % removed for all non-reflex and adaptation trials as well specifically
   % to be used in DFPAR analysis alone. As DFPAR and LDPAR data are
   % utlized/compared together LDPAR study analysis will also maintain this
   % policy. ( only for DFPAR/LDPAR study )
   % Edit 3 : editted the script to contain mic onset/ headphone onset and
   % common offset based on headphone onset and offset marked in
   % spectragram. Utlized delay to align the mic onset. Saving the arrays
   % of f0 contours from mic onset onwards for channel 1 for post hoc
   % analysis
%% Edit (1,2,3,4) Hasini Weerathunge 08/02/2018 - Pert error caused to to pert length not being 850ms - 
% editted script to fill the remained with -Inf 
% so the main arrays can still be concatenated without errors - files with
% such missinf daya will be discarded in the reflex-vector egenration
% scripts step
%% LDPAR folder pathways 

close all
clear all

%participant = input('Subject ID? ', 's'); %input subject ID number
prompt={'Subjet ID:','Amplification:'};
    name='Reflex Details';
    numlines=1;
    defaultanswer={'LDPAR','5'};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    participant = (answer{1});
    amplification = str2num(answer{2});

FolderName = uigetdir('R:\SteppLab2\SPEECHDATA\Voice\DFPAR_LDPAR','Select Data Folder'); %choose data folder; make sure there are "mic" and "perturbed_mic" subfolders in this main data folder
cd(FolderName) %navigate to folder
if ~exist('mic', 'dir')
mkdir mic %make folder to store mic files (ch1)
else
end

if ~exist('headphones', 'dir')
mkdir headphones %make folder to store headphone files (ch2)
else 
end

micFolder = [FolderName '\mic'];
headphonesFolder = [FolderName '\headphones'];

DestFolder = uigetdir('R:\SteppLab3\Projects\Voice\DFPAR_LDPAR\AuditoryPerturbationAnalysis\Analysis_Scripts\Data_temp','Select Destination Folder'); %choose folder to save data and figures in
cd(DestFolder)

if ~exist(participant, 'dir')
mkdir(participant) %make folder for participant
else 
end

PartFolder = [DestFolder '\' participant]; %set a path for the participant folder

cd(PartFolder)

%% defines which experiment analyzing % Hasini edit - %rename participant folder to appropriate script type e.g: Adapt_Shift_Up_5dB (folder naming)
choiceType = 'Reflex';

% choiceType = questdlg('Is this Adaptation or Reflex?', ...
%     'Script Type', ...
%     'Adaptation', 'Reflex','Reflex');

switch choiceType
    case 'Adaptation'
        scripttype = 1; %adaptation script
        
        experimentType = questdlg('Is this Shift-up or Shift-down?', ...
            'Experiment Type', ...
            'Shift-Up','Shift-Down','Control','Control');
        switch experimentType
            case  'Shift-Up'
                savetype = '_Shift-Up';
                Folder = ['Adapt_Shift_Up' '_' num2str(amplification) 'dB']; %rename participant folder to appropriate script type
                
                
                if ~exist(Folder, 'dir')
                    mkdir(Folder)  %shift up directory
                else
                end
                PartFolder = [PartFolder '\' Folder]; %rename participant folder to appropriate script type
                
            case  'Shift-Down'
                savetype = '_Shift-Down';
                Folder = ['Adapt_Shift_Down' '_' num2str(amplification) 'dB']; %rename participant folder to appropriate script type
               
                if ~exist(Folder, 'dir')
                    mkdir(Folder)  %shift down directory
                else
                end
                
                PartFolder = [PartFolder '\' Folder]; %rename participant folder to appropriate script type
            case  'Control'
                savetype = '_Control';
                Folder = ['Adapt_Control' '_' num2str(amplification) 'dB']; %rename participant folder to appropriate script type
               
                if ~exist(Folder, 'dir')
                    mkdir(Folder)  %control directory
                else
                end
                
                PartFolder = [PartFolder '\' Folder]; %rename participant folder to appropriate script type
                
        end
        
    case 'Reflex'
        scripttype = 2; %reflex script
        experimentType ='Shift-up'; % Hasini edit - removed choice for now as we do shift up only
%         experimentType = questdlg('Is this Shift-up or Shift-down?', ...
%             'Experiment Type', ...
%             'Shift-up', 'Shift-down','Shift-up');
        switch experimentType
            case  'Shift-up'
                savetype = 'Reflex_Shift_up';
                Folder = ['Reflex_Shift_up' '_' num2str(amplification) 'dB']; %rename participant folder to appropriate script type
               
                if ~exist(Folder, 'dir')
                    mkdir(Folder)  
                else
                end
                
                PartFolder = [PartFolder '\' Folder]; %rename participant folder to appropriate script type
            case  'Shift-down'
                savetype = 'Reflex_Shift_down';
                Folder = ['Reflex_Shift_down' '_' num2str(amplification) 'dB']; %rename participant folder to appropriate script type
               
                if ~exist(Folder, 'dir')
                    mkdir(Folder)
                else
                end
                
                PartFolder = [PartFolder '\' Folder]; %rename participant folder to appropriate script type
           case  'reflexC'
                savetype = '_reflexC';
                
                if ~exist('reflexC', 'dir')
                    mkdir reflexC 
                else
                end
                
                PartFolder = [PartFolder '\' 'reflexC']; %rename participant folder to appropriate script type
        end
end

savefilename = strcat('dataOutput', savetype, '.mat');



cd(PartFolder)
if ~exist('mic', 'dir')
mkdir mic %make folder to store mic files (ch1)
else
end

if ~exist('headphones', 'dir')
mkdir headphones %make folder to store headphone files (ch2)
else 
end

txtmicFolder = [PartFolder '\mic'];
txtheadphonesFolder = [PartFolder '\headphones'];


%Find total number of wav files`
d = dir([FolderName, '\Trial*.wav']); %set the directory
fnames = sort_nat({d.name}); %sort directory alphabetically in order


%% splits and saves ch1 (mic signal) and ch2 (headphones) output for all files
ch1filenames_save = [];
ch2filenames_save = [];

splitwavfiles = input('Do you need to split wav files? (1=yes, 0=no) ', 's');  %type yes if you need to run through splitting wav files if this is your first time analyzing this subject
splitwavfiles = str2num(splitwavfiles);
if splitwavfiles ==1
    for j = 1:length(fnames) %creates loop, for 1 until the last file
        filename = strcat(FolderName,'\',fnames{j}); %creates full file name using path and wav name
        [Y, Fs] = audioread(filename); %loads wav file
        ch1filename = strcat(FolderName,'\mic','\',fnames{j}); % makes new file path to "mic" folder
        ch2filename = strcat(FolderName,'\headphones','\',fnames{j}); % makes new file path to "mic" folder
        
        
        micdata  = Y(:,1);
        headphonedata = Y(:,2);
        
        
        %% uncomment to split audio data
        audiowrite(ch1filename,micdata,Fs) %creates audio file
        audiowrite(ch2filename,headphonedata,Fs) %creates audio file
        
        ch1filenames_save = [ch1filenames_save; {ch1filename}];
        ch2filenames_save = [ch2filenames_save; {ch2filename}];
        
    end
else
    %do nothing
end


%go back to destination folder and go one up to all analysis scripts. If the
%analysis folder is changed then THIS needs to be changed to navigate to
%the correct folder for the scripts
cd(DestFolder); %go to destination folder
cd ../ %go up one to folder where praat interface is

praatfiles = input('Do you need to read in/analyze the praat files to create a txt file? (1=yes, 0=no) ', 's');  %type yes if you need to run through the praat script if this is your first time analyzing this subject
praatfiles = str2num(praatfiles);

voicingOnset = input('voicing onset? (default is .02):');
if isempty(voicingOnset)
    voicingOnset = 0.02;
end

if praatfiles ==1
    %% make sure you have the scripts you need for praat interface to work
    pi_dir = 'praat_interface' ;
    if ~exist( pi_dir , 'dir' )
        error( 'folder ''praat_interface'' not found in directory' )
    end
    p_fn = fullfile( pi_dir , 'praat.exe' ) ;
    if ~exist( p_fn , 'file' )
        error( 'file ''praat.exe'' not found' )
    end
    sp_fn = fullfile( pi_dir , 'sendpraat.exe' ) ;
    if ~exist( sp_fn , 'file' )
        error( 'file ''sendpraat.exe'' not found' )
    end
    gp_fn = fullfile( pi_dir , 'get_f0_EET.praat' ) ;
    if ~exist( gp_fn , 'file' )
        error( 'file ''get_f0_EET.praat'' not found' )
    end
    
    
    %% initilization settings for praat
    
    
    voicingThreshold = input('voicing threshold for praat? (default is 0.45): ');
    if isempty(voicingThreshold)
        voicingThreshold = 0.45;
    end
    pitchFloor = input('pitch floor for praat? (default is 75 Hz):');
    if isempty(pitchFloor)
        pitchFloor = 75;
    end
    pitchCeiling = input('pitch ceiling for praat? (default is 500 Hz):');
    if isempty(pitchCeiling)
        pitchCeiling = 500;
    end
    
    maxformant = input('max formant Hz for praat? (default is 5500 Hz):');
    if isempty(maxformant)
        maxformant = 5500;
    end
    
    numformant = input('number of formants for praat? (default is 5):');
    if isempty(numformant)
        numformant = 5;
    end
    
    
    newmic = [micFolder, '\']; %add a slash to the mic folder
    ext = '.wav'; %extension of files
    

    %Build DOS calls to control praat
    call2 = sprintf( '%s praat "execute %s %s %s %s %f %f %f %f %f"' , ...
        fullfile( pwd , sp_fn ) , ... %location of sendpraat.exe
        fullfile( pwd , gp_fn ) , ... %location of praat script wrote (get_fo_EET.praat)
        newmic , ... %where wav files are located
        ext , ... %extension of files
        txtmicFolder, ... %destination where .txt files will stay
        voicingThreshold,  ... %voicing threshold indicated
        pitchFloor, ... %pitch floor indicated
        pitchCeiling, ... %pitch ceiling indicated
        maxformant, ... %max formant value in hx
        numformant ... %number of formants to calculate
        ) ;
    
    
    
    [ s , r ] = dos( call2 ) ;
    if s ~= 0
        dos( [ fullfile( pwd , p_fn ) ' &' ] ) ;
        [ s , r ] = dos( call2 ) ;
        if s ~= 0
            disp( r )
            error( [ 'ERROR: get_f0_EET.praat failed to run. Make sure '...
                'there are no spaces in the selected filepaths.'])
        end
    end
    
    
    %% headphones
    newheadphones = [headphonesFolder, '\']; %add a slash to the headphones folder
    ext = '.wav'; %extension of files
    
    %Build DOS calls to control praat
    call3 = sprintf( '%s praat "execute %s %s %s %s %f %f %f %f %f"' , ...
        fullfile( pwd , sp_fn ) , ... %location of sendpraat.exe
        fullfile( pwd , gp_fn ) , ... %location of praat script wrote (get_fo_EET.praat)
        newheadphones , ... %where wav files are located
        ext , ... %extension of files
        txtheadphonesFolder, ... %destination where .txt files will stay
        voicingThreshold,  ... %voicing threshold indicated
        pitchFloor, ... %pitch floor indicated
        pitchCeiling, ... %pitch ceiling indicated
        maxformant, ... %max formant value in hx
        numformant ... %number of formants to calculate
        ) ;
    
    %strrep(fnames(:,1),'.','_'), ...
    
    [ s , r ] = dos( call3 ) ;
    if s ~= 0
        dos( [ fullfile( pwd , p_fn ) ' &' ] ) ;
        [ s , r ] = dos( call3 ) ;
        if s ~= 0
            disp( r )
            error( [ 'ERROR: get_f0_EET.praat failed to run. Make sure '...
                'there are no spaces in the selected filepaths.'])
        end
    end
else
    %do nothing
end

%set up for reading in .txt files
field1 = 'time';
field2 = 'F0';
field3 = 'Int';
field4 = 'F1';
field5 = 'F2';
field6 = 'F3';


cd(txtmicFolder); %start with the mic file to get the list of names
%create a strucutre with all the praat data
txtlist = dir(fullfile(pwd, '*.txt'));
txtnames = {txtlist.name}';
txtnames = sort_nat(txtnames);
txtFileNames = char(strrep(txtnames, '.txt', ''));
micFile =1; %start with the mic trial

%set up structure for saving
%keepdata.keepcode = [{'trial'},{'1 if keep was checked'}, {'1 if redo was checked'}, {'reason why redo was checked'} ]; %set up for structure
%keepdata.average = [{'trial'}, {'1 if keep was checked'}, {'1 if redo was checked'}, {'reason why redo was checked'}, {'length of shifted reflex: steady portion(shifted reflex only)'},{'length of FULL shifted reflex'}, {'delay'}, {'length from voice onset detection to start of reflex'}, {'voice onset time (shifted reflex only)'},  {'avg Int(ch1)'}, {'avg Int(ch2)'},{'avg F0(ch1)'}, {'avg F0(ch2)'},{'avg F1(ch1)'}, {'avg F1(ch2)'},{'avg F2(ch1)'}, {'avg F2(ch2)'}];
%keepdata.average = [{'trial'}, {'1 if keep was checked'}, {'1 if redo was checked'}, {'reason why redo was checked'}, {'delay'}, {'length from voice onset detection to start of reflex'}, {'voice onset time (shifted reflex only)'},  {'avg Int(ch1)'}, {'avg Int(ch2)'},{'avg F0 basleine'},{'avg F0(ch1)'}, {'avg F0(ch2)'},{'avg F1(ch1)'}, {'avg F1(ch2)'},{'avg F2(ch1)'}, {'avg F2(ch2)'}];
keepdata.average = [{'trial'}, {'1 if keep was checked'}, {'1 if redo was checked'}, {'reason why redo was checked'}, {'delay'}, {'length from voice onset detection to start of reflex'}, {'voice onset time (shifted reflex only)'},{'voice onset time Ch1 (control trials only)'}, {'voice onset time Ch2 (control trials only)'},{'voice offset time (control trials only)'},  {'avg Int(ch1)'}, {'avg Int(ch2)'},{'avg F0 basleine'},{'avg F0(ch1)'}, {'avg F0(ch2)'},{'avg F1(ch1)'}, {'avg F1(ch2)'},{'avg F2(ch1)'}, {'avg F2(ch2)'}];% Edit 3 Hasini : new variables to save
   
cd(PartFolder)
if exist( savefilename )
load(savefilename)
currenttrial = keepdata.average(end,1); %read in last file name
currenttrial = strsplit(currenttrial{1,1},'_'); %split name by underscores
kk = str2num(currenttrial{2})+1; %trial # is second part (assumes file names: trial_##)
else
kk =1;
end



for fileNum = kk: length(txtnames)
    
    close all %close open figures
    
    
    %move to the folder where the wav file is stored
    cd(FolderName) %go to the folder with the wav files
    list = dir(fullfile(pwd, 'trial*')); %list the information for all the trial
    trialnames = {list.name}'; %pull out the trial names
    trialnames = sort_nat(trialnames); %puts them in the same order as the txt file names
    
    
    %pull in the wav file
    FiletoPlay = trialnames{fileNum};
    [Y, Fs] = audioread(FiletoPlay);
    ch1_wav = Y(:,1);
    ch2_wav = Y(:,2);
    
    cd(txtmicFolder); %need to repeat this for the loop
    
    %pull in an individual file and place it in a table
    T = readtable(char(txtnames(fileNum,:)),'Delimiter',' ');
    T.Properties.VariableNames = {'Time_s' 'F0_Hz' 'Intensity_dB' 'F1_Hz' 'F2_Hz' 'F3_Hz'};
    
    %pull out 5 measures so can convert them to numbers
    time_ch1= T.Time_s;
    F0_ch1 = str2num(char(strrep(T.F0_Hz, '--undefined--', '0')));
    Int_ch1 = str2num(char(strrep(T.Intensity_dB, '--undefined--', '0')));
    F1_ch1 = str2num(char(strrep(T.F1_Hz, '--undefined--', '0')));
    F2_ch1 = str2num(char(strrep(T.F2_Hz, '--undefined--', '0')));
    F3_ch1 = str2num(char(strrep(T.F3_Hz, '--undefined--', '0')));
    
    
    trialname_update = strrep(txtnames, '-', 'NEG'); %needed because tables won't read in negative label
    trialname_update  = strrep(trialname_update, '.txt', ' ');
    trialnum = char(strtrim(trialname_update(fileNum,:)));
    data.mic.(trialnum) = struct(field1,time_ch1,field2,F0_ch1,field3,Int_ch1,field4,F1_ch1,field5,F2_ch1,field6,F3_ch1);
    
    timeVec=[0:1/44100:(length(Y)-1)/44100]'; % create a time vector
    
    %will use ch 1 intensity to set window for ch 1 and ch 2
    Intrepeat = round(length(timeVec)/length(Int_ch1));
    newInt_ch1 = repelem(Int_ch1, Intrepeat);
    if length(newInt_ch1)> length(timeVec)
        Intrepeat = Intrepeat - 1;
        newInt_ch1 = repelem(Int_ch1, Intrepeat);
    end
    zervec = zeros(length(timeVec)-length(newInt_ch1), 1);    %have to make the sizes match the length of the spectrogram window
    Int_use_ch1 = [newInt_ch1; zervec]; %have to add some zeros at the end to make the files the same length
    
    F0repeat = round(length(timeVec)/length(F0_ch1));
    newF0_ch1 = repelem(F0_ch1,F0repeat );
    if length(newF0_ch1)> length(timeVec)
        F0repeat = F0repeat - 1;
        newF0_ch1 = repelem(F0_ch1,F0repeat );
    end
    F0_use_ch1 = [newF0_ch1; zervec];
    
    
    F1repeat = round(length(timeVec)/length(F1_ch1));
    newF1_ch1 = repelem(F1_ch1,F1repeat );
    if length(newF1_ch1)> length(timeVec)
        F1repeat = F1repeat - 1;
        newF1_ch1 = repelem(F1_ch1,F1repeat );
    end
    F1_use_ch1 = [newF1_ch1; zervec];
    F1_use_ch1 = F1_use_ch1/1000; %needed to fit in the scale of the spectrogram
    
    F2repeat = round(length(timeVec)/length(F2_ch1));
    newF2_ch1 = repelem(F2_ch1, F2repeat);
    if length(newF2_ch1)> length(timeVec)
        F2repeat = F2repeat - 1;
        newF2_ch1 = repelem(F2_ch1, F2repeat);
    end
    F2_use_ch1 = [newF2_ch1; zervec];
    F2_use_ch1 = F2_use_ch1/1000; %needed to fit in the scale of the spectrogram
    
    F3repeat = round(length(timeVec)/length(F3_ch1));
    newF3_ch1 = repelem(F3_ch1, F3repeat);
    if length(newF3_ch1)> length(timeVec)
        F3repeat = F3repeat - 1;
        newF3_ch1 = repelem(F3_ch1, F3repeat);
    end
    F3_use_ch1 = [newF3_ch1; zervec];
    F3_use_ch1 = F3_use_ch1/1000; %needed to fit in the scale of the spectrogram
    
    
    %plot the intensity figure so can deteremine where voicing starts
    intFigure = figure();
    set(intFigure,'units','normalized','outerposition',[0 0 1 1]) %makes figure full screen
    plot(timeVec, Int_use_ch1)
    [intx, inty] = ginput(1);
    
    startVoicing = intx(1);
    [c intStartIndex_Voicing] = min(abs(timeVec - startVoicing));
    plotTime = timeVec(intStartIndex_Voicing:end);
    plotF0_ch1 = F0_use_ch1(intStartIndex_Voicing:end);
    plotF1_ch1 = F1_use_ch1(intStartIndex_Voicing:end);
    plotF2_ch1 = F2_use_ch1(intStartIndex_Voicing:end);
    PlotF3_ch1 = F3_use_ch1(intStartIndex_Voicing:end);
    plotWav_ch1 = ch1_wav(intStartIndex_Voicing:end);
    plotWav_ch2 = ch2_wav(intStartIndex_Voicing:end);
    
    close(intFigure)
    
    %make the analysis figure
    mFigure = figure();
    set(mFigure,'units','normalized','outerposition',[0 0 1 1]) %makes figure full screen
    
    %% setup textboxes
    
    delaytextbx_name = annotation('textbox');
    delaytextbx_name.Position = [.78 .81 .12 .025];
    delaytextbx_name.String = 'Delay between channels (ms)';
    
    delaytextbx = annotation('textbox');
    delaytextbx.Position = [.9 .81 .05 .025];
    
    F0textbx_name = annotation('textbox');
    F0textbx_name.Position = [.52 .1 .1 .025];
    F0textbx_name.String = 'Average F0 (Hz)';
    
    F0textbx_ch1 = annotation('textbox');
    F0textbx_ch1.Position = [.62 .1 .1 .025];
    
    F0textbx_ch2 = annotation('textbox');
    F0textbx_ch2.Position = [.72 .1 .1 .025];
    
    F1textbx_name = annotation('textbox');
    F1textbx_name.Position = [.52 .07 .1 .025];
    F1textbx_name.String = 'Average F1 (Hz)';
    
    F1textbx_ch1 = annotation('textbox');
    F1textbx_ch1.Position = [.62 .07 .1 .025];
    
    F1textbx_ch2 = annotation('textbox');
    F1textbx_ch2.Position = [.72 .07 .1 .025];
    
    F2textbx_name = annotation('textbox');
    F2textbx_name.Position = [.52 .04 .1 .025];
    F2textbx_name.String = 'Average F2 (Hz)';
    
    F2textbx_ch1 = annotation('textbox');
    F2textbx_ch1.Position = [.62 .04 .1 .025];
    
    F2textbx_ch2 = annotation('textbox');
    F2textbx_ch2.Position = [.72 .04 .1 .025];
    
    %set up checkboxes
    keeepcheckbox = uicontrol('style','checkbox','String', 'keep?', 'units', 'normalized', 'Position', [.85 .12 .08 .02]);
    redocheckbox = uicontrol('style','checkbox','String', 'redo offline?', 'units', 'normalized', 'Position', [.85 .09 .08 .02]);
    
    %set up field to fill in why you clicked redo
    whyredobox = uicontrol('style','edit','String', '0', 'units', 'normalized', 'Position', [.85 .04 .1 .05]);
    
    %set of variables to be passed between nested functions
    selectTime_ch1 = [];
    selectF0_ch1 = [];
    norm_selectF0_ch1 =[];

    selectF1_ch1 = [];
    selectF2_ch1 = [];
    selectInt_ch1= [];
    selectTime_ch2 = [];
    selectF0_ch2 = [];
    selectF1_ch2 = [];
    selectF2_ch2 = [];
    selectInt_ch2= [];
     F0calc_ch1 = [];
    F1calc_ch1 = [];
    F2calc_ch1 = [];
    Intcalc_ch1= [];
    F0calc_ch2= [];
    F1calc_ch2= [];
    F2calc_ch2= [];
    Intcalc_ch2= [];
    timex_delay =[];
    timey_delay = [];
    delaycalc = [];
    delaycalcsec = [];
    keepvalue =[];
    redovalue = [];
    whyredo = [];
    F01 = [];
    F02 = [];
    F11 = [];
    F21 = [];
    F12 = [];
    F22 = [];
    reflexShift = [];
    
    %for reflex only
     selectTime_ch1Mean = [];
    selectF0_ch1Mean = [];
    selectF1_ch1Mean = [];
    selectF2_ch1Mean = [];
    selectInt_ch1Mean= [];
    selectTime_ch2Mean = [];
    selectF0_ch2Mean = [];
    selectF1_ch2Mean = [];
    selectF2_ch2Mean = [];
    selectInt_ch2Mean= [];
    baselineTime_ch1 = [];
    baselineF0_ch1 = [];
    
    norm_baselineF0_ch1 =[]; 
    
    baselineF1_ch1 = [];
    baselineF2_ch1 = [];
    baselineInt_ch1= [];
    baselineTime_ch2 = [];
    baselineF0_ch2 = [];
    baselineF1_ch2 = [];
    baselineF2_ch2 = [];
    baselineInt_ch2= [];
     timeVoiceOnset = [];
    lengthreflex = []; %not used in sustained
    lengthVoiceOnset_Reflex = [];
    lengthfullreflex = []; %not used in sustained
     baselineF0_ch1_mean = [];
     baselineF0_ch2_mean = [];
    
   % reflexTime_ch1 = [];
   % reflexF0_ch1 = [];
   % reflexF1_ch1 = [];
   % reflexF2_ch1 = [];
   % reflexInt_ch1= [];
   % reflexTime_ch2 = [];
   % reflexF0_ch2 = [];
   % reflexF1_ch2 = [];
   % reflexF2_ch2 = [];
    %reflexInt_ch2= [];
    timeVoiceOnset_ch1 = [];
    timeVoiceOffset=[];
    timeVoiceOnset_ch2=[];
   
    
    %% set up pushbuttons to show/hide formants and pitch
    PushButtonFormant_ch1 = uicontrol(gcf, 'Style', 'push', 'String', 'show formants (ch1)',  'Callback',@showformants_ch1, 'units', 'normalized','Position' , [.8 .66 .1 .025] );
    PushButtonPitch_ch1 = uicontrol(gcf, 'Style', 'push', 'String', 'show pitch (ch1)',  'Callback',@showpitch_ch1, 'units', 'normalized','Position' , [.8 .63 .1 .025] );
    PushButtonFormant_ch2 = uicontrol(gcf, 'Style', 'push', 'String', 'show formants (ch2)', 'Callback',@showformants_ch2, 'units', 'normalized','Position' , [.8 .39 .1 .025] );
    PushButtonPitch_ch2 = uicontrol(gcf, 'Style', 'push', 'String', 'show pitch (ch2)',  'Callback',@showpitch_ch2, 'units', 'normalized','Position' , [.8 .36 .1 .025] );
    
    
    PushButtonFormantHide_ch1 = uicontrol(gcf, 'Style', 'push', 'String', 'hide formants (ch1)', 'Callback',@hideformants_ch1, 'units', 'normalized','Position' , [.9 .66 .1 .025] );
    PushButtonPitchHide_ch1 = uicontrol(gcf, 'Style', 'push', 'String', 'hide pitch (ch1)', 'Callback',@hidepitch_ch1, 'units', 'normalized','Position' , [.9 .63 .1 .025] );
    PushButtonFormantHide_ch2 = uicontrol(gcf, 'Style', 'push', 'String', 'hide formants (ch2)',  'Callback',@hideformants_ch2, 'units', 'normalized','Position' , [.9 .39 .1 .025] );
    PushButtonPitchHide_ch2 = uicontrol(gcf, 'Style', 'push', 'String', 'hide pitch (ch2)', 'Callback',@hidepitch_ch2, 'units', 'normalized','Position' , [.9 .36 .1 .025] );
    
    %set up push buttons for measuring
    PushButtonspect = uicontrol(gcf, 'Style', 'push', 'String', 'Spectrogram Timing: select in Ch 2', 'Visible', 'Off', 'Callback',@callbackspect, 'units', 'normalized', 'Position' , [.3 .04 .2 .05] );
    PushButtondelay = uicontrol(gcf, 'Style', 'push', 'String', 'Press when ready to measure delay (zoom in first- start with ch1)', 'Callback',@callbackdelay, 'units', 'normalized', 'Position' , [.02 .04 .2 .05]  );
    PushButtonnext = uicontrol(gcf, 'Style', 'push', 'String', 'Next Trial', 'BackgroundColor', [1 0 0], 'Callback',@callbacknext, 'units', 'normalized','Position' , [.95 .04 .05 .05] );
    
    %set up push buttons for playing sound
    PushButtonplay_ch1 = uicontrol(gcf, 'Style', 'push', 'String', 'play ch1', 'BackgroundColor', [0 1 0],  'Callback',@playpitch_ch1, 'units', 'normalized','Position' , [.15 .95 .1 .025] );
    PushButtonplay_ch2 = uicontrol(gcf, 'Style', 'push', 'String', 'play ch2', 'BackgroundColor', [0 1 0], 'Callback',@playpitch_ch2, 'units', 'normalized','Position' , [.15 .79 .1 .025] );
    
    
    %% make subplots
    
    %make subplot for channel 1 waveform
    subplot(4,1,1)
    plot(plotTime, plotWav_ch1);
    axis tight
    xlabel('time(s)') % x-axis label
    ylabel('amplitude') % y-axis label
    sub1 = subplot(4,1,1); %pull in information about current subplot
    pnew_sub1 = [.02 .85 .9 .1]; %make a new position
    set(sub1, 'Position', pnew_sub1)
    set(sub1, 'DefaulttextInterpreter', 'none')
    title (sprintf ('channel1: microphone %s', trialnum))
    hold on
    
    
    %make subplot for channel 2 waveform
    subplot(4,1,2)
    plot(plotTime, plotWav_ch2);
    axis tight
    title ('channel2 :headphones')
    xlabel('time(s)') % x-axis label
    ylabel('amplitude') % y-axis label
    sub2 = subplot(4,1,2); %pull in information about current subplot
    pnew_sub2 = [.035 .69 .9 .1]; %make a new position
    set(sub2, 'Position', pnew_sub2)
    hold on
    
    
    %spectrogram settings
    win_len = round(0.007 * Fs); %this calculates the number of frames in analysis window (0.005 s window length in Praat wideband spectrogram settings - changed to .007)
    win_overlap = round(0.006*Fs); %this is the amount of window overlap for analysis (in number of samples - multiply seconds by sampling frequency);
    nfft = 1024; %uses nfft sampling points to calculate the discrete Fourier transform,  higher number will clarify glottal pulses( factor of 2 increases)
    
    subplot(4,1,3)
    spectrogram(ch1_wav, win_len, win_overlap, nfft, Fs, 'yaxis');
    starttime_spec = timeVec(intStartIndex_Voicing); %based on intensity calculations above, this is the start time we want in the window
    sub3 = subplot(4,1,3); %pull in information about current subplot
    pnew = [.038 .42 .9 .2]; %make a new position
    set(sub3, 'Position', pnew)
    ylimsub3 = get(sub3, 'Ylim'); %get current y limit
    newylim3 = [ylimsub3(1) 5]; %set max to 5000 Hz, keep min same
    set(sub3, 'Ylim', newylim3)
    xlimsub3 = get(sub3, 'Xlim'); %get current x limit
    newxlim3 = [starttime_spec xlimsub3(2)]; %set length of window to be the same as what you will plot over it
  
    set(sub3, 'Xlim', newxlim3)
    colormap('gray') %sets and customizes color map
    map = colormap;
    colormap(1-map)
    yyaxis right
    ylim([60 400])
    hold on
    
    
    %% switch to the files where the headphones data are saved
    cd(txtheadphonesFolder); %need to repeat this for the loop
    
    %pull in an individual file and place it in a table
    T = readtable(char(txtnames(fileNum,:)),'Delimiter',' ');
    T.Properties.VariableNames = {'Time_s' 'F0_Hz' 'Intensity_dB' 'F1_Hz' 'F2_Hz' 'F3_Hz'};
    
    %pull out 5 measures so can convert them to numbers
    time_ch2= T.Time_s;
    F0_ch2 = str2num(char(strrep(T.F0_Hz, '--undefined--', '0')));
    Int_ch2 = str2num(char(strrep(T.Intensity_dB, '--undefined--', '0')));
    F1_ch2 = str2num(char(strrep(T.F1_Hz, '--undefined--', '0')));
    F2_ch2 = str2num(char(strrep(T.F2_Hz, '--undefined--', '0')));
    F3_ch2 = str2num(char(strrep(T.F3_Hz, '--undefined--', '0')));
    
 
    trialnum = char(strtrim(trialname_update(fileNum,:)));
    data.headphones.(trialnum) = struct(field1,time_ch2,field2,F0_ch2,field3,Int_ch2,field4,F1_ch2,field5,F2_ch2,field6,F3_ch2);
    
    %int calc is just used for saving later
    Intrepeat = round(length(timeVec)/length(Int_ch2));
    newInt_ch2 = repelem(Int_ch2, Intrepeat);
    if length(newInt_ch2)> length(timeVec)
        Intrepeat = Intrepeat - 1;
        newInt_ch2 = repelem(Int_ch2, Intrepeat);
    end
    zervec = zeros(length(timeVec)-length(newInt_ch2), 1);    %have to make the sizes match the length of the spectrogram window
    Int_use_ch2 = [newInt_ch2; zervec]; %have to add some zeros at the end to make the files the same length
    
    F0repeat = round(length(timeVec)/length(F0_ch2));
    newF0_ch2 = repelem(F0_ch2,F0repeat );
    if length(newF0_ch2)> length(timeVec)
        F0repeat = F0repeat - 1;
        newF0_ch2 = repelem(F0_ch2,F0repeat );
    end
    F0_use_ch2 = [newF0_ch2; zervec];
    
    
    F1repeat = round(length(timeVec)/length(F1_ch2));
    newF1_ch2 = repelem(F1_ch2,F1repeat );
    if length(newF1_ch2)> length(timeVec)
        F1repeat = F1repeat - 1;
        newF1_ch2 = repelem(F1_ch2,F1repeat );
    end
    F1_use_ch2 = [newF1_ch2; zervec];
    F1_use_ch2 = F1_use_ch2/1000; %needed to fit in the scale of the spectrogram
    
    F2repeat = round(length(timeVec)/length(F2_ch2));
    newF2_ch2 = repelem(F2_ch2, F2repeat);
    if length(newF2_ch2)> length(timeVec)
        F2repeat = F2repeat - 1;
        newF2_ch2 = repelem(F2_ch2, F2repeat);
    end
    F2_use_ch2 = [newF2_ch2; zervec];
    F2_use_ch2 = F2_use_ch2/1000; %needed to fit in the scale of the spectrogram
    
    F3repeat = round(length(timeVec)/length(F3_ch2));
    newF3_ch2 = repelem(F3_ch2, F3repeat);
    if length(newF3_ch2)> length(timeVec)
        F3repeat = F3repeat - 1;
        newF3_ch2 = repelem(F3_ch2, F3repeat);
    end
    F3_use_ch2 = [newF3_ch2; zervec];
    F3_use_ch2 = F3_use_ch2/1000; %needed to fit in the scale of the spectrogram
    
    
    %use intensity start from channel1
    plotTime = timeVec(intStartIndex_Voicing:end);
    plotF0_ch2 = F0_use_ch2(intStartIndex_Voicing:end);
    plotF1_ch2 = F1_use_ch2(intStartIndex_Voicing:end);
    plotF2_ch2 = F2_use_ch2(intStartIndex_Voicing:end);
    PlotF3_ch2 = F3_use_ch2(intStartIndex_Voicing:end);
    
    
    %% spectrogram of channel2
    subplot(4,1,4)
    spectrogram(Y(:,2), win_len, win_overlap, nfft, Fs, 'yaxis');
    sub4 = subplot(4,1,4); %pull in information about current subplot
    pnew_sub4 = [.038 .16 .9 .2]; %make a new position
    set(sub4, 'Position', pnew_sub4)
    ylimsub4 = get(sub4, 'Ylim'); %get current y limit
    newylim4 = [ylimsub4(1) 5]; %set max to 5000 Hz, keep min same
    set(sub4, 'Ylim', newylim4)
    xlimsub4 = get(sub4, 'Xlim'); %get current x limit
    newxlim4 = [starttime_spec xlimsub4(2)]; %set length of window to be the same as what you will plot over it
    set(sub4, 'Xlim', newxlim4)
    colormap('gray') %sets and customizes color map
    map = colormap;
    colormap(1-map)
    yyaxis right
    ylim([60 400])
    hold on
    
    %% if a reflex script, need to calculate/save other information
    if scripttype == 2
        
        Reflexcond = questdlg('Is a shifted Reflex trial?', ...
            'Reflex Type', ...
            'Yes','No','Yes');
        
        switch Reflexcond
            case 'Yes'
                reflexShift =1;
                axes(sub1)
                
                I = find(ch1_wav>voicingOnset,1,'first');
                timeVoiceOnset = I/44100;
                timeVoiceOnset_ch1= 0; % Edit 3 Hasini : start saving onsets -none for reflex
                timeVoiceOnset_ch2 = 0; % Edit 3 Hasini : start saving onsets -none for reflex
                timeVoiceOffset = 0; % Edit 3 Hasini : start saving offsets -none for reflex
        
                
                axes(sub3)
                newxlim3_reflex = [starttime_spec starttime_spec*6]; %set length of full window so you can zoom in to where reflex is
                %newxlim3_reflex = [.7 1.7]; %zoom in to .5 to 1.5 so can better see reflex
                set(sub3, 'Xlim', newxlim3_reflex)
                yyaxis right
                F01 = plot(sub3, timeVec,F0_use_ch1, 'r', 'LineWidth', 1.5);
                hold on
                
                
                axes(sub4)
                newxlim4_reflex = [starttime_spec starttime_spec*6]; %set length of full window so you can zoom in to where reflex is
              
                %newxlim4_reflex = [.7 1.7]; %zoom in to .5 to 1.5 so can better see reflex
                set(sub4, 'Xlim', newxlim4_reflex)
                yyaxis right
                F02 = plot(sub4, timeVec,F0_use_ch2, 'r', 'LineWidth', 1.5);
                hold on
            case 'No'
                reflexShift =2;
                timeVoiceOnset = 0; %if it is not a shifted reflex trial
                
        end
              
    else %if it is an adaptation trial
        reflexShift =0;
        timeVoiceOnset = 0;
        
        
    end
    
    waitfor(mFigure)
%% Save data

    %keepdata.average = [ keepdata.average; {trialnum}, keepvalue, redovalue,whyredo, delaycalc, lengthVoiceOnset_Reflex,  timeVoiceOnset, Intcalc_ch1, Intcalc_ch2, baselineF0_ch1_mean, F0calc_ch1, F0calc_ch2, F1calc_ch1, F1calc_ch2, F2calc_ch1, F2calc_ch2];
    keepdata.average = [ keepdata.average; {trialnum}, keepvalue, redovalue,whyredo,delaycalc, lengthVoiceOnset_Reflex,  timeVoiceOnset,timeVoiceOnset_ch1,timeVoiceOnset_ch2,timeVoiceOffset, Intcalc_ch1, Intcalc_ch2,  baselineF0_ch1_mean, F0calc_ch1, F0calc_ch2, F1calc_ch1, F1calc_ch2, F2calc_ch1, F2calc_ch2]; % Edit 3 Hasini : new variables to save
    keepdata.All_Ch1data.(trialnum) = [selectTime_ch1, selectInt_ch1,selectF0_ch1, selectF1_ch1*1000,selectF2_ch1*1000];
    keepdata.All_Ch2data.(trialnum) = [selectTime_ch2, selectInt_ch2,selectF0_ch2,selectF1_ch2*1000,selectF2_ch2*1000];
          
    if   reflexShift == 1
        keepdata.F0_Ch1_Baseline.(trialnum) = [baselineTime_ch1, baselineInt_ch1,baselineF0_ch1, norm_baselineF0_ch1]; %baseline time, intensity, F0 values (Hz), F0 values (normalized cents relative to mean of baseline period for that trial)
        keepdata.Formants_Ch1_Baseline.(trialnum) = [baselineTime_ch1, baselineInt_ch1,baselineF1_ch1*1000,baselineF2_ch1*1000]; %baseline time, intensity, F1 values (Hz), F2 values (Hz)
        keepdata.F0_Ch1_Pert.(trialnum) =[selectTime_ch1, selectInt_ch1, selectF0_ch1, norm_selectF0_ch1]; %time in sec (where 0 is onset of perturbation), ch1 intensity, F0 values (Hz), F0 values (normalized cents relative to mean baseline period for that trial)

    %% Make vectors for plotting
        % save vectors of Ch1 baseline period and time after perturbation for plotting
        % in Hz and in normalized cents (both relative to mean of baseline period)
        
        % currently not saving all the vectors for calculations, but can
        % always run those again after keepdata has been generated
        % previously by uncommenting lines where it says keepdata. and
        % commenting the lines below without keepdata (ie.861 - 864)

        %% baseline time vector save sampling every 10 ms
        
        %100 ms before perturbation onset (10 samples away because sample every 10 ms), based on Behroozmand, Korzyukov, Sattler, Larson 2012 
   
        
        kdbt = keepdata.F0_Ch1_Baseline;
        FN = fieldnames(kdbt);
        baseline_time = [];

        for trial = 1:numel(FN)
            trialName = FN{trial};
            kdbt_Temp = kdbt.(trialName);
           baseline_time = [baseline_time, kdbt_Temp(1:440:end,1)]; %takes time every 10 ms

        end
        
        baseline_time_Avg = mean(baseline_time, 2);
        keepdata.Baseline_times = baseline_time;
        
        clear kdbt
        clear FN
        clear kdbt_Temp
        clear trial
        clear trialName
        
        %% baseline 100 ms ch1 ALL TRIALS
        kdbc1d = keepdata.F0_Ch1_Baseline;
        FN = fieldnames(kdbc1d);
        baseline_F0_ch1 = [];

        for trial = 1:numel(FN)
            trialName = FN{trial};
            kdbc1d_Temp = kdbc1d.(trialName);
            baselinedata_trial= kdbc1d_Temp(:,3);
            baselinedata_trial_zero = any(baselinedata_trial==0); %remove all zero (excluded) vectors for averaging
            baselinedata_trial = baselinedata_trial(:, ~baselinedata_trial_zero);
            baseline_F0_ch1 = [baseline_F0_ch1, baselinedata_trial];
        end
        
        baseline_F0_ch1_Avg = mean(baseline_F0_ch1, 2);
        keepdata.F0_Ch1_Baseline_AllReflexTrials = baseline_F0_ch1;
        keepdata.F0_Ch1_Baseline_Avg = baseline_F0_ch1_Avg;
%          F0_Ch1_Baseline_AllReflexTrials = baseline_F0_ch1;
%          F0_Ch1_Baseline_Avg = baseline_F0_ch1_Avg;

        clear kdbc1d
        clear FN
        clear kdbc1d_Temp
        clear trial
        clear trialName

        %% baseline NORMALIZED 100 ms ch1 - INCLUDED TRIALS ONLY every 10 ms

        kdfc = keepdata.F0_Ch1_Baseline;
        FN = fieldnames(kdfc);
        baseline_F0_Norm_ch1 = [];

        for trial = 1:numel(FN)
            trialName = FN{trial};
            kdfc_Temp = kdfc.(trialName);
            
            baselinedata_trialnorm= kdfc_Temp(:,4);
            baselinedata_trialnorm_zero = any(baselinedata_trialnorm==0); %remove all zero (excluded) vectors for averaging
            baselinedata_trialnorm = baselinedata_trialnorm(:, ~baselinedata_trialnorm_zero);
        
            baseline_F0_Norm_ch1 = [baseline_F0_Norm_ch1, baselinedata_trialnorm];

        end

        baseline_F0_Norm_ch1_Avg = mean(baseline_F0_Norm_ch1, 2);
         keepdata.F0_Norm_Ch1_Baseline_AllReflexTrials = baseline_F0_Norm_ch1;
         keepdata.F0_Norm_Ch1_Baseline_Avg = baseline_F0_Norm_ch1_Avg;
         
         keepdata.F0_Norm_Ch1_Baseline_10ms_AllReflexTrials = baseline_F0_Norm_ch1(1:440:end,:);
         keepdata.F0_Norm_Ch1_Baseline_10ms_Avg = baseline_F0_Norm_ch1_Avg(1:440:end);
         
         
         keepdata.F0_Norm_Ch1_Baseline_10ms_AllReflexTrials = keepdata.F0_Norm_Ch1_Baseline_10ms_AllReflexTrials(1:10,:);
         keepdata.F0_Norm_Ch1_Baseline_10ms_Avg = keepdata.F0_Norm_Ch1_Baseline_10ms_Avg(1:10);

%          F0_Norm_Ch1_Baseline_AllReflexTrials = baseline_F0_Norm_ch1;
%          F0_Norm_Ch1_Baseline_Avg = baseline_F0_Norm_ch1_Avg;

        clear kdfc
        clear FN
        clear kdfc_Temp
        clear trial
        clear trialName

        %% perturbation 850 ms ch1 ALL trials 
        
        kdpert = keepdata.F0_Ch1_Pert;
        FN = fieldnames(kdpert);
        Pert_F0_ch1 = [];

        for trial = 1:numel(FN)
        trialName = FN{trial};
        kdpert_Temp = kdpert.(trialName);
        pertdata_trial= kdpert_Temp(:,3); 
         
        % edit1 (start): Hasini Weerathunge 08/02/2018
        if (length(pertdata_trial) < Fs*0.850) 
        mInfvec = -Inf(Fs*0.850+1-length(pertdata_trial), 1); % fill remainder of pert trial after pert array (850ms) with -Inf if array is too short
        pertdata_trial = [pertdata_trial; mInfvec]; %have to add some -Inf at the end to make the files the same length
        end
        % edit1 (end): Hasini Weerathunge 08/02/2018
         
        Pert_F0_ch1 = [Pert_F0_ch1 pertdata_trial];

        end

        Pert_F0_ch1_Avg = mean(Pert_F0_ch1, 2);
%         keepdata.F0_Ch1_Pert_AllReflexTrials = Pert_F0_ch1;
%         keepdata.F0_Ch1_Pert_Avg = Pert_F0_ch1_Avg;
        
        F0_Ch1_Pert_AllReflexTrials = Pert_F0_ch1;
        F0_Ch1_Pert_Avg = Pert_F0_ch1_Avg;

        clear kdpert
        clear FN
        clear kdpert_Temp
        clear trial
        clear trialName

        
        
        
        
        %% perturbation 850 ms ch1 every 10 ms ALL trials
        
        kdpert2 = keepdata.F0_Ch1_Pert;
        FN = fieldnames(kdpert2);
        Pert_F0_ch1_10ms = [];

        for trial = 1:numel(FN)
        trialName = FN{trial};
        kdpert2_Temp = kdpert2.(trialName);
        pertdata_trial_10ms= kdpert2_Temp(:,3); 
        
        % edit2 (start): Hasini Weerathunge 08/02/2018
        if (length(pertdata_trial_10ms) < Fs*0.850) 
        mInfvec = -Inf(Fs*0.850+1-length(pertdata_trial_10ms), 1); % fill remainder of pert trial after pert array (850ms) with -Inf if array is too short
        pertdata_trial_10ms = [pertdata_trial_10ms; mInfvec]; %have to add some -Inf at the end to make the files the same length
        end
        % edit2 (end): Hasini Weerathunge 08/02/2018
        
        
        Pert_F0_ch1_10ms = [Pert_F0_ch1_10ms  pertdata_trial_10ms(1:440:end)]; %every 10 ms matching audapter analysis script

        end

        Pert_F0_ch1_10ms_Avg = mean(Pert_F0_ch1_10ms, 2);
%         keepdata.F0_Ch1_Pert_10ms_AllReflexTrials = Pert_F0_ch1_10ms;
%         keepdata.F0_F0_Ch1_Pert_10ms_Avg = Pert_F0_ch1_10ms_Avg;

         F0_Ch1_Pert_10ms_AllReflexTrials = Pert_F0_ch1_10ms;
         F0_Ch1_Pert_10ms_Avg = Pert_F0_ch1_10ms_Avg;

        clear kdpert2
        clear FN
        clear kdpert2_Temp
        clear trial
        clear trialName

        %% perturbation NORMALIZED 850 ms ch1 ONLY INCLUDED TRIALS
         % and vectors with samples every 10 ms
        
        kdpnorm = keepdata.F0_Ch1_Pert;
        FN = fieldnames(kdpnorm);
        Pert_F0_Norm_ch1 = [];

        for trial = 1:numel(FN)
        trialName = FN{trial};
        kdpnorm_Temp = kdpnorm.(trialName);
        
        pertdata_trialnorm= kdpnorm_Temp(:,4);
        
        % edit3 (start): Hasini Weerathunge 08/02/2018
        if (length(pertdata_trialnorm) < Fs*0.850) 
        mInfvec = -Inf(Fs*0.850+1-length(pertdata_trialnorm), 1); % fill remainder of pert trial after pert array (850ms) with -Inf if array is too short
        pertdata_trialnorm = [pertdata_trialnorm; mInfvec]; %have to add some -Inf at the end to make the files the same length
        end
        % edit3 (end): Hasini Weerathunge 08/02/2018
        
               
        pertdata_trialnorm_zero = any(pertdata_trialnorm==0); %remove all zero (excluded) vectors for averaging
        pertdata_trialnorm = pertdata_trialnorm(:, ~pertdata_trialnorm_zero);
    
        % edit4 (end): Hasini Weerathunge 08/02/2018
        %pertdata_trialnorm_Inf = any(pertdata_trialnorm==-Inf); %remove all -Inf (excluded) vectors for averaging
        %pertdata_trialnorm = pertdata_trialnorm(:, ~pertdata_trialnorm_Inf);
        % edit4 (end): Hasini Weerathunge 08/02/2018
        
        Pert_F0_Norm_ch1 = [Pert_F0_Norm_ch1 pertdata_trialnorm];

        end

        Pert_F0_Norm_ch1_Avg = mean(Pert_F0_Norm_ch1, 2);
%         keepdata.F0_Norm_Ch1_Pert_AllReflexTrials = Pert_F0_Norm_ch1;
%         keepdata.F0_Norm_Ch1_Pert_Avg = Pert_F0_Norm_ch1_Avg;
%         
        keepdata.F0_Norm_Ch1_Pert_AllReflexTrials = Pert_F0_Norm_ch1;
        keepdata.F0_Norm_Ch1_Pert_Avg = Pert_F0_Norm_ch1_Avg;
%        
        %FINAL every 10 ms sample versions of same vectors of included trials
        
        keepdata.F0_Norm_Ch1_Pert_10ms_AllReflexTrials = keepdata.F0_Norm_Ch1_Pert_AllReflexTrials(1:440:end,:);
        keepdata.F0_Norm_Ch1_Pert_10ms_Avg = keepdata.F0_Norm_Ch1_Pert_Avg(1:440:end);
        clear kdpnorm
        clear FN
        clear kdpnorm_Temp
        clear trial
        clear trialName
        
     
    
    
    
%% Run this to exclude trials where variability in baseline is > 15 cents

% Remove trials with variability limit chosen and save new vectors.
% First save original vectors, calculate new averages for the reduced trial set
% Make final vectors with the new trial set and means
% Create a final vector with the averaged baseline and perturbation data
% combined together
% Make a .txt file with all the excluded trials and reason and save in the
% folder, then resave the .mat keepdata file to include all new vectors
% 
%        
% keepdata.Original_All_F0_Norm_Ch1_Baseline_10ms_AllReflexTrials=keepdata.F0_Norm_Ch1_Baseline_10ms_AllReflexTrials;
% keepdata.Original_All_F0_Norm_Ch1_Pert_10ms_AllReflexTrials = keepdata.F0_Norm_Ch1_Pert_10ms_AllReflexTrials; 
% 
% variability check and removal
% 
%
%
%
% keepdata.F0_Norm_Ch1_Baseline_10ms_Avg=mean(keepdata.F0_Norm_Ch1_Baseline_10ms_AllReflexTrials,2);
% keepdata.F0_Norm_Ch1_Pert_10ms_Avg=mean(keepdata.F0_Norm_Ch1_Pert_10ms_AllReflexTrials,2);
% 
% 
% 
% 
% keepdata.FINAL_Baseline_AllReflexTrials=keepdata.F0_Norm_Ch1_Baseline_10ms_AllReflexTrials;
% keepdata.FINAL_Baseline_Avg=mean(keepdata.FINAL_Baseline_AllReflexTrials, 2);
% 
% 
% keepdata.FINAL_Pert_AllReflexTrials = keepdata.F0_Norm_Ch1_Pert_10ms_AllReflexTrials; 
% keepdata.FINAL_Pert_Avg = mean(keepdata.FINAL_Pert_AllReflexTrials, 2);
% 
% keepdata.FINAL_Baseline_and_Pert = [keepdata.FINAL_Baseline_Avg; keepdata.FINAL_Pert_Avg];
% 
% 
% 

               


    
    

    end
    
cd(PartFolder);


if exist( savefilename )
save(savefilename, 'keepdata'); 
else
keepdata.AllCh1Ch2_datalables = [{'Time(s)'}, {'intensity(dB)'}, {'F0(Hz)'}, {'F1(Hz)'}, {'F2(Hz)'}];
keepdata.F0_Ch1_Baseline_datalabels = [{'Baseline Time(s)'}, {'Baseline intensity(dB)'}, {'Baseline F0(Hz)'}, {'Baseline F0 norm(cents)'}];
keepdata.F0_Ch1_Pert_datalabels = [{'Time(s)'}, {'intensity(dB)'}, {'F0(Hz)'}, {'F0 norm(cents)'}];


%% save full keepdata

save(savefilename, 'keepdata');   

%% Run offline for regeneration of vectors:

%%create baseline normalized mean every 10 ms

% keepdata.F0_Norm_Ch1_Baseline_10ms_Avg=mean(keepdata.F0_Norm_Ch1_Baseline_10ms_AllReflexTrials,2);

% create perturbation normalized mean every 10 ms

% keepdata.F0_Norm_Ch1_Pert_10ms_Avg=mean(keepdata.F0_Norm_Ch1_Pert_10ms_AllReflexTrials,2);

% keepdata.FINAL_Full_Norm_Avg_Vector =[keepdata.F0_Norm_Ch1_Baseline_10ms_Avg; keepdata.F0_Norm_Ch1_Pert_10ms_Avg];
% plot(keepdata.FINAL_Full_Norm_Avg_Vector)

%%

end

end
%% functions for showing/hiding formants, pitch, delay, 
    function showformants_ch1(PushButtonFormant_ch1, EventData)
        axes(sub3)
        
        yyaxis left
        
        F11=  plot(sub3, timeVec,F1_use_ch1, 'b', 'LineWidth', 1.5);
        F21= plot(sub3, timeVec,F2_use_ch1, 'g', 'LineWidth', 1.5);
        % plot(sub3, timeVec,F3_use_ch1, 'm', 'LineWidth', 1.5)
        hold on
        
    end

    function hideformants_ch1(PushButtonFormantHide_ch1, EventData)
        axes(sub3)
        delete(F11)
        delete(F21)
        
    end

    function showformants_ch2(PushButtonFormant_ch2, EventData)
        axes(sub4)
        yyaxis left
        F12 =  plot(sub4, timeVec,F1_use_ch2, 'b', 'LineWidth', 1.5);
        F22 = plot(sub4, timeVec,F2_use_ch2, 'g', 'LineWidth', 1.5);
        hold on
        
    end

    function hideformants_ch2(PushButtonFormantHide_ch2, EventData)
        axes(sub4)
        delete(F12)
        delete(F22)
        
    end
    function showpitch_ch2(PushButtonPitch_ch2, EventData)
        axes(sub4)
        yyaxis right
        
        
        F02 = plot(sub4, timeVec,F0_use_ch2, 'r', 'LineWidth', 1.5);
        hold on
        
    end
    function hidepitch_ch2(PushButtonPitchHide_ch2, EventData)
        %plot the pitch on a secondary y axis
        axes(sub4)
        delete(F02);
        
    end
    function showpitch_ch1(PushButtonPitch_ch1, EventData)
        %plot the pitch on a secondary y axis
        axes(sub3)
        yyaxis right
        
        
        F01 = plot(sub3, timeVec,F0_use_ch1, 'r', 'LineWidth', 1.5);
        
        hold on
    end
    function hidepitch_ch1(PushButtonPitchHide_ch1, EventData)
        %plot the pitch on a secondary y axis
        axes(sub3)
        delete(F01)
    end

    function callbackspect(PushButtonspect, EventData)
        
        if reflexShift == 1 %if a shifted reflex trial
            [x_spect,y_spect] = ginput(2); %only doing timing in output

            timex_spect = x_spect;
            timey_spect = y_spect;

            headphonestart = timex_spect(1); %first click is start of steady state perturbation
            %headphoneend = timex_spect(2);

            headphonestart_fullreflex = timex_spect(2);%second click is when perturbation is starting, first point not at baseline
            %fullreflexend = timex_spect(4);

            %micstart = timex_spect(1) - delaycalc; % start of sustained reflex, when it reaches steadystate
            %micstart_fullreflex = timex_spect(2) - delaycalc; %start of full reflex, first point when it isn't baseline
            micstart = timex_spect(1); % start of sustained reflex, when it reaches steadystate # editted and removed delay alignment (edit hasini 05/10/2018)
            micstart_fullreflex = timex_spect(2) ; %start of full reflex, first point when it isn't baseline # editted and removed delay alignment (edit hasini 05/10/2018)
            
            lengthVoiceOnset_Reflex = micstart_fullreflex - timeVoiceOnset;  %time from voice onset detection to start of reflex

            headphonebaselinestart = headphonestart_fullreflex - .1; %calculating 100 ms prior to onset of perturbation as baseline
            micbaselinestart = micstart_fullreflex - .1;

            headphoneperturbend = headphonestart + .1; %calculated 100 ms from the start of sustained perturbation
            micperturbend = micstart + 0.85; %saving first 850 ms from start of perturbation

            [c headphonebaselineStartindex] = min(abs(timeVec - headphonebaselinestart)); %100 ms before start of perturbation
            [c headphonebaselineEndindex] = min(abs(timeVec - headphonestart_fullreflex)); %start of perturbation

            [c micbaselineStartindex] = min(abs(timeVec - micbaselinestart)); %same as headphones above, except shifted for microphone trials #(comment hasini 05/10/2018)<- not shifted anymore
            [c, micbaselineEndindex] = min(abs(timeVec - micstart_fullreflex));


            [c headphoneStartindex] = min(abs(timeVec - headphonestart)); %start of steady state perturbation
            [c headphoneEndindex] = min(abs(timeVec - headphoneperturbend)); %850 ms (100ms??) after start of steady state perturbation

            [c micStartindex] = min(abs(timeVec - micstart)); %same as headphones above, except shifted for microphone trials #(comment hasini 05/10/2018)<- not shifted anymore
            [c micEndindex] = min(abs(timeVec - micperturbend));




            %ch1, 100 ms before start of perturbation
            baselineTime_ch1 = timeVec(micbaselineStartindex:micbaselineEndindex);
            baselineF0_ch1 = F0_use_ch1(micbaselineStartindex:micbaselineEndindex);
            baselineF1_ch1 = F1_use_ch1(micbaselineStartindex:micbaselineEndindex);
            baselineF2_ch1 = F2_use_ch1(micbaselineStartindex:micbaselineEndindex);
            baselineInt_ch1 = Int_use_ch1(micbaselineStartindex:micbaselineEndindex);
            baselineF0_ch1_mean = mean(baselineF0_ch1);

            %vectors for analysis
            norm_baselineF0_ch1= 1200*log2(baselineF0_ch1./mean(baselineF0_ch1));

            %ch2
            baselineTime_ch2 = timeVec(micbaselineStartindex:micbaselineEndindex);
            baselineF0_ch2 = F0_use_ch2(micbaselineStartindex:micbaselineEndindex);
            baselineF1_ch2 = F1_use_ch2(micbaselineStartindex:micbaselineEndindex);
            baselineF2_ch2 = F2_use_ch2(micbaselineStartindex:micbaselineEndindex);
            baselineInt_ch2 = Int_use_ch2(micbaselineStartindex:micbaselineEndindex);
            baselineF0_ch2_mean = mean(baselineF0_ch2);


                       % save information from ch 1 and ch 2
            %channel1
            selectTime_ch1 = timeVec(micStartindex:micEndindex);
            selectF0_ch1 = F0_use_ch1(micStartindex:micEndindex);
            selectF1_ch1 = F1_use_ch1(micStartindex:micEndindex);
            selectF2_ch1 = F2_use_ch1(micStartindex:micEndindex);
            selectInt_ch1 = Int_use_ch1(micStartindex:micEndindex);
            
            %vectors for analysis
            norm_selectF0_ch1= 1200*log2(selectF0_ch1./mean(baselineF0_ch1));
            
      
           %%
           %%

            %channel2
            selectTime_ch2 = timeVec(headphoneStartindex:headphoneEndindex);
            selectF0_ch2 = F0_use_ch2(headphoneStartindex:headphoneEndindex);
            selectF1_ch2 = F1_use_ch2(headphoneStartindex:headphoneEndindex);
            selectF2_ch2 = F2_use_ch2(headphoneStartindex:headphoneEndindex);
            selectInt_ch2 = Int_use_ch2(headphoneStartindex:headphoneEndindex);
            
               %calculate all means only for 300-500 ms after perturbation steady
            %state start

            headphone_meanstart = headphonestart + .15; %150 ms after start of steady state perturbation
            headphone_meanend = headphone_meanstart + .15 ; % 150 ms after start of mean measurement (300 ms after start of steady state pert)

            mic_meanstart = micstart + .15; %150 ms after start of steady state perturbation
            mic_meanend = mic_meanstart + .15 ; % 150 ms after start of mean measurement (300 ms after start of steady state pert)


            [c headphoneMeanStartindex] = min(abs(timeVec - headphone_meanstart)); 
            [c headphoneMeanEndindex] = min(abs(timeVec - headphone_meanend)); 

            [c micMeanStartindex] = min(abs(timeVec - mic_meanstart)); 
            [c micMeanEndindex] = min(abs(timeVec - mic_meanend));

            %index that will be saved for mean calculations for shifted trials
                    %channel1
            selectTime_ch1Mean = timeVec(micMeanStartindex:micMeanEndindex);
            selectF0_ch1Mean = F0_use_ch1(micMeanStartindex:micMeanEndindex);
            selectF1_ch1Mean = F1_use_ch1(micMeanStartindex:micMeanEndindex);
            selectF2_ch1Mean = F2_use_ch1(micMeanStartindex:micMeanEndindex);
            selectInt_ch1Mean = Int_use_ch1(micMeanStartindex:micMeanEndindex);

                        %channel2
            selectTime_ch2Mean = timeVec(headphoneMeanStartindex:headphoneMeanEndindex);
            selectF0_ch2Mean = F0_use_ch2(headphoneMeanStartindex:headphoneMeanEndindex);
            selectF1_ch2Mean = F1_use_ch2(headphoneMeanStartindex:headphoneMeanEndindex);
            selectF2_ch2Mean = F2_use_ch2(headphoneMeanStartindex:headphoneMeanEndindex);
            selectInt_ch2Mean = Int_use_ch2(headphoneMeanStartindex:headphoneMeanEndindex);


            F0calc_ch1 = mean(selectF0_ch1Mean);
            F1calc_ch1 = mean(selectF1_ch1Mean*1000);
            F2calc_ch1 = mean(selectF2_ch1Mean*1000);
            Intcalc_ch1 = mean(selectInt_ch1Mean);

            F0calc_ch2 = mean(selectF0_ch2Mean);
            F1calc_ch2 = mean(selectF1_ch2Mean*1000);
            F2calc_ch2 = mean(selectF2_ch2Mean*1000);
            Intcalc_ch2 = mean(selectInt_ch2Mean);


             %what to display

            F0textbx_ch1.EdgeColor = [1 0 0];
            F0textbx_ch1.String = {F0calc_ch1};

            F1textbx_ch1.EdgeColor = [1 0 0];
            F1textbx_ch1.String = {F1calc_ch1};

            F2textbx_ch1.EdgeColor = [1 0 0];
            F2textbx_ch1.String = {F2calc_ch1};

            F0textbx_ch2.EdgeColor = [1 0 0];
            F0textbx_ch2.String = {F0calc_ch2};

            F1textbx_ch2.EdgeColor = [1 0 0];
            F1textbx_ch2.String = {F1calc_ch2};

            F2textbx_ch2.EdgeColor = [1 0 0];
            F2textbx_ch2.String = {F2calc_ch2};

            axes(sub3) %select subplot3
            line([mic_meanstart mic_meanstart], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
            line([mic_meanend mic_meanend], ylim, 'color', [1 1 1], 'LineWidth', 1.5)

            axes(sub4) %select subplot4
            line([headphone_meanstart headphone_meanstart], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
            line([headphone_meanend headphone_meanend], ylim, 'color', [1 1 1], 'LineWidth', 1.5)

        else %if a non-shifted reflex trial or an adaptation trial
            [x_spect,y_spect] = ginput(2); %only doing timing in output

            timex_spect = x_spect;
            timey_spect = y_spect;

            headphonestart = timex_spect(1); %marked start of trial
            headphoneend = timex_spect(2); %marked end of trial
            
            micstart = timex_spect(1)-delaycalc; %Edit 3 Hasini uncommented: shift the mic back the delay time, relative to the headphones. 
            %micend = timex_spect(2)-delaycalc;
            %micstart = timex_spect(1); %Edit 3 Hasini commented #(edit hasini 05/10/2018)<- mic signal is not shifted back relative to the headphone signal. This fix is specifically for DFPAR study only as the delays significantly huge! range [20-150ms]
            micend = timex_spect(2); % #(edit hasini 05/10/2018)<- mic signal is not shifted back relative to the headphone signal. This fix is specifically for DFPAR study only as the delays significantly huge! range [20-150ms]
           
            timeVoiceOnset_ch1= micstart; % Edit 3 Hasini : start saving onsets
            timeVoiceOnset_ch2 =headphonestart; % Edit 3 Hasini : start saving onsets
            timeVoiceOffset = headphoneend; % Edit 3 Hasini : start saving offsets
        
            %find indicies that relate to the above values
            [c headphoneStartindex] = min(abs(timeVec - headphonestart)); %find the vector point closests to the calculated start/end time
            [c headphoneEndindex] = min(abs(timeVec - headphoneend));


            [c micStartindex] = min(abs(timeVec - micstart));
            [c micEndindex] = min(abs(timeVec - micend));
            
             %channel1 ( select headphone voice onset to headphone voice
            %offsetfor mean f0 calculation ) %Edit 3 Hasini
            selectTime_ch1_mean = timeVec(headphoneStartindex:micEndindex);
            selectF0_ch1_mean = F0_use_ch1(headphoneStartindex:micEndindex);
            selectF1_ch1_mean = F1_use_ch1(headphoneStartindex:micEndindex);
            selectF2_ch1_mean = F2_use_ch1(headphoneStartindex:micEndindex);
            selectInt_ch1_mean = Int_use_ch1(headphoneStartindex:micEndindex);
                   
            
            % save information from ch 1 and ch 2
            %channel1
            selectTime_ch1 = timeVec(micStartindex:micEndindex);
            selectF0_ch1 = F0_use_ch1(micStartindex:micEndindex);
            selectF1_ch1 = F1_use_ch1(micStartindex:micEndindex);
            selectF2_ch1 = F2_use_ch1(micStartindex:micEndindex);
            selectInt_ch1 = Int_use_ch1(micStartindex:micEndindex);

            %channel2
            selectTime_ch2 = timeVec(headphoneStartindex:headphoneEndindex);
            selectF0_ch2 = F0_use_ch2(headphoneStartindex:headphoneEndindex);
            selectF1_ch2 = F1_use_ch2(headphoneStartindex:headphoneEndindex);
            selectF2_ch2 = F2_use_ch2(headphoneStartindex:headphoneEndindex);
            selectInt_ch2 = Int_use_ch2(headphoneStartindex:headphoneEndindex);
            
            
            % What tp display % Edit 3 Hasini
            F0calc_ch1 = mean(selectF0_ch1_mean);      %F0calc_ch1 = mean(selectF0_ch1);
            F1calc_ch1 = mean(selectF1_ch1_mean*1000); %F1calc_ch1 = mean(selectF1_ch1*1000);
            F2calc_ch1 = mean(selectF2_ch1_mean*1000); %F2calc_ch1 = mean(selectF2_ch1*1000);
            Intcalc_ch1 = mean(selectInt_ch1_mean);    %Intcalc_ch1 = mean(selectInt_ch1);
           
            F0calc_ch2 = mean(selectF0_ch2);
            F1calc_ch2 = mean(selectF1_ch2*1000);
            F2calc_ch2 = mean(selectF2_ch2*1000);
            Intcalc_ch2 = mean(selectInt_ch2);

             %what to display

            F0textbx_ch1.EdgeColor = [1 0 0];
            F0textbx_ch1.String = {F0calc_ch1};

            F1textbx_ch1.EdgeColor = [1 0 0];
            F1textbx_ch1.String = {F1calc_ch1};

            F2textbx_ch1.EdgeColor = [1 0 0];
            F2textbx_ch1.String = {F2calc_ch1};

            F0textbx_ch2.EdgeColor = [1 0 0];
            F0textbx_ch2.String = {F0calc_ch2};

            F1textbx_ch2.EdgeColor = [1 0 0];
            F1textbx_ch2.String = {F1calc_ch2};

            F2textbx_ch2.EdgeColor = [1 0 0];
            F2textbx_ch2.String = {F2calc_ch2};



            axes(sub3) %select subplot3
            line([micstart micstart], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
            line([headphonestart headphonestart], ylim, 'color', [1 1 1], 'LineWidth', 1.5) % Edit 3 Hasini
            line([micend micend], ylim, 'color', [1 1 1], 'LineWidth', 1.5)

            axes(sub4) %select subplot4
            line([headphonestart headphonestart], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
            line([headphoneend headphoneend], ylim, 'color', [1 1 1], 'LineWidth', 1.5)


            %these are not used for these trials, so set them to 0 
            %  lengthreflex = 0; %if not a reflex script, don't change this
            lengthVoiceOnset_Reflex = 0;
            %lengthfullreflex =0;
            baselineF0_ch1_mean = 0;
            baselineF0_ch2_mean = 0;
        
        end
        
    end
    function callbackdelay (PushButtondelay, EventData)
        
        [delayx,delayy] = ginput(2);
        
        timex_delay = delayx;
        timey_delay = delayy;
        
        %display delay
        ch1_delaystart = timex_delay(1);
        ch2_delaystart = timex_delay(2);
        delaycalc = ch2_delaystart - ch1_delaystart;
        delaycalcsec = delaycalc*1000;
        delaytextbx.EdgeColor = [1 0 0];
        delaytextbx.String = {delaycalcsec};
        
        set(PushButtonspect, 'Visible', 'on');
        
        if reflexShift == 1
            axes(sub1)
            line([timeVoiceOnset timeVoiceOnset], ylim, 'color', [1 0 0], 'LineWidth', 1.5) %plot the voicing onset
        end
        
    end
    function callbacknext (PushButtonnext, EventData)
        
        %%for saving the information
        keepvalue = {get(keeepcheckbox, 'Value')};
        redovalue = {get(redocheckbox, 'Value')};
        whyredo = {get(whyredobox, 'string')};
        close all
    end

    function playpitch_ch1 (PushButtonplay_ch1, EventData)
        cd(micFolder)
        sound(ch1_wav,Fs)
    end
    function playpitch_ch2 (PushButtonplay_ch1, EventData)
        cd(headphonesFolder)
        sound(ch2_wav,Fs)
    end




% keepdata.ch1 and keepdata.ch2 information columnes 1-5
%1: Time selected to analyze signal
%2: Intensity values over selected signal
%3: F0 values over selected signal
%4: F1 values over selected signal
%5: F2 values over selected signal


end




function [cs,index] = sort_nat(c,mode)
%sort_nat: Natural order sort of cell array of strings.
% usage:  [S,INDEX] = sort_nat(C)
%
% where,
%    C is a cell array (vector) of strings to be sorted.
%    S is C, sorted in natural order.
%    INDEX is the sort order such that S = C(INDEX);
%
% Natural order sorting sorts strings containing digits in a way such that
% the numerical value of the digits is taken into account.  It is
% especially useful for sorting file names containing index numbers with
% different numbers of digits.  Often, people will use leading zeros to get
% the right sort order, but with this function you don't have to do that.
% For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
% will give you
%
%       {'file1.txt'  'file10.txt'  'file2.txt'}
%
% whereas, sort_nat will give you
%
%       {'file1.txt'  'file2.txt'  'file10.txt'}
%
% See also: sort

% Version: 1.4, 22 January 2011
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Set default value for mode if necessary.
if nargin < 2
    mode = 'ascend';
end

% Make sure mode is either 'ascend' or 'descend'.
modes = strcmpi(mode,{'ascend','descend'});
is_descend = modes(2);
if ~any(modes)
    error('sort_nat:sortDirection',...
        'sorting direction must be ''ascend'' or ''descend''.')
end

% Replace runs of digits with '0'.
c2 = regexprep(c,'\d+','0');

% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';

% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c,'\d+','match','start','end');

% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
    num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
    num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end

% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);

% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);

% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;

% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;

% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);

% Sort rows of composite matrix and use index to sort c in ascending or
% descending order, depending on mode.
[unused,index] = sortrows(comp);
if is_descend
    index = index(end:-1:1);
end
index = reshape(index,size(c));
cs = c(index);
end