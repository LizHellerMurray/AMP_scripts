function MainGui_Adaptation_SustainedReflex()
%%things to change for individual study
%FolderName (line 18): where data is stored
%DestFolder (line 33):where you want analysis to be stored
% Any naming convention you want changed for the analysis folders, starting at line 50
% Numbers starting at line 951, if you want the saved average to be different. currenlty average is 150 ms -300 ms after start of steady state perturbation. note: saves all data from 200 ms before to 1000 ms after perturbation, so can alcaulte average offline if wanted as well

%%OF NOTE: for reflex trial, the first click should be at start of where the perturbation
%has reached steady state and second should be the start of the
%perturbation (i.e., the first point it deviates from baseline) 
%05/11/2018 EHM
%Delay calculated but not being used for time alignment, delay between mic
%and headphones is left as is

%12/13/2018 EHM edits
% editted so voice onset is selected by the user during the selection of
% delay
%editted so cannot select analysis portion until delay/voice onset are
%check
%during adaptation --> microphone is saved from the voice onset to end of
%trial - shown in figure. times selected in headphones are selected for
%analysis
%during reflex --> baseline, perturbation, and end shown in figure

close all
clear all

participant = input('Subject ID? ', 's'); %input subject ID numberd

%% change for individual studies
FolderName = uigetdir('R:\SteppLab2\SPEECHDATA\Voice\VH_Pitch_Adaptation_and_Reflex\Adult_Pediatric_Control','Select Data Folder'); %choose data folder; make sure there are "mic" and "perturbed_mic" subfolders in this main data folder
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

%% change for individual studies
DestFolder = uigetdir('R:\SteppLab3\Projects\Voice\SchoolAge_Studies\Perturbation_Study\Analysis\AnalysisFolder','Select Destination Folder'); %choose folder to save data and figures in
cd(DestFolder)

if ~exist(participant, 'dir')
mkdir(participant) %make folder for participant
else 
end

PartFolder = [DestFolder '\' participant]; %set a path for the participant folder

cd(PartFolder)

%% defines which experiment analyzing
choiceType = questdlg('Is this Adaptation or Reflex?', ...
    'Script Type', ...
    'Adaptation', 'Reflex','Reflex');

switch choiceType
    case 'Adaptation'
        scripttype = 1; %adaptation script
        
        experimentType = questdlg('Is this Shift-up or Shift-down?', ...
            'Experiment Type', ...
            'Shift-Up','Shift-Down','Control','Control');
        switch experimentType
            case  'Shift-Up'
                savetype = '_Shift-Up';
                
                if ~exist('Adapt_Shift_Up', 'dir')
                    mkdir Adapt_Shift_Up  %shift up directory
                else
                end
                PartFolder = [PartFolder '\' 'Adapt_Shift_Up']; %rename participant folder to appropriate script type
                
            case  'Shift-Down'
                savetype = '_Shift-Down';
                
                if ~exist('Adapt_Shift_Down', 'dir')
                    mkdir Adapt_Shift_Down  %shift down directory
                else
                end
                
                PartFolder = [PartFolder '\' 'Adapt_Shift_Down']; %rename participant folder to appropriate script type
            case  'Control'
                savetype = '_Control';
                
                if ~exist('Adapt_Control', 'dir')
                    mkdir Adapt_Control  %control directory
                else
                end
                
                PartFolder = [PartFolder '\' 'Adapt_Control']; %rename participant folder to appropriate script type
                
        end
        
    case 'Reflex'
        scripttype = 2; %reflex script
        
        %%UPDATE THIS TO WHAT YOU WANT YOUR REFLEX FOLDERS SAVED AS
        experimentType = questdlg('Which reflex?', ...
            'Which Reflex?', ...
            'Reflex#1', 'Reflex#2','Reflex#1');
        switch experimentType
            case  'Reflex#1'
                savetype = 'Reflex#1';
                
                if ~exist('Reflex#1', 'dir')
                    mkdir Reflex#1
                else
                end
                
                PartFolder = [PartFolder '\' 'Reflex#1']; %rename participant folder to appropriate script type
                %  firstreflex = 1;
            case  'Reflex#2'
                savetype = 'Reflex#2';
                
                if ~exist('Reflex#2', 'dir')
                    mkdir Reflex#2
                else
                end
                
                PartFolder = [PartFolder '\' 'Reflex#2']; %rename participant folder to appropriate script type
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
d = dir([FolderName, '\trial*.wav']); %set the directory
fnames = sort_nat({d.name}); %sort directory alphabetically in order


%% splits and saves ch1 (mic signal) and ch2 (headphones) output for all files
ch1filenames_save = [];
ch2filenames_save = [];

splitwavfiles = input('Do you need to split wav files(1=yes, 0=no) ', 's');  %type yes if you need to run through splitting wav files if this is your first time analyzing this subject
%splitwavfiles = 0; 
str2num(splitwavfiles);
if splitwavfiles ==1
    for j = 1:length(fnames) %creates loop, for 1 until the last file
        filename = strcat(FolderName,'\',fnames{j}); %creates full file name using path and wav name
        [Y, Fs] = audioread(filename); %loads wav file
        ch1filename = strcat(FolderName,'\mic','\',fnames{j}); % makes new file path to "mic" folder
        ch2filename = strcat(FolderName,'\headphones','\',fnames{j}); % makes new file path to "mic" folder
        
        
        micdata  = Y(:,1);
        headphonedata = Y(:,2);
        
   
        audiowrite(ch1filename,micdata,Fs) %creates audio file
        audiowrite(ch2filename,headphonedata,Fs) %creates audio file
        
        ch1filenames_save = [ch1filenames_save; {ch1filename}];
        ch2filenames_save = [ch2filenames_save; {ch2filename}];
        
    end
else
end

%go back to destination folder and go one up to all analysis scripts. If the
%analysis folder is changed then THIS needs to be changed to navigate to
%the correct folder for the scripts
cd(DestFolder); %go to destination folder
cd ../ %go up one to folder where praat interface is

praatfiles = input('Do you need to read in/analyze the praat files to create a txt file? (1=yes, 0=no) ', 's');  %type yes if you need to run through the praat script if this is your first time analyzing this subject
%praatfiles = 0;
str2num(praatfiles);

% voicingOnset = 0.005; %%input('voicing onset? (default is .005):');
% if isempty(voicingOnset)
%     voicingOnset = 0.005;
% end

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
    gp_fn = fullfile( pi_dir , 'get_analysis.praat' ) ;
    if ~exist( gp_fn , 'file' )
        error( 'file ''get_analysis.praat'' not found' )
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
        fullfile( pwd , gp_fn ) , ... %location of praat script wrote (get_analysis.praat)
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
            error( [ 'ERROR: get_analysis.praat failed to run. Make sure '...
                'there are no spaces in the selected filepaths.'])
        end
    end
    
    
    %% headphones
    newheadphones = [headphonesFolder, '\']; %add a slash to the headphones folder
    ext = '.wav'; %extension of files
    
    %Build DOS calls to control praat
    call3 = sprintf( '%s praat "execute %s %s %s %s %f %f %f %f %f"' , ...
        fullfile( pwd , sp_fn ) , ... %location of sendpraat.exe
        fullfile( pwd , gp_fn ) , ... %location of praat script wrote (get_analysis.praat)
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
            error( [ 'ERROR: get_analysis.praat failed to run. Make sure '...
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
keepdata.average = [{'trial'}, {'1 if keep was checked'}, {'1 if redo was checked'}, {'reason why redo was checked'}, {'delay'}, {'length from voice onset detection to start of reflex'}, {'voice onset time (shifted reflex only)'},  {'avg Int(ch1)'}, {'avg Int(ch2)'},{'avg F0 basleine'},{'avg F0(ch1)'}, {'avg F0(ch2)'},{'avg F1(ch1)'}, {'avg F1(ch2)'},{'avg F2(ch1)'}, {'avg F2(ch2)'}];

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
    keeepcheckbox = uicontrol('style','checkbox','String', 'Keep?', 'units', 'normalized', 'Position', [.85 .12 .08 .02]);
    redocheckbox = uicontrol('style','checkbox','String', 'Do Not Use?', 'units', 'normalized', 'Position', [.85 .09 .08 .02]);
    
    %set up field to fill in why you clicked redo
    whyredobox = uicontrol('style','edit','String', '0', 'units', 'normalized', 'Position', [.85 .04 .1 .05]);
    
    %set of variables to be passed between nested functions
    selectTime_ch1 = [];
    selectF0_ch1 = [];
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
    PushButtondelay = uicontrol(gcf, 'Style', 'push', 'String', 'Press when ready to select Voice Onset in both channels (zoom in first- start with ch1)', 'Callback',@callbackdelay, 'units', 'normalized', 'Position' , [.02 .04 .2 .05]  );
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
    starttime_spec = timeVec(intStartIndex_Voicing); %based on intensity calculations above, this is the start time we want in teh window
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
        
        %all trials are shifted in pediatric script
        Reflexcond = questdlg('Is a shifted Reflex trial?', ...
            'Reflex Type', ...
            'Yes','No','Yes');
        
       switch Reflexcond
           case 'Yes'
                reflexShift =1;
                axes(sub1)
                
%                 I = find(ch1_wav>voicingOnset,1,'first');
%                 timeVoiceOnset = I/44100;
%                 
%                 %if timeVoiceOnset is empty, timevoiceOnset =0
%                    if isempty(timeVoiceOnset)
%                        timeVoiceOnset = 0;
%                    end
                
                
                axes(sub3)
                newxlim3_reflex = [.7 1.7]; %zoom in to .5 to 1.5 so can better see reflex
                set(sub3, 'Xlim', newxlim3_reflex)
                yyaxis right
                F01 = plot(sub3, timeVec,F0_use_ch1, 'r', 'LineWidth', 1.5);
                hold on
                
                
                axes(sub4)
                
                newxlim4_reflex = [.7 1.7]; %zoom in to .5 to 1.5 so can better see reflex
                set(sub4, 'Xlim', newxlim4_reflex)
                yyaxis right
                F02 = plot(sub4, timeVec,F0_use_ch2, 'r', 'LineWidth', 1.5);
                hold on
           case 'No'
               reflexShift =2;
             %  timeVoiceOnset = 0; %if it is not a shifted reflex trial
              axes(sub3)

                yyaxis right
                F01 = plot(sub3, timeVec,F0_use_ch1, 'r', 'LineWidth', 1.5);
                hold on
                
                
                axes(sub4)

                yyaxis right
                F02 = plot(sub4, timeVec,F0_use_ch2, 'r', 'LineWidth', 1.5);
                hold on
                
       end
              
    else %if it is an adaptation trial
        reflexShift =0;
     %   timeVoiceOnset = 0;
        
        axes(sub3)
        yyaxis right
        F01 = plot(sub3, timeVec,F0_use_ch1, 'r', 'LineWidth', 1.5);
        hold on
        
        
        axes(sub4)
        yyaxis right
        F02 = plot(sub4, timeVec,F0_use_ch2, 'r', 'LineWidth', 1.5);
        hold on
         
        
        
    end
    delaycalc=0;
    waitfor(mFigure)
    %% saving
   % keepdata.keepcode = [keepdata.keepcode; {trialnum}, keepvalue, redovalue,whyredo];
 %   keepdata.average = [ keepdata.average; {trialnum}, keepvalue, redovalue,whyredo, lengthreflex,lengthfullreflex, delaycalc, lengthVoiceOnset_Reflex,  timeVoiceOnset, Intcalc_ch1, Intcalc_ch2, F0calc_ch1, F0calc_ch2, F1calc_ch1, F1calc_ch2, F2calc_ch1, F2calc_ch2];
      keepdata.average = [ keepdata.average; {trialnum}, keepvalue, redovalue,whyredo, delaycalc, lengthVoiceOnset_Reflex,  timeVoiceOnset, Intcalc_ch1, Intcalc_ch2, baselineF0_ch1_mean, F0calc_ch1, F0calc_ch2, F1calc_ch1, F1calc_ch2, F2calc_ch1, F2calc_ch2];
    keepdata.ch1data.(trialnum) = [selectTime_ch1, selectInt_ch1,selectF0_ch1, selectF1_ch1*1000,selectF2_ch1*1000];
    keepdata.ch2data.(trialnum) = [selectTime_ch2, selectInt_ch2,selectF0_ch2,selectF1_ch2*1000,selectF2_ch2*1000];
    if   reflexShift == 1
        keepdata.baselinech1data.(trialnum) = [baselineTime_ch1, baselineInt_ch1,baselineF0_ch1, baselineF1_ch1*1000,baselineF2_ch1*1000];
        keepdata.baselinech2data.(trialnum) = [baselineTime_ch2, baselineInt_ch2,baselineF0_ch2,baselineF1_ch2*1000,baselineF2_ch2*1000];
        
    end
    
 cd(PartFolder);


if exist( savefilename )
save(savefilename, 'keepdata'); 
else
keepdata.datalables = [{'Time(s)'}, {'intensity(dB)'}, {'F0(Hz)'}, {'F1(Hz)'}, {'F2(Hz)'}];
save(savefilename, 'keepdata');   
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
        
   %edited 12/13/2018 to clean up redundancies and comments added
        if reflexShift == 1 %if a shifted reflex trial
        [x_spect,y_spect] = ginput(2); %only doing timing in output
        
        timex_spect = x_spect; %only using "X"
        timey_spect = y_spect;
        
        %same point in time in microphone and headphones, 
        shiftStart = timex_spect(1); %first click is start of steady state perturbation      
        baselineEnd = timex_spect(2);%second click is when perturbation is starting, first point not at baseline
                
        lengthVoiceOnset_Reflex = shiftStart - timeVoiceOnset;  %time from voice onset detection to start of reflex
        
        baselineStart = baselineEnd - .2; %calculating 200 ms prior to onset of perturbation as baseline
        shiftEnd = shiftStart + 1; %saving first 1000 ms from start of perturbation    

        
        %find relevant indicaties
        [c baselineStartindex] = min(abs(timeVec - baselineStart)); %200 ms before start of perturbation
        [c baselineEndindex] = min(abs(timeVec - baselineEnd)); %end of baseline period, when begins deviating from baseline
        
        [c Startindex] = min(abs(timeVec - shiftStart)); %start of steady state perturbation
        [c Endindex] = min(abs(timeVec - shiftEnd)); %1000 ms after start of steady state perturbation

        
      
        % save baseline
        %ch1, 200 ms before start of perturbation
        baselineTime_ch1 = timeVec(baselineStartindex:baselineEndindex);
        baselineF0_ch1 = F0_use_ch1(baselineStartindex:baselineEndindex);
        baselineF1_ch1 = F1_use_ch1(baselineStartindex:baselineEndindex);
        baselineF2_ch1 = F2_use_ch1(baselineStartindex:baselineEndindex);
        baselineInt_ch1 = Int_use_ch1(baselineStartindex:baselineEndindex);
        baselineF0_ch1_mean = mean(baselineF0_ch1); %save mean baseline f0
            
         %ch2
        baselineTime_ch2 = timeVec(baselineStartindex:baselineEndindex);
        baselineF0_ch2 = F0_use_ch2(baselineStartindex:baselineEndindex);
        baselineF1_ch2 = F1_use_ch2(baselineStartindex:baselineEndindex);
        baselineF2_ch2 = F2_use_ch2(baselineStartindex:baselineEndindex);
        baselineInt_ch2 = Int_use_ch2(baselineStartindex:baselineEndindex);
        baselineF0_ch2_mean = mean(baselineF0_ch2);
        
        
        % save reflex information from ch 1 and ch 2, 0 to 1000 ms
        %channel1
        selectTime_ch1 = timeVec(Startindex:Endindex);
        selectF0_ch1 = F0_use_ch1(Startindex:Endindex);
        selectF1_ch1 = F1_use_ch1(Startindex:Endindex);
        selectF2_ch1 = F2_use_ch1(Startindex:Endindex);
        selectInt_ch1 = Int_use_ch1(Startindex:Endindex);
     
                    %channel2
        selectTime_ch2 = timeVec(Startindex:Endindex);
        selectF0_ch2 = F0_use_ch2(Startindex:Endindex);
        selectF1_ch2 = F1_use_ch2(Startindex:Endindex);
        selectF2_ch2 = F2_use_ch2(Startindex:Endindex);
        selectInt_ch2 = Int_use_ch2(Startindex:Endindex);
        
               
        %calculate means  for 150-300 ms after perturbation steady
         meanstart = shiftStart + .15; %150 ms after start of steady state perturbation
         meanend = meanstart + .15 ; % 150 ms after start of mean measurement (300 ms after start of steady state pert)
        
        
        [c MeanStartindex] = min(abs(timeVec - meanstart)); 
        [c MeanEndindex] = min(abs(timeVec - meanend)); 

        %index that will be saved for mean calculations 
                %channel1
        selectTime_ch1Mean = timeVec(MeanStartindex:MeanEndindex);
        selectF0_ch1Mean = F0_use_ch1(MeanStartindex:MeanEndindex);
        selectF1_ch1Mean = F1_use_ch1(MeanStartindex:MeanEndindex);
        selectF2_ch1Mean = F2_use_ch1(MeanStartindex:MeanEndindex);
        selectInt_ch1Mean = Int_use_ch1(MeanStartindex:MeanEndindex);
     
                    %channel2
        selectTime_ch2Mean = timeVec(MeanStartindex:MeanEndindex);
        selectF0_ch2Mean = F0_use_ch2(MeanStartindex:MeanEndindex);
        selectF1_ch2Mean = F1_use_ch2(MeanStartindex:MeanEndindex);
        selectF2_ch2Mean = F2_use_ch2(MeanStartindex:MeanEndindex);
        selectInt_ch2Mean = Int_use_ch2(MeanStartindex:MeanEndindex);
        
        %saving all means
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
        line([baselineStart baselineStart], ylim, 'color', [1 1 1], 'LineWidth', 1.5) %start of baseline
        line([shiftStart shiftStart], ylim, 'color', [1 1 1], 'LineWidth', 1.5) %start of perturbation
        line([shiftEnd shiftEnd], ylim, 'color', [1 1 1], 'LineWidth', 1.5) %start of perturbation
        
        axes(sub4) %select subplot4
        line([baselineStart baselineStart], ylim, 'color', [1 1 1], 'LineWidth', 1.5) %start of baseline
        line([shiftStart shiftStart], ylim, 'color', [1 1 1], 'LineWidth', 1.5) %start of perturbation
        line([shiftEnd shiftEnd], ylim, 'color', [1 1 1], 'LineWidth', 1.5) %start of perturbation
        
        else %if a non-shifted reflex trial or an adaptation trial
        
            [x_spect,y_spect] = ginput(2); %only doing timing in output
        
        timex_spect = x_spect; %only using "x"
        timey_spect = y_spect;
  
        
        shiftStart = timex_spect(1); %marked start of trial
        shiftEnd = timex_spect(2); %marked end of trial
       
        if scripttype == 1 %adaptation script
        % only for saving full mic trial, not used for analysis
       
          micstart = timex_spect(1)-delaycalc; %shift the mic back the delay time, relative to the headphones. %%only for saving full mic trial, not used for analysis
          [c micStartindex] = min(abs(timeVec - micstart));
          [c Endindex] = min(abs(timeVec - shiftEnd));
          %channel1 for SAVING (shifted)
         saveTime_ch1 = timeVec(micStartindex:Endindex);
        saveF0_ch1 = F0_use_ch1(micStartindex:Endindex);
        saveF1_ch1 = F1_use_ch1(micStartindex:Endindex);
        saveF2_ch1 = F2_use_ch1(micStartindex:Endindex);
        saveInt_ch1 = Int_use_ch1(micStartindex:Endindex);
         end
         
           
             % %find indicies that relate to the above values
        [c Startindex] = min(abs(timeVec - shiftStart)); %find the vector point closests to the calculated start/end time
        [c Endindex] = min(abs(timeVec - shiftEnd));
        
       
         
           % save information from ch 1 and ch 2
        %channel1 for ANALYSIS
        selectTime_ch1 = timeVec(Startindex:Endindex);
        selectF0_ch1 = F0_use_ch1(Startindex:Endindex);
        selectF1_ch1 = F1_use_ch1(Startindex:Endindex);
        selectF2_ch1 = F2_use_ch1(Startindex:Endindex);
        selectInt_ch1 = Int_use_ch1(Startindex:Endindex);
        
        
        
         %channel2 for ANALYSIS AND SAVING
        selectTime_ch2 = timeVec(Startindex:Endindex);
        selectF0_ch2 = F0_use_ch2(Startindex:Endindex);
        selectF1_ch2 = F1_use_ch2(Startindex:Endindex);
        selectF2_ch2 = F2_use_ch2(Startindex:Endindex);
        selectInt_ch2 = Int_use_ch2(Startindex:Endindex);
        
        F0calc_ch1 = mean(selectF0_ch1);
        F1calc_ch1 = mean(selectF1_ch1*1000);
        F2calc_ch1 = mean(selectF2_ch1*1000);
        Intcalc_ch1 = mean(selectInt_ch1);
        
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
        line([shiftStart shiftStart], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
        line([shiftEnd shiftEnd], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
        
        if scripttype == 1
        line([micstart micstart], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
        end
        
        axes(sub4) %select subplot4
        line([shiftStart shiftStart], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
        line([shiftEnd shiftEnd], ylim, 'color', [1 1 1], 'LineWidth', 1.5)
        
        
        %these are not used for these trials, so set them to 0 
        %  lengthreflex = 0; %if not a reflex script, don't change this
           lengthVoiceOnset_Reflex = 0;
           %lengthfullreflex =0;
            baselineF0_ch1_mean = 0;
               baselineF0_ch2_mean = 0;
        
        end
        
    end
    function callbackdelay (PushButtondelay, EventData)
        %% calculate delay at voice onset!
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
        
        
        timeVoiceOnset = ch1_delaystart; 
          axes(sub1) %select subplot3
        line([timeVoiceOnset, timeVoiceOnset], ylim, 'color', [1 0 0], 'LineWidth', 1.5) %plot the voicing onset
         set(PushButtonspect, 'Visible', 'on');

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