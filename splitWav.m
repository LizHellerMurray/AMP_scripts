function splitWav(fnames, FolderName) 
%this function saves the audapter signal in and signal out for the F1
%reflex analysis function

for j = 1:length(in) %creates loop, for 1 until the last file
    filename = strcat(FolderName,'\',in{j}); %creates full file name using path and mat name
    load(filename) %loads mat file
    ch1filename = strcat(FolderName,'\mic','\',in{j}); % makes new file path to "mic" folder
    ch1filename = strrep(ch1filename, 'mat','wav'); %changes file name from .mat to .wav 
    y = data.signalIn; %y is signal in data
    Fs = 16000; %sampling frequency
    audiowrite(ch1filename,y,Fs) %creates audio file
end

%same thing but for audapter signal out
for i = 1:length(in) %creates loop, for 1 until the last file
    filename = strcat(FolderName,'\',in{i}); 
    load(filename)
    ch2filename = strcat(FolderName,'\perturbed_mic','\',in{i});
    ch2filename = strrep(ch2filename, 'mat','wav');
    y = data.signalOut;
    Fs = 16000;
    audiowrite(ch2filename,y,Fs)
end
end