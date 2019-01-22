% This is an automated script for generation of ms1.mat and ms2.mat files 
% from input mzXML files. raw files need to be converted in mzXML files 
% using, for example, MSConvert (http://proteowizard.sourceforge.net/tools.shtml).
% When converting files make sure to untick zlib compression. Matlab is unable to process compressed mzXML files. 
% 
% No need to add any additional user information, only run this 
% script in the same folder as the all mzXML files which should be converted. 
% 
% Input: non-compressed mzXML files
% Output: ms1.mat and ms2.mat for each original mzXML files



% Import of results MS2
results = dir('*.mzXML');
for o = 1: length(results)
    tic
    file = results(o).name;
    
    %MS 2 export starts here
    mzXMLstruct = mzxmlread(file,'Levels',2);


% Export Ms2 files with precuros mass and retention time 
MS2 = cell(size(mzXMLstruct.scan,1),7);
for k = 1:size(mzXMLstruct.scan,1)
    MSR = zeros(length(mzXMLstruct.scan(k).peaks.mz)/2,2);
    for i = 1:(length(mzXMLstruct.scan(k).peaks.mz)/2)
    MSR(i,1)= mzXMLstruct.scan(k).peaks.mz((i*2)-1,1);
    MSR(i,2)= mzXMLstruct.scan(k).peaks.mz(i*2,1);
    end

    MS2{k,1}(:,1) = MSR(:,1); 
    MS2{k,1}(:,2) = MSR(:,2);
% %     MS2{k,1}(:,2) = log10(MSR(:,2)); % for log transofmration 
% %     MS2{k,1} = mspeaks(MSR(:,1),MSN(:,1),'DENOISING',false,'HEIGHTFILTER',0.1); % Peak picking DONT DO IT
    MS2{k,2} = mzXMLstruct.scan(k).num; % Number of scan
    MS2{k,3} = mzXMLstruct.scan(k).precursorMz.value; % Precursor peak
    MS2{k,4} = mzXMLstruct.scan(k).precursorMz.precursorScanNum; % Precursor scan number
    MS2{k,5} = (str2double(mzXMLstruct.scan(k).retentionTime(3:length(mzXMLstruct.scan(k).retentionTime)-1)))./60; % Retention time in minutes
    MS2{k,6} = mzXMLstruct.scan(k).precursorMz.precursorIntensity; % Precursor intenzity 
    MS2{k,7} = mzXMLstruct.scan(k).polarity; % Polarity of scan 
    clear MSR
end
 %Export
    readouts_MS2 = char('spectrum','scan number','precursor peak value','precursor scan number','precursor retention time','precursor intensity','polarity');
    save(sprintf('%s.ms2.mat',results(o).name(1:end-6)),'MS2','readouts_MS2')
    clearvars -except results o file
    %MS2 expoert ends here
    
    %MS1 export starts here
       mzXMLstruct = mzxmlread(file,'Levels',1); % MS level during import 
    
% Export Ms1 files with precuros mass and retention time 
 MS1 = cell(size(mzXMLstruct.scan,1),4);
for k = 1:size(mzXMLstruct.scan,1)    
    MSR = zeros(length(mzXMLstruct.scan(k).peaks.mz)/2,2);
    for i = 1:(length(mzXMLstruct.scan(k).peaks.mz)/2)
    MSR(i,1)= mzXMLstruct.scan(k).peaks.mz((i*2)-1,1);
    MSR(i,2)= mzXMLstruct.scan(k).peaks.mz(i*2,1);
    end
    MSP = mspeaks(MSR(:,1),MSR(:,2),'DENOISING',false); % Peak picking
    MSP = sortrows(MSP,2,'descend'); %Sort by intensity 
    if size(MSP,1) >= 1000 
    MS1{k,1} = sortrows(MSP(1:1000,:),1,'ascend'); % Filter 1000 most intense peaks and sort by m/z
    else
    MS1{k,1} = sortrows(MSP,1,'ascend'); % If there is lower number of peaks take them all    
    end
    MS1{k,2} = mzXMLstruct.scan(k).num; % Scan number
    MS1{k,3} = (str2double(mzXMLstruct.scan(k).retentionTime(3:length(mzXMLstruct.scan(k).retentionTime)-1)))./60; % Retention time in minures
    MS1{k,4} = mzXMLstruct.scan(k).polarity; % Polarity of scan 
    clear MSR   
end
  MS1_read = char('spectra','scan number','retention time','polarity');
  save(sprintf('%s.ms1.mat',results(o).name(1:end-6)),'MS1','MS1_read')
  clearvars -except results o
 toc   
    
end
disp('No errors!')

