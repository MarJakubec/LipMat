% This is an automated script for MS1 analysis. Script will load all 
% *for_MS1.mat files in the working folder. It will then compare retention
% time of identified species. If two (or more) lipid species have same m/z
% (within ppm error) and retention time (RT + offset), these species will
% be labeled as nondistinguishable, as they elute as one peak. From this 
% list, only the highest scoring lipid is taken for further analysis. 
% Files with all nondistinguishable lipid species (*_Non_distinguishable.csv)
% is generated. Individual lipid species elution areas are then analyzed and
% files *MS1_area.csv is generated. For five most abundant lipid species in
% each headgroup, elution chromatogram is generated as *_CHROMATOGRAM.png 
% file. For each identified phospholipid also bar plot with five most 
% abundant fatty acid is generated. The user can change retention time 
% offset. From our experience, 2 minutes will allow to observe the whole 
% peak in mixed lipid samples However for sharper peaks, the lower 
% value is recommended. 
% 
% Input: *for_MS1.mat and *MS1.mat files
% Output: *_Non_distinguishable.csv, *MS1_area.csv, Chromatogram pplot 
% and fatty acid bar plot for each phospholipid headgroup
% 



%------Begining of user input-------------
ppm = 5; %ppm error for MS1 search
ret_offset = 2; %additional offset for MS1 range in minutes 
%------End of user input -----------------


%Import
files = dir('*for_MS1.mat');

for o = 1: size(files,1) 
    load(files(o).name);
    load(sprintf('%s.ms1.mat',files(o).name(1:end-13)));
    
    newfolder = sprintf('%s',files(o).name (1:end-13));
    addpath(fullfile(pwd, newfolder));
    
    lipids = out.Lipid_headgroup;   
    un_lipids = unique(lipids);
    out_work = table2struct(out);
    
for k = 1:size(un_lipids,1)
    ind = strcmp(out.Lipid_headgroup,un_lipids{k});
    out_temp = out_work(ind);
% Check for non-distinguisable lipids 
same_hits = [];
  m = 0;
  [out_temp(:).Tested] = deal(0);
for i = 1:size(out_temp,1) 
un_test = [out_temp.Mean_mz]' >= out_temp(i).Mean_mz-out_temp(i).Mean_mz*ppm/1E6 &...
[out_temp.Mean_mz]' <= out_temp(i).Mean_mz+out_temp(i).Mean_mz*ppm/1E6 &...
[out_temp.Mean_of_retention_time]' >= out_temp(i).Mean_of_retention_time - out_temp(i).STD_of_retention_time - ret_offset &...
[out_temp.Mean_of_retention_time]' <= out_temp(i).Mean_of_retention_time + out_temp(i).STD_of_retention_time + ret_offset;
 
if sum(un_test) == 1 && out_temp(i).Tested == 0
 m = m+1;
 out_unique(m) = out_temp(i);
elseif sum(un_test) > 1 && out_temp(i).Tested == 0
 m = m+1;
 out_unique(m) = out_temp(i);
 same_hits = [same_hits; out_temp(un_test)]; %Generting output of nondistingusihable species
 [out_temp(un_test).Tested] = deal(1); 
    
end
clear un_test
end

MS1_to_plot = rmfield(out_unique,'Tested');

 for j = 1:size(out_unique,2)
  m = 0;  
   low_ret = out_unique(j).Mean_of_retention_time-out_unique(j).STD_of_retention_time - 0.5;
   high_ret = out_unique(j).Mean_of_retention_time+out_unique(j).STD_of_retention_time + 0.5;
 for i = 1:size(MS1,1)
     if MS1{i,3} > low_ret && MS1{i,3} <high_ret
         m = m+1;
        temp = MS1{i,1}(MS1{i,1}(:,1) > out_unique(j).Mean_mz-out_unique(j).Mean_mz*ppm/1E6 &...
        MS1{i,1}(:,1) < out_unique(j).Mean_mz+out_unique(j).Mean_mz*ppm/1E6,:);       
     MS1_to_plot(j).plot(m,1) = MS1{i,3};
     MS1_to_plot(j).plot(m,2) = sum(temp(:,2));
     end

 end
     MS1_to_plot(j).plot(MS1_to_plot(j).plot(:,2) == 0,:) = [];
     if size(MS1_to_plot(j).plot(:,2),1) >= 11
         
     MS1_to_plot(j).plot(:,3) = abs(sgolayfilt(MS1_to_plot(j).plot(:,2),3,11));
     else
     MS1_to_plot(j).plot(:,3) = MS1_to_plot(j).plot(:,2);
     end
     
     MS1_to_plot(j).area = sum(MS1_to_plot(j).plot(:,3));
 end
 %sort MS1 struct according to most abundant Hits
 cells = struct2cell(MS1_to_plot); 
 sortvals = cells(12,1,:); 
 mat = cell2mat(sortvals); 
 mat = squeeze(mat); 
 [sorted,ix] = sort(mat); 
 MS1_to_plot = MS1_to_plot(flipud(ix)); 
 
 %Plot first five most abundant hits----------------

 if size(MS1_to_plot,2) >= 5
     num_plot = 5;
 else
     num_plot = size(MS1_to_plot,2);
 end
 
     figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
       hold on
     for i = 1:num_plot
       h = subplot(5,1,i);
       plot(MS1_to_plot(i).plot(:,1),MS1_to_plot(i).plot(:,3),'LineWidth',2)
       xlim([0 MS1{end,3}])
       ylim([0 1.1*max(MS1_to_plot(1).plot(:,3))])
       xticks([])
       xpos = get(h,'Position');
       dim = [0.6 xpos(2)-0.1 0.2 0.2];
       text_box = sprintf('%s Area=%.4e',MS1_to_plot(i).Lipid_species,MS1_to_plot(i).area);
       annotation('textbox',dim,'String',text_box,'FitBoxToText','on');
     end
       xticks(0:5:round(MS1{end,3},1))
       ax1 = axes('Position',[0 0 1 1],'Visible','off');
       descr = 'Intensity';
       descr2 ='Time [min]';
       descr3 = sprintf('%d most abundant %s lipids',i,un_lipids{k});
%        axes(ax1); % sets ax1 to current axes
       hy = text(0.1,0.38,descr,'FontSize', 20,'rotation', 90);
       hx = text(0.5,0.05,descr2,'FontSize', 20);
       hw = text(0.5,0.95,descr3,'FontSize', 24);

       print(sprintf('%s/%s_CHROMATOGRAM.png',fullfile(pwd, newfolder),un_lipids{k}),'-dpng');
       close all 
%End of plot of most abundant hits
%Exporting files of MS1 peak area and possible non-distinguishable peaks
rm_fields ={'plot','Sum_of_score','Lipid_headgroup'};
 MS1_export = rmfield(MS1_to_plot,rm_fields);
 MS1_exp_cells = struct2cell(MS1_export); 
 MS1_exp_cells = squeeze(MS1_exp_cells);
 MS1_out = cell2table(MS1_exp_cells','VariableNames',...
        {'Lipid_species' 'Mean_mz' 'ppm_of_mz' 'Number_of_hits' 'Mean_score' 'STD_Score' ...
        'Mean_of_retention_time' 'STD_of_retention_time' 'Area'});
 writetable(MS1_out,sprintf('%s/%s_MS1_Area.csv',fullfile(pwd, newfolder),un_lipids{k}));

 rm_fields ={'Sum_of_score','Lipid_headgroup','Tested'};
 if ~isempty(same_hits)
 MS1_export = rmfield(same_hits,rm_fields);
 MS1_exp_cells = struct2cell(MS1_export); 
 MS1_exp_cells = squeeze(MS1_exp_cells);
 MS1_out = cell2table(MS1_exp_cells','VariableNames',...
        {'Lipid_species' 'Mean_mz' 'ppm_of_mz' 'Number_of_hits' 'Mean_score' 'STD_Score' ...
        'Mean_of_retention_time' 'STD_of_retention_time'});
 writetable(MS1_out,sprintf('%s/%s_Non_distinguishable.csv',fullfile(pwd, newfolder),un_lipids{k}));
 end
 
 %End of exporting files

%Distribution of FA acids
if ~(strcmp(un_lipids{k},'CL') || strcmp(un_lipids{k},'SM'))
    for i = 1:size(MS1_to_plot,2)
 Str = MS1_to_plot(i).Lipid_species;
 Index = strfind(Str, ':');
 MS1_to_plot(i).FA1 = Str(Index(1)-2:Index(1)+1);
 FA1{i} = Str(Index(1)-2:Index(1)+1);
 MS1_to_plot(i).FA2 = Str(Index(2)-2:Index(2)+1);
 FA2{i} = Str(Index(2)-2:Index(2)+1);
    end
elseif strcmp(un_lipids{k},'SM')
    for i = 1:size(MS1_to_plot,2)
 Str = MS1_to_plot(i).Lipid_species;
 Index = strfind(Str, ':');
 MS1_to_plot(i).FA1 = Str(Index(1)-3:Index(1)+1);
 FA1{i} = Str(Index(1)-3:Index(1)+1);
 MS1_to_plot(i).FA2 = Str(Index(2)-2:Index(2)+1);
 FA2{i} = Str(Index(2)-2:Index(2)+1);
    end
elseif strcmp(un_lipids{k},'CL') 
    for i = 1:size(MS1_to_plot,2)
 Str = MS1_to_plot(i).Lipid_species;
 Index = strfind(Str, ':');
 MS1_to_plot(i).FA1 = Str(Index(1)-2:Index(1)+1);
 FA1{i} = Str(Index(1)-2:Index(1)+1);
 MS1_to_plot(i).FA2 = Str(Index(2)-2:Index(2)+1);
 FA2{i} = Str(Index(2)-2:Index(2)+1);
 MS1_to_plot(i).FA3 = Str(Index(3)-2:Index(3)+1);
 FA3{i} = Str(Index(3)-2:Index(3)+1);
 MS1_to_plot(i).FA4 = Str(Index(4)-2:Index(4)+1);
 FA4{i} = Str(Index(4)-2:Index(4)+1); 
    end
    
end
% Unique FA list
if isfield(MS1_to_plot,'FA3') && isfield(MS1_to_plot,'FA4')
FAList = [{MS1_to_plot.FA1} {MS1_to_plot.FA2} {MS1_to_plot.FA3} {MS1_to_plot.FA4}];
else
FAList = [{MS1_to_plot.FA1} {MS1_to_plot.FA2}];
end
FAList = unique(FAList)';
FAList(strcmp(FAList,'/0:0')) = [];

%Sum of area
[FAList{:,2}] = deal(0);
if isfield(MS1_to_plot,'FA3') && isfield(MS1_to_plot,'FA4')
for i = 1:size(FAList,1)
    Index = find(contains(FA1,FAList{i}));
    if ~isempty(Index)
        for j = 1:size(Index,1)
            FAList{i,2} = FAList{i,2} + MS1_to_plot(Index(j)).area/4;          
        end
    end
    Index = find(contains(FA2,FAList{i}));
    if ~isempty(Index)
        for j = 1:size(Index,1)
            FAList{i,2} = FAList{i,2} + MS1_to_plot(Index(j)).area/4;          
        end
    end 
        Index = find(contains(FA3,FAList{i}));
    if ~isempty(Index)
        for j = 1:size(Index,1)
            FAList{i,2} = FAList{i,2} + MS1_to_plot(Index(j)).area/4;          
        end
    end
        Index = find(contains(FA4,FAList{i}));
    if ~isempty(Index)
        for j = 1:size(Index,1)
            FAList{i,2} = FAList{i,2} + MS1_to_plot(Index(j)).area/4;          
        end
    end
end    
else
for i = 1:size(FAList,1)
    Index = find(contains(FA1,FAList{i}));
    if ~isempty(Index)
        for j = 1:size(Index,1)
            FAList{i,2} = FAList{i,2} + MS1_to_plot(Index(j)).area/2;          
        end
    end
    Index = find(contains(FA2,FAList{i}));
    if ~isempty(Index)
        for j = 1:size(Index,1)
            FAList{i,2} = FAList{i,2} + MS1_to_plot(Index(j)).area/2;          
        end
    end 
end
end

log_area = ([FAList{:,2}])';
perc = log_area/sum(log_area)*100;

FAList = [FAList num2cell(perc)];
FAList = flipud(sortrows(FAList,3));

%Bar plot of most common FA chains? 
 if size(FAList,1) >= 5
     num_bar = 5;
 else
     num_bar= size(FAList,1);
 end
 
bar_p = [FAList{1:num_bar,3}];
name = [FAList(1:num_bar,1)];

figure('units','normalized','outerposition',[0 0 0.60 1],'visible','off')
bar(bar_p)
xticklabels(name)
ylabel('Abundance [%]')
title(sprintf('%s Most abundant FA chains',un_lipids{k}))
set(gca,'FontSize',24)
print(sprintf('%s/%s_FA ABUNDANCE.png',fullfile(pwd, newfolder),un_lipids{k}),'-dpng');
close all
%End of bar plot 

%Export of FA abundance
FA_out = cell2table(FAList,'VariableNames',...
        {'FA_Chain' 'Total_area' 'Abundance'});
 writetable(FA_out,sprintf('%s/%s_FA_Abundance.csv',fullfile(pwd, newfolder),un_lipids{k}));

    clearvars -except o files ppm ret_offset to_scan un_lipids k newfolder lipids out_work out MS1
end
    clearvars -except o files ppm ret_offset
end


['No errors!']
    