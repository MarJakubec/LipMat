% This is an automated script for further filtering of MS2 hits and plotting
% of the results. The script will load all *MS2_hits.mat files in the working
% folder and filter the results according to their score. User input here is 
% score value which is in the default set to 20. Higher score than 30 indicates
% very reliable identification. Score less than 20 is very unreliable 
% identification. The script will generate *_MS2_identified.xls table with
% a table of identified lipid species, sorted by score. It will also plot 
% fragmentation for all identified lipids. In case of more hits for one lipid 
% species, it will plot only highest scoring lipid species. Figures are 
% numbered based on position in .xls table. Both .xls file and figures will
% be in a new folder with the name of analyzed files. It will also prepare 
% files for futher MS1 analysis. 
% 
% Input: *MS2_hits.mat
% Output: *_for_MS1.mat, *_MS2_identified.xls and lipid fragmentation plots 



%------Begining of user input-------------
cut_off = 20; % Cutt off score, only higher score will be ploted and exported. 
%------End of user input -----------------


% Import _*filt_P_positive.mat files 
files = dir('*MS2_hits.mat');
for o = 1:size(files,1)
load(files(o).name);
% 
% Filtering of MS2 peaks 
MS2_hits = MS2_hits([MS2_hits(:).max_score] >= cut_off);%!Cut off !


% if ~all(isnan(MS2_hits(1).all_hits(:,1)))
if isempty(MS2_hits)
    continue
end

    finds = [{MS2_hits.structure}' {MS2_hits.adduct}' {MS2_hits.lipid}'];
    
    for i = 1:size(finds,1)
        if strcmp(finds{i,3},'SM')
        str_finds{i} = sprintf('%s %s d%d:%d/%d:%d',finds{i,3},finds{i,2},finds{i,1}(1),...
            finds{i,1}(2),finds{i,1}(3),finds{i,1}(4));
        elseif strcmp(finds{i,3},'CL')
        str_finds{i} = sprintf('%s %s %d:%d/%d:%d/%d:%d/%d:%d',finds{i,3},finds{i,2},finds{i,1}(1),...
            finds{i,1}(2),finds{i,1}(3),finds{i,1}(4),finds{i,1}(5),finds{i,1}(6),finds{i,1}(7),finds{i,1}(8));    
        else
        str_finds{i} = sprintf('%s %s %d:%d/%d:%d',finds{i,3},finds{i,2},finds{i,1}(1),...
            finds{i,1}(2),finds{i,1}(3),finds{i,1}(4));
        end
    end
    un_finds = unique(str_finds)'; 
    
    for j = 1:size(un_finds,1)
    un_finds{j,2} = mean([MS2_hits(strcmp(str_finds,un_finds{j,1})).Precursor_peak_mz]); % m/z mean
    un_finds{j,3} = std([MS2_hits(strcmp(str_finds,un_finds{j,1})).Precursor_peak_mz])/un_finds{j,2}*1E6; % ppm of m/z
    un_finds{j,4} = sum(strcmp(str_finds,un_finds{j,1})); % Number of hits 
    un_finds{j,5} = mean([MS2_hits(strcmp(str_finds,un_finds{j,1})).max_score]); % Mean Score 
    un_finds{j,6} = std([MS2_hits(strcmp(str_finds,un_finds{j,1})).max_score]);% STD of score 
    un_finds{j,7} = sum([MS2_hits(strcmp(str_finds,un_finds{j,1})).max_score]); % Sum of score
    un_finds{j,8} = mean([MS2_hits(strcmp(str_finds,un_finds{j,1})).MS1_RT]); % Mean of RT
    un_finds{j,9} = std([MS2_hits(strcmp(str_finds,un_finds{j,1})).MS1_RT]); % standard deviation of RT
    un_finds{j,10} = MS2_hits((strcmp(str_finds,un_finds{j,1}))).lipid; % Lipid headgroup
    end
   
    out = sortrows(un_finds,[4,5,6],'descend');
    out = cell2table(out,'VariableNames',...
        {'Lipid_species' 'Mean_mz' 'ppm_of_mz' 'Number_of_hits' 'Mean_score' 'STD_Score' 'Sum_of_score' ...
        'Mean_of_retention_time' 'STD_of_retention_time' 'Lipid_headgroup'});
    % Adduct quantitative export
    
    % Plot of MS2 plots
        newfolder = sprintf('%s',files(o).name (1:end-13));
        mkdir(newfolder);
        addpath(fullfile(pwd, newfolder)); %current folder is returned by pwd. Combined with newfolder gives full path.           
        writetable(out,sprintf('%s/%s_MS2_identified.xls',fullfile(pwd, newfolder),files(o).name (1:end-12)));
        save(sprintf('%s_for_MS1.mat',files(o).name (1:end-12)),'out')
        
        for i = 1:size(out,1)
            filter = MS2_hits(strcmp(str_finds,out{i,1}));
            to_plot =filter([filter.max_score] == max([filter.max_score]));
            if size(to_plot,2) >1
                to_plot = to_plot([to_plot.MS1_intensity] == max([to_plot.MS1_intensity]));
            end            
            to_plot.hits_frag = unique(to_plot.hits_frag,'rows');
            width = 1; 
            figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
            h{1} = subplot(1,3,[1 2]);
            bar(to_plot.MS2_spectrum(:,1),to_plot.MS2_spectrum(:,2),'FaceColor','k',...
                 'EdgeColor','k','LineWidth',0.5,'BarWidth',width/min(diff(sort(to_plot.MS2_spectrum(:,1)))));
            hold on
            if size(to_plot.hits_frag,1)>1
            bar(to_plot.hits_frag(:,1),to_plot.hits_frag(:,2),'FaceColor','r','EdgeColor','r',...
                'LineWidth',0.5,'BarWidth',width/min(diff(sort(to_plot.hits_frag(:,1)))));
            else
            bar(to_plot.hits_frag(:,1),to_plot.hits_frag(:,2),'FaceColor','r','EdgeColor','r',...
                'LineWidth',0.5);
            end
            
            set(gca, 'YScale', 'log');
            xlabel('m/z');
            ylabel('Intensity'); 
            title(sprintf('%s',cell2mat(out{i,1})))
            set(gca,'FontSize', 16);
                        for j = 1:size(to_plot.hits_frag,1)
                            txt = sprintf('%0.3f',to_plot.hits_frag(j,1));
                            text((to_plot.hits_frag(j,1)-8),(to_plot.hits_frag(j,2)+500),txt,'FontSize', 10, 'Rotation',90)
                        end      

            set(h{1}, 'Position', [0.05, 0.100, 0.6, 0.85]);


              c = ismember(round(to_plot.frag,1), round(to_plot.hits_frag(:,1),1));

              to_print = to_plot.frag_name(c);
              to_print_mz = to_plot.frag(c);
              dim = [.68 .4 .5 .5];
              str = [sprintf('%s %.3f m/z',cell2mat(out{i,1}), to_plot.predict_mass)];
            for  j = 1:size(to_print,2)
                str = [str,sprintf('\n %s = %.3f ',to_print{j},to_print_mz(j))];
            end
            annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12)

print(sprintf('%s/LIPID_ID_%d.png',fullfile(pwd, newfolder),i),'-dpng');
       close all
        end
    
    
    
    
    

     
 clearvars -except o files   
end
clear all 
close all 
['No errors!']

