% This is automated scripts for scoring of MS2 files according to the 
% prepared lipid fragmentation libraries. The script will load all 
% ms2.mat and all lipid libraries (matrix.mat files) present in the 
% working folder and produce the main table with MS2 scoring section. 
% This table will be saved as [NAME_of_FILE]_MS2_hits.mat. In case 
% of no hits, file [NAME_of_FILE]_MS2_NO_hits.mat will be generated. 
% NO_hits.mat is used by the further script as a guideline which files 
% should be skipped, so please don’t delete it. There is no need to 
% keep positive and negative MS run separated as the script is checking
% for compatibility of libraries and files. So its possible to run all
% MS2 files and all libraries simultaneously. For more information about 
% scoring, please read [once published, relevant paper link will be here].
% 
% Input: *ms2.mat and *matrix.mat files
% Output: *MS2_hits.mat or *MS2_NO_hits.mat



% Import or proccesed results
results = dir('*ms2.mat');
matrix_dir = dir('*_Matrix.mat');
ppm1 = 10; %for MS1 scan 
ppm2 = 20; % for MS2 fragments 
colheadings = {'MS2_spectrum' 'MS2_scan' 'Precursor_peak_mz' 'MS1_scan' 'MS1_RT' 'MS1_intensity' 'polarity'};
for o = 1:size(results,1)
load(results(o).name);

m = 0;
for p = 1:size(matrix_dir,1) % Loop for matrix will use all matrix files in folder
    load(matrix_dir(p).name);
    
if ~all([MS2{:,7}] == polarity) % check polarity of tested files and library
fprintf(['Non matching polarities found in ',results(o).name,' and ',matrix_dir(p).name,...
    '\n Continuing with next iteration of loop \n'])
continue 
end


lipid_scan = unique([P_matrix(:).mass]); %find individual exact masses 

for i = 1:size(MS2,1)

    %test with i = 582
    % PA test 2358  
   ind_test = lipid_scan(lipid_scan < MS2{i,3}+MS2{i,3}*ppm1/1E6 &...
      lipid_scan > MS2{i,3}-MS2{i,3}*ppm1/1E6);
  
    if ~isempty(ind_test)
        P_to_test = P_matrix([P_matrix(:).mass] < MS2{i,3}+MS2{i,3}*ppm1/1E6 &...
        [P_matrix(:).mass] > MS2{i,3}-MS2{i,3}*ppm1/1E6);
        for j = 1:size(P_to_test,2)
            hit = [];
            for k = 1:size(P_to_test(j).frag,2)
            temp = MS2{i,1}(MS2{i,1}(:,1) > P_to_test(j).frag(k)-P_to_test(j).frag(k)*ppm2/1E6 &...
            MS2{i,1}(:,1) < P_to_test(j).frag(k)+P_to_test(j).frag(k)*ppm2/1E6,:);
            hit =[hit;temp];
            end
         P_to_test(j).Hits = hit;
         if isempty(P_to_test(j).Hits)
             continue 
         end         
         %Scoring section - Peak score 
            TheoreticalPeaks = size(P_to_test(j).frag,2);
            ExperimentalPeaks = size(MS2{i,1},1);
            Num_Hits = size(P_to_test(j).Hits,1);
%             Bins = 67703;
           P_to_test(j).Peak_Score = -log (hygepdf(Num_Hits,67703,TheoreticalPeaks,ExperimentalPeaks));
           % Intensity score  
            GotHigher = 1;
            intensityAll = MS2{i,1}(:,2)';
%             tic
            if size(intensityAll,2) >= 6
                for k=1:500
            intensityAll = intensityAll(randperm(length(intensityAll)));
            intensityShort = mat2cell(intensityAll, 1, [Num_Hits,length(intensityAll) - Num_Hits]);
            if sum(intensityShort{1}) > sum(P_to_test(j).Hits(:,2))
            GotHigher = GotHigher +1;
            end
                end
                P_to_test(j).Int_Score = -log(GotHigher/500);
            else    
            all_comb = combnk(intensityAll,Num_Hits);
            for k = 1:size(all_comb,1)
            if sum(all_comb(k,:)) > sum(P_to_test(j).Hits(:,2))
            GotHigher = GotHigher +1;
            end
            end
                P_to_test(j).Int_Score = -log(GotHigher/size(all_comb,1));
            end
        % end for intensity score 
%            toc
            %Sum of both score
            P_to_test(j).Score = P_to_test(j).Peak_Score + P_to_test(j).Int_Score;
        end
        P_to_test = P_to_test(~cellfun(@isempty,{P_to_test.Hits})); % delete entries without hits  
        if ~isempty(P_to_test)
            m = m+1;
        temp = cell2struct(MS2(i,:),colheadings,2);
        MS2_hits(m).MS2_spectrum = temp.MS2_spectrum;
        MS2_hits(m).MS2_scan = temp.MS2_scan;
        MS2_hits(m).Precursor_peak_mz = temp.Precursor_peak_mz;
        MS2_hits(m).MS1_scan = temp.MS1_scan;
        MS2_hits(m).MS1_RT = temp.MS1_RT;
        MS2_hits(m).MS1_intensity = temp.MS1_intensity;
        MS2_hits(m).polarity = temp.polarity;
        
        MS2_hits(m).all_hits = P_to_test;
        MS2_hits(m).max_score = max([P_to_test.Score]);
        MS2_hits(m).lipid = lipid_name;
        
        % filter maximum hits 
        temp_N = find([MS2_hits(m).all_hits(:).Score] == max([MS2_hits(m).all_hits(:).Score]));
        MS2_hits(m).structure = MS2_hits(m).all_hits(temp_N).structure ;
        MS2_hits(m).adduct = MS2_hits(m).all_hits(temp_N).adduct ;
        MS2_hits(m).predict_mass = MS2_hits(m).all_hits(temp_N).mass;
        MS2_hits(m).frag = MS2_hits(m).all_hits(temp_N).frag;
        MS2_hits(m).hits_frag = MS2_hits(m).all_hits(temp_N).Hits; 
        MS2_hits(m).frag_name = MS2_hits(m).all_hits(temp_N).frag_name; 
        
        end
    end


end




clearvars -except results o ppm* colheadings MS2_hits matrix_dir m p MS2 readouts_MS2
        
end
     
if exist('MS2_hits','var')
save(sprintf('%s_MS2_hits.mat',(results(o).name(1:size(results(o).name,2)-8))),...
    'MS2_hits','ppm*')
else
    No_hits = {'No hits found!'}
save(sprintf('%s_MS2_NO_hits.mat',(results(o).name(1:size(results(o).name,2)-8))),...
    'No_hits')   
end
   
clearvars -except results o ppm* colheadings  matrix_dir p  

 
end

['No errors!']
