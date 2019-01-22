% CL Negative Fragmentation Library
% This script will build CL_neg_Matrix.mat file which will contain
% prediction of CL lipid fragmentation 
% Warning!!! CL have four fatty acid chains, high number of combination
% will lead to exponentional increase of time needed for library creation. 

%-------------------User settings for library generation 
min_FA = 14; % Minimum number of C in fatty acid chain
max_FA = 24; % Maximum number of C in fatty acid chain
double_bonds= 1; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

polarity = ('-'); 

%Names of fragments 
% negative fragments  
fragments_non{1} = {'glycerophosphate',...
    'FA1 as acid','FA2 as acid','FA3 as acid','FA4 as acid',...
    'FA1 glycerol phosphate - water','FA2 glycerol phosphate - water',...
    'FA3 glycerol phosphate - water','FA4 glycerol phosphate - water',...
    'FA1 glycerol phosphate','FA2 glycerol phosphate',...
    'FA3 glycerol phosphate','FA4 glycerol phosphate',...
    'FA1 FA2 glycerol phosphate','FA3 FA4 glycerol phosphate',...
    'loss of FA1 FA2 glycerol phosphate','loss of FA3 FA4 glycerol phosphate',...
    'loss of FA1 FA2 glycerol', 'loss of FA3 FA4 glycerol',...
    'loss of FA1','loss of FA2','loss of FA3','loss of FA4'}; 

%FA list generation
FA_list = combvec(min_FA:max_FA,0:double_bonds)'; % first is number of carbons, second is number double bonds 
m = 1;
k = 0;
for i = 1:size(FA_list,1)
    for j =1:size(FA_list,1)-k
        FA_list2(m,[1 2]) = FA_list(i,:);
        FA_list2(m,[3 4]) = FA_list(j+k,:);
        m= m + 1;
    end
    k = k+1;
end

m = 1;
k = 0;
for i = 1:size(FA_list2,1)
    for j =1:size(FA_list2,1)-k
        FA_list3(m,[1 2 3 4]) = FA_list2(i,:);
        FA_list3(m,[5 6 7 8]) = FA_list2(j+k,:);
        m= m + 1;
    end
    k = k+1;
end
FA_list = FA_list3;


% Exact masses
% Atoms:

H1_atom = 1.00783;
O16_atom = 15.99491;
C12_atom = 12.00000;
N14_atom = 14.00307;
P31_atom = 30.97376;
Li7_atom = 7.01600;
Na23_atom = 22.98977;
K39_atom = 38.96371;
Cl35_atom = 34.96885;
electron_neg = 0.00055;

% Lipid pieces:
glycerol_residue = 74.03678;		
phosphate_dehydro = 79.96632;
hydroxyl = 17.00274;
water = 18.01056;

% Ions and adducts:
H_pos = 1.00728;
Na_pos = 22.98922;
K_pos = 38.96316;
Li_pos = 7.01545;
NH4_pos = 18.03383;
Cl_35_neg = 34.96940;
deprotanation_neg = -1.00728; 

formate_neg = 44.99820;
acetate_neg = 59.01385;
OH_neg = 17.00329;

% adduct = {'-H','-Cl','-OH','-OOCH','-OOCCH3'}; 
adduct = {'-H'};
adduct_mass =[-1.00728];

% Building search matrix for non metal adducts 
m = 1;
for o = 1: size(adduct,2)

for i = 1:(size(FA_list,1))
    
        P_matrix(m).structure(:,[1:8]) = FA_list(i,:);
        FA1 = (C12_atom*P_matrix(m).structure(1) + (1*O16_atom) +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2)) -1));
        FA2 = (C12_atom*P_matrix(m).structure(3) + (1*O16_atom) +...
            H1_atom*(2*(P_matrix(m).structure(3)-P_matrix(m).structure(4)) -1));
        FA3 = (C12_atom*P_matrix(m).structure(5) + (1*O16_atom) +...
            H1_atom*(2*(P_matrix(m).structure(5)-P_matrix(m).structure(6)) -1));
        FA4 = (C12_atom*P_matrix(m).structure(7) + (1*O16_atom) +...
            H1_atom*(2*(P_matrix(m).structure(7)-P_matrix(m).structure(8)) -1));
        
        P_matrix(m).mass = FA1 + FA2 + FA3 +FA4 +2*phosphate_dehydro+3*glycerol_residue+4*hydroxyl -3*water + deprotanation_neg;
        P_matrix(m).adduct = adduct{o};  
        
        P_matrix(m).frag(1) = phosphate_dehydro + glycerol_residue + deprotanation_neg; %glycerophosphate
        
        P_matrix(m).frag(2) = FA1 + hydroxyl + deprotanation_neg; %FA1 as acid
        P_matrix(m).frag(3) = FA2 + hydroxyl + deprotanation_neg; %FA2 as acid
        P_matrix(m).frag(4) = FA3 + hydroxyl + deprotanation_neg; %FA3 as acid
        P_matrix(m).frag(5) = FA4 + hydroxyl + deprotanation_neg; %FA4 as acid
        
        P_matrix(m).frag(6) = FA1 +  phosphate_dehydro + hydroxyl + glycerol_residue - water + deprotanation_neg; %FA1 glycerol phosphate - water
        P_matrix(m).frag(7) = FA2 +  phosphate_dehydro + hydroxyl + glycerol_residue - water + deprotanation_neg; %FA2 glycerol phosphate - water
        P_matrix(m).frag(8) = FA3 +  phosphate_dehydro + hydroxyl + glycerol_residue - water + deprotanation_neg; %FA3 glycerol phosphate - water
        P_matrix(m).frag(9) = FA4 +  phosphate_dehydro + hydroxyl + glycerol_residue - water + deprotanation_neg; %FA4 glycerol phosphate - water

        P_matrix(m).frag(10) = FA1 +  phosphate_dehydro + hydroxyl + glycerol_residue + deprotanation_neg; %FA1 glycerol phosphate
        P_matrix(m).frag(11) = FA2 +  phosphate_dehydro + hydroxyl + glycerol_residue + deprotanation_neg; %FA2 glycerol phosphate
        P_matrix(m).frag(12) = FA3 +  phosphate_dehydro + hydroxyl + glycerol_residue + deprotanation_neg; %FA3 glycerol phosphate
        P_matrix(m).frag(13) = FA4 +  phosphate_dehydro + hydroxyl + glycerol_residue + deprotanation_neg; %FA4 glycerol phosphate
        
        P_matrix(m).frag(14) = FA1 + FA2 + 2*hydroxyl + phosphate_dehydro + glycerol_residue - water + deprotanation_neg; %lFA1 FA2 glycerol phosphate
        P_matrix(m).frag(15) = FA3 + FA4 + 2*hydroxyl + phosphate_dehydro + glycerol_residue - water + deprotanation_neg; %lFA3 FA4 glycerol phosphate
        
        P_matrix(m).frag(16) = FA3 + FA4 +2*hydroxyl + phosphate_dehydro + 2*glycerol_residue - 2*water + deprotanation_neg; %loss of FA1 FA2 glycerol phosphate
        P_matrix(m).frag(17) = FA1 + FA2 +2*hydroxyl + phosphate_dehydro + 2*glycerol_residue - 2*water + deprotanation_neg; %loss of FA3 FA4 glycerol phosphate
        
        P_matrix(m).frag(18) = FA3 + FA4 + 2* hydroxyl + 2*phosphate_dehydro+2*glycerol_residue - 2*water + deprotanation_neg; %loss of FA1 FA2 glycerol 
        P_matrix(m).frag(19) = FA1 + FA2 + 2* hydroxyl + 2*phosphate_dehydro+2*glycerol_residue - 2*water + deprotanation_neg; %loss of FA3 FA4 glycerol 
        
        P_matrix(m).frag(20) = FA2 + FA3 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 3*water + deprotanation_neg; %loss of FA1
        P_matrix(m).frag(21) = FA1 + FA3 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 3*water + deprotanation_neg; %loss of FA2
        P_matrix(m).frag(22) = FA1 + FA2 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 3*water + deprotanation_neg; %loss of FA3
        P_matrix(m).frag(23) = FA1 + FA2 + FA3 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 3*water + deprotanation_neg; %loss of FA4
      
        P_matrix(m).frag_name = fragments_non{1};
        
        %Replace repeating values with NaN
        temp =  P_matrix(m).frag; 
        P_matrix(m).frag = [];
        [n, bin] = histc(temp, unique(temp));
        [tempx, ci, cii] = unique(bin,'first');
        for q = 1:size(ci,1)
        tempy(ci(q)) = temp(ci(q));
        end 
        tempy(tempy ==0) = NaN;
        P_matrix(m).frag = tempy;
        
       m = m +1; 
       clear FA1 FA2 FA3 FA4
       clear temp tempx tempy ci cii n bin
end

end


lipid_name = ('CL');
save('CL_neg_Matrix.mat','P_matrix','lipid_name','polarity');


