% CL Positive Fragmentation Library
% This script will build CL_pos_Matrix.mat file which will contain
% prediction of CL lipid fragmentation 
% Warning!!! CL have four fatty acid chains, high number of combination
% will lead to exponentional increase of time needed for library creation. 

%-------------------User settings for library generation 
min_FA = 14; % Minimum number of C in fatty acid chain
max_FA = 24; % Maximum number of C in fatty acid chain
double_bonds= 1; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

polarity = ('+'); 

%Names of fragments 
% protonated, + non metal  
fragments_non{1} = {'FA1 as ketene','FA2 as ketene',...
    'FA3 as ketene','FA4 as ketene','FA1 + glycerol',...
    'FA2 + glycerol','FA3 + glycerol','FA4 + glycerol',...
    'FA1 FA2 + glycerol','FA3 FA4 + glycerol'}; 

% metal 
fragments_add{1} = {'FA1 as ketene','FA2 as ketene',...
    'FA3 as ketene','FA4 as ketene','FA1 + glycerol',...
    'FA2 + glycerol','FA3 + glycerol','FA4 + glycerol',...
    'FA1 FA2 + glycerol','FA3 FA4 + glycerol',...
    'FA1 FA2 glycerol phosphate','FA3 FA4 glycerol phosphate',...
    'loss of FA1 as carboxylic acid','loss of FA2 as carboxylic acid',...
    'loss of FA3 as carboxylic acid','loss of FA4 as carboxylic acid',...
    'loss of FA1 as ketene','loss of FA2 as ketene',...
    'loss of FA3 as ketene','loss of FA4 as ketene'}; % Adduct


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


adduct = {'+H','+NH4','+TEAC','-H+2NH4','-H+2TEAC'}; 
adduct_mass =[1.00728 18.03383 130.159024 35.060373 259.310773];

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
        
        P_matrix(m).mass = FA1 + FA2 + FA3 +FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 4*hydroxyl - 3*water + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};      
        P_matrix(m).frag(1) = FA1 - H1_atom + H_pos; %FA1 as ketene
        P_matrix(m).frag(2) = FA2 - H1_atom + H_pos; %FA2 as ketene
        P_matrix(m).frag(3) = FA3 - H1_atom + H_pos; %FA3 as ketene
        P_matrix(m).frag(4) = FA4 - H1_atom + H_pos; %FA4 as ketene
        P_matrix(m).frag(5) = FA1 - H1_atom + glycerol_residue + H_pos; %lFA1 and glycerol
        P_matrix(m).frag(6) = FA2 - H1_atom + glycerol_residue + H_pos; %lFA2 and glycerol
        P_matrix(m).frag(7) = FA3 - H1_atom + glycerol_residue + H_pos; %lFA3 and glycerol
        P_matrix(m).frag(8) = FA4 - H1_atom + glycerol_residue + H_pos; %lFA4 and glycerol
        P_matrix(m).frag(9) = FA1 + FA2 - 2*H1_atom + glycerol_residue + H_pos; %FA1 FA2 and glycerol
        P_matrix(m).frag(10) = FA3 + FA4 - 2*H1_atom + glycerol_residue + H_pos; %FA3 FA4 and glycerol
        
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
        clear temp tempx tempy ci cii n bin
        clear FA1 FA2 FA3 FA4
end

end

% Continue building matrix for metal adducts

adduct = {'+Na','+K','+Li','-H+2Na','-H+2K','-H+2Li'}; 
adduct_mass =[22.98922 38.96316 7.01545 44.971165 76.919041 13.023635];


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
        
        P_matrix(m).mass = FA1 + FA2 + FA3 +FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 4*hydroxyl - 3*water + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};      
        P_matrix(m).frag(1) = FA1 - H1_atom + H_pos; %FA1 as ketene
        P_matrix(m).frag(2) = FA2 - H1_atom + H_pos; %FA2 as ketene
        P_matrix(m).frag(3) = FA3 - H1_atom + H_pos; %FA3 as ketene
        P_matrix(m).frag(4) = FA4 - H1_atom + H_pos; %FA4 as ketene
        
        P_matrix(m).frag(5) = FA1 - H1_atom + glycerol_residue + H_pos; %lFA1 and glycerol
        P_matrix(m).frag(6) = FA2 - H1_atom + glycerol_residue + H_pos; %lFA2 and glycerol
        P_matrix(m).frag(7) = FA3 - H1_atom + glycerol_residue + H_pos; %lFA3 and glycerol
        P_matrix(m).frag(8) = FA4 - H1_atom + glycerol_residue + H_pos; %lFA4 and glycerol
        
        P_matrix(m).frag(9) = FA1 + FA2 - 2*H1_atom + glycerol_residue + H_pos; %FA1 FA2 and glycerol
        P_matrix(m).frag(10) = FA3 + FA4 - 2*H1_atom + glycerol_residue + H_pos; %FA3 FA4 and glycerol
        
        P_matrix(m).frag(11) = FA1 + FA2 -2*H1_atom + water + 2*glycerol_residue + 2*phosphate_dehydro + adduct_mass(o); %FA1 FA2 glycerol phosphate
        P_matrix(m).frag(12) = FA3 + FA4 -2*H1_atom + water + 2*glycerol_residue + 2*phosphate_dehydro + adduct_mass(o); %FA3 FA4 glycerol phosphate
        
        P_matrix(m).frag(13) = FA2 + FA3 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 3*water + adduct_mass(o); %loss of FA1 as carboxylic acid 
        P_matrix(m).frag(14) = FA1 + FA3 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 3*water + adduct_mass(o); %loss of FA2 as carboxylic acid 
        P_matrix(m).frag(15) = FA1 + FA2 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 3*water + adduct_mass(o); %loss of FA3 as carboxylic acid 
        P_matrix(m).frag(16) = FA1 + FA2 + FA3 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 3*water + adduct_mass(o); %loss of FA4 as carboxylic acid 
        
        P_matrix(m).frag(17) = FA2 + FA3 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 2*water + adduct_mass(o); %loss of FA1 as ketene
        P_matrix(m).frag(18) = FA1 + FA3 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 2*water + adduct_mass(o); %loss of FA2 as ketene
        P_matrix(m).frag(19) = FA1 + FA2 + FA4 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 2*water + adduct_mass(o); %loss of FA3 as ketene
        P_matrix(m).frag(20) = FA1 + FA2 + FA3 + 2*phosphate_dehydro + 3*glycerol_residue + 3*hydroxyl - 2*water + adduct_mass(o); %loss of FA4 as ketene
        
        
        P_matrix(m).frag_name = fragments_add{1};
        
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
        clear temp tempx tempy ci cii n bin
        clear FA1 FA2 FA3 FA4
end
end

lipid_name = ('CL');
save('CL_pos_Matrix.mat','P_matrix','lipid_name','polarity');


