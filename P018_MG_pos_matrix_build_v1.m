% MG Positive Fragmentation Library
% This script will build MG_pos_Matrix.mat file which will contain
% prediction of MG lipid fragmentation 

clear all
close all
%-------------------User settings for library generation 
min_FA = 10; % Minimum number of C in fatty acid chain
max_FA = 30; % Maximum number of C in fatty acid chain
double_bonds= 5; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

polarity = ('+'); 

%Names of fragments 
% protonated, + non metal  
fragments_non{1} = {'Precursor ion [M+H]+ with loss of H2O',...
    'FA1 as ketene',...
    'FA1 as ketene with loss of H2O',...
    'FA1 as carboxylic acid'}; 

FA_list = combvec(min_FA:max_FA,0:double_bonds)'; % first is number of carbons, second is number double bonds 

% Atoms:
H1_atom = 1.007825;
O16_atom = 15.994915;
C12_atom = 12.00000;
N14_atom = 14.003074;
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
tg_head = 3*O16_atom + 3*C12_atom + 5*H1_atom;
NH3 = N14_atom + 3*H1_atom;
% Ions and adducts:
H_pos = 1.00728;
Na_pos = 22.98922;
K_pos = 38.96316;
Li_pos = 7.01545;
NH4_pos = 18.03383;


adduct = {'+H','+NH4'}; 
adduct_mass =[1.00728 18.03383];

% Building search matrix for non metal adducts 
m = 1;
for o = 1: size(adduct,2)

for i = 1:(size(FA_list,1))
    
        P_matrix(m).structure(:,[1:2]) = FA_list(i,:);
        FA1_ketene = (C12_atom*P_matrix(m).structure(1) + (1*O16_atom) +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2)) -2));

        
        P_matrix(m).mass = FA1_ketene + glycerol_residue+ water + adduct_mass(o);

        P_matrix(m).adduct = adduct{o};      
        
        P_matrix(m).frag(1) = P_matrix(m).mass - adduct_mass(o) - water + H_pos; %'Precursor ion [M+H]+ with loss of H2O'
        
        P_matrix(m).frag(2) = FA1_ketene + H_pos; %FA1 as ketene
        
        P_matrix(m).frag(3) = FA1_ketene - water + H_pos; %FA1 as ketene with loss of water
        
        P_matrix(m).frag(4) = FA1_ketene + H_pos + water; %FA as carboxylic acid
        
        P_matrix(m).frag_name = fragments_non{1};
        
        %Replace repeating values with NaN
        %prevents doubling score for same mz fragments
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
        clear FA1_ketene 
end

end


lipid_name = ('MG');
save('MG_pos_Matrix.mat','P_matrix','lipid_name','polarity');


