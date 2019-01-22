% PC Negative Fragmentation Library
% This script will build PC_neg_Matrix.mat file which will contain
% prediction of PC lipid fragmentation 
%-------------------User settings for library generation 
min_FA = 10; % Minimum number of C in fatty acid chain
max_FA = 30; % Maximum number of C in fatty acid chain
double_bonds= 5; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

% PC fragmentation pathways  
polarity = ('-'); 
%Names of fragments 
% metal adducts
fragments_add{1} = {'FA1 carboxylate anion','FA2 carboxylate anion','loss of C5H12N; loss of FA1 as carboxylic acid',...
    'loss of C5H12N; loss of FA2 as carboxylic acid','loss of C5H12N; loss of FA1 as ketene','loss of C5H12N; loss of FA2 as ketene',...
    'loss of C3H10N; loss of FA1 as carboxylic acid','loss of C3H10N; loss of FA2 as carboxylic acid',...
    'loss of C3H10N; loss of FA1 as ketene','loss of C3H10N; loss of FA2 as ketene',...
    'demethylation; loss of FA1 as carboxylic acid','demethylation; loss of FA2 as carboxylic acid',...
    'demethylation; loss of FA1 as ketene','demethylation; loss of FA2 as ketene',...
    'loss of FA1 as carboxylic acid','loss of FA2 as carboxylic acid','loss of FA1 as ketene',...
    'loss of FA2 as ketene'}; % Adduct
fragments_add{2} = {'FA1 carboxylate anion','loss of C5H12N; loss of FA1 as carboxylic acid',...
    'loss of C5H12N; loss of FA1 as ketene',...
    'loss of C3H10N; loss of FA1 as carboxylic acid',...
    'loss of C3H10N; loss of FA1 as ketene',...
    'demethylation; loss of FA1 as carboxylic acid',...
    'demethylation; loss of FA1 as ketene',...
    'loss of FA1 as carboxylic acid','loss of FA1 as ketene'}; % Lyso adduct

% = generate FA list 
FA_list = combvec(min_FA:max_FA,0:double_bonds)'; % first is number of carbons, second is number double bonds 

% fragment library masses

PC_group = 183.06604;


% lipid pieaces
glycerol_residue = 74.03678;
neck = 135.99254;				%% is piece of head, is it?
ethane = 26.01565;				%% -2xH in fact
demethylation = 45.05785;		%%	hmmm...
lyso_demethylation = 72.037; %% loss during demethylation 
water = 18.01056;
H1_atom = 1.00783;
O16_atom = 15.99491;
C12_atom = 12.00000;
N14_atom = 14.00307;
P31_atom = 30.97376;
Li7_atom = 7.01600;
Na23_atom = 22.98977;
K39_atom = 38.96371;
Cl35_atom = 34.96885;
loss_C5H12N = 143.11;
loss_C3H10N = 117.095;
% % Ions and adducts:
% % 
deprotanation_neg = -1.00728;	
Cl_35_neg = 34.96940;	
formate_neg = 44.99820;	
acetate_neg = 59.01385;
OH_neg = 17.00329;	


% Building search matrix for H+ 
m = 1;
k= 0;

% Continue adding non metal adducts 
adduct = {'-H','-OH','-Cl','-CHO2','-CH3CO2'}; 
adduct_mass =[-1.00728 17.00329 34.96940 44.99820 59.01385];

for o = 1: size(adduct,2)
k= 0;

for i = 1:(size(FA_list,1))
    for j = 1:(size(FA_list,1)-k)
        P_matrix(m).structure(:,[1 2]) = FA_list(i,:);
        P_matrix(m).structure(:,[3 4]) = FA_list(j+k,:);  
        F1_ketene = (C12_atom*P_matrix(m).structure(1) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2)) -2));
        F2_ketene = (C12_atom*P_matrix(m).structure(3) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(3)-P_matrix(m).structure(4)) -2));
        P_matrix(m).mass = F1_ketene + F2_ketene + glycerol_residue + PC_group + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};
        
        P_matrix(m).frag(1) = F1_ketene + OH_neg; %F1 carboxylate anion
        P_matrix(m).frag(3) = F2_ketene + neck -water + OH_neg; %loss of C5H12N; loss of FA1 as carboxylic acid
        P_matrix(m).frag(5) = F2_ketene + neck + OH_neg; % loss of C5H12N; loss of FA1 as ketene
        P_matrix(m).frag(7) = F2_ketene + neck + ethane - water + OH_neg; %loss of C3H10N; loss of FA1 as carboxylic acid
        P_matrix(m).frag(9) = F2_ketene + neck + ethane + OH_neg;%loss of C3H10N; loss of FA1 as ketene
        P_matrix(m).frag(11) = F2_ketene + neck + ethane + demethylation - water + OH_neg; %demethylation loss of FA1 as carboxylic acid
        P_matrix(m).frag(13) = F2_ketene + neck + ethane + demethylation + OH_neg; %demethylation loss of FA1 as ketene
        P_matrix(m).frag(15) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(17) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
         if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4]))
        P_matrix(m).frag(2) = NaN; %F2 carboxylate anion
        P_matrix(m).frag(4) = NaN; %loss of C5H12N; loss of FA2 as carboxylic acid
        P_matrix(m).frag(6) = NaN; % loss of C5H12N; loss of FA2 as ketene
        P_matrix(m).frag(8) = NaN; %loss of C3H10N; loss of FA2 as carboxylic acid
        P_matrix(m).frag(10) = NaN;%loss of C3H10N; loss of FA2 as ketene
        P_matrix(m).frag(12) = NaN; %demethylation loss of FA2 as carboxylic acid
        P_matrix(m).frag(14) = NaN; %demethylation loss of FA2 as ketene
        P_matrix(m).frag(16) = NaN; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(18) = NaN; %loss of FA2 as ketene      
        else
        P_matrix(m).frag(2) = F2_ketene + OH_neg; %F2 carboxylate anion
        P_matrix(m).frag(4) = F1_ketene + neck -water + OH_neg; %loss of C5H12N; loss of FA2 as carboxylic acid
        P_matrix(m).frag(6) = F1_ketene + neck + OH_neg; % loss of C5H12N; loss of FA2 as ketene
        P_matrix(m).frag(8) = F1_ketene + neck + ethane - water + OH_neg; %loss of C3H10N; loss of FA2 as carboxylic acid
        P_matrix(m).frag(10) = F1_ketene + neck + ethane + OH_neg;%loss of C3H10N; loss of FA2 as ketene
        P_matrix(m).frag(12) = F1_ketene + neck + ethane + demethylation - water + OH_neg; %demethylation loss of FA2 as carboxylic acid
        P_matrix(m).frag(14) = F1_ketene + neck + ethane + demethylation + OH_neg; %demethylation loss of FA2 as ketene
        P_matrix(m).frag(16) = P_matrix(m).mass -F2_ketene - water; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(18) = P_matrix(m).mass -F2_ketene; %loss of FA2 as ketene       
        end
        P_matrix(m).frag_name = fragments_add{1};
        m = m +1; 
    end
    k = k +1;
end

% Building search matrix for Lyso-Metal adducts+
for i = 1:(size(FA_list,1))
        P_matrix(m).structure(:,[1 2]) = FA_list(i,:);
        P_matrix(m).structure(:,[3 4]) = [0 0];
        F1_ketene = (C12_atom*P_matrix(m).structure(1) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2)) -2));
        P_matrix(m).mass = F1_ketene + glycerol_residue + PC_group + adduct_mass(o);
        
        P_matrix(m).adduct = sprintf('Lyso %s',adduct{o});
        
        P_matrix(m).frag(1) = F1_ketene + OH_neg; %F1 carboxylate anion
        P_matrix(m).frag(2) = P_matrix(m).mass - F1_ketene -water -loss_C5H12N -adduct_mass(o) + OH_neg; %loss of C5H12N; loss of FA1 as carboxylic acid
        P_matrix(m).frag(3) = P_matrix(m).mass - F1_ketene -loss_C5H12N -adduct_mass(o) + OH_neg; % loss of C5H12N; loss of FA1 as ketene
        P_matrix(m).frag(4) = P_matrix(m).mass - F1_ketene -water -loss_C3H10N -adduct_mass(o) + OH_neg; %loss of C3H10N; loss of FA1 as carboxylic acid
        P_matrix(m).frag(5) = P_matrix(m).mass - F1_ketene -loss_C3H10N -adduct_mass(o) + OH_neg;%loss of C3H10N; loss of FA1 as ketene
        P_matrix(m).frag(6) = P_matrix(m).mass - F1_ketene -water -lyso_demethylation  -adduct_mass(o) + OH_neg; %demethylation loss of FA1 as carboxylic acid
        P_matrix(m).frag(7) = P_matrix(m).mass - F1_ketene -lyso_demethylation -adduct_mass(o) + OH_neg; %demethylation loss of FA1 as ketene
        P_matrix(m).frag(8) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(9) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag_name = fragments_add{2};
        m = m+1;      
end

end

lipid_name = ('PC');
save('PC_neg_Matrix.mat','P_matrix','lipid_name','polarity');


