% PI Negative Fragmentation Library
% This script will build PI_neg_Matrix.mat file which will contain
% prediction of PI lipid fragmentation 
%-------------------User settings for library generation 
min_FA = 10; % Minimum number of C in fatty acid chain
max_FA = 30; % Maximum number of C in fatty acid chain
double_bonds= 5; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!


% PI fragmentation pathways  
polarity = ('-'); 
%Names of fragments 

% metal adducts
fragments_add{1} = {'FA1 carboxylate anion','FA2 carboxylate anion','glycerophosphoinositol anion',...
    'glycerophosphoinositol anion - water','glycerophosphoinositol anion - 2x water',...
    'loss of inositol and FA1 as carboxylic acid','loss of inositol and FA2 as carboxylic acid',...
    'loss of inositol and FA1 as ketene','loss of inositol and FA2 as ketene',...
    'loss of FA1 as carboxylic acid','loss of FA2 as carboxylic acid','loss of FA1 as ketene',...
    'loss of FA2 as ketene','PO3','H2PO4','glycerol phosphate - water',...
    'glycerol phosphate','phosphoinositol - 2x water','phosphoinositol - water',...
    'phosphoinositol','glycerophosphoinositol - 3x water'}; % Adduct
fragments_add{2} = {'FA1 carboxylate anion','glycerophosphoinositol anion',...
    'glycerophosphoinositol anion - water','glycerophosphoinositol anion - 2x water',...
    'loss of inositol and FA1 as carboxylic acid','loss of inositol and FA1 as ketene',...
    'loss of FA1 as carboxylic acid','loss of FA1 as ketene',...
    'PO3','H2PO4','glycerol phosphate - water',...
    'glycerol phosphate','phosphoinositol - 2x water','phosphoinositol - water',...
    'phosphoinositol','glycerophosphoinositol - 3x water'}; % Lyso adduct

% = generate FA list 
FA_list = combvec(min_FA:max_FA,0:double_bonds)'; % first is number of carbons, second is number double bonds 

% fragment library masses
PI_group = 260.02971;
glycerol_residue = 74.03678;	
inositol = 162.05282;	
glycerophosphoinositol_neg = 333.05921;	
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
        P_matrix(m).mass = F1_ketene + F2_ketene + glycerol_residue + PI_group + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};
        
        P_matrix(m).frag(1) = F1_ketene + OH_neg; %F1 carboxylate anion
        P_matrix(m).frag(3) = glycerophosphoinositol_neg; %glycerophosphoinositol 
        P_matrix(m).frag(4) = glycerophosphoinositol_neg - water; %glycerophosphoinositol
        P_matrix(m).frag(5) = glycerophosphoinositol_neg -2*water; %glycerophosphoinositol
        P_matrix(m).frag(6) = P_matrix(m).mass - F1_ketene - water - inositol; %loss of inositol and FA1 as carboxylic acid
        P_matrix(m).frag(8) = P_matrix(m).mass - F1_ketene - inositol; % loss of inositol and FA1 as ketene
        P_matrix(m).frag(10) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(12) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(14) = 78.9585; %PO3
        P_matrix(m).frag(15) = 96.9691; %H2PO4
        P_matrix(m).frag(16) = 152.995; %G3P - H2O
        P_matrix(m).frag(17) = 171.006; %G3P 
        P_matrix(m).frag(18) = 223.001; %phosphoinositol - 2H2O
        P_matrix(m).frag(19) = 241.011; %phosphoinositol - H2O
        P_matrix(m).frag(20) = 259.022; %phosphoinositol 
        P_matrix(m).frag(21) = 279.027; %glycerophosphoinositol - 3 water        
         if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4]))
        P_matrix(m).frag(2) = NaN; %F2 carboxylate anion
        P_matrix(m).frag(7) = NaN; %loss of inositol and FA2 as carboxylic acid
        P_matrix(m).frag(9) = NaN; % loss of inositol and FA2 as ketene
        P_matrix(m).frag(11) = NaN; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(13) = NaN; %loss of FA2 as ketene      
        else
        P_matrix(m).frag(2) = F2_ketene + OH_neg; %F2 carboxylate anion
        P_matrix(m).frag(7) = P_matrix(m).mass - F2_ketene - water - inositol; %loss of inositol and FA2 as carboxylic acid
        P_matrix(m).frag(9) = P_matrix(m).mass - F2_ketene - inositol; % loss of inositol and FA2 as ketene
        P_matrix(m).frag(11) = P_matrix(m).mass -F2_ketene - water; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(13) = P_matrix(m).mass -F2_ketene; %loss of FA2 as ketene    
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
        P_matrix(m).mass = F1_ketene + glycerol_residue + PI_group + adduct_mass(o);
        
        P_matrix(m).adduct = sprintf('Lyso %s',adduct{o});
        
        P_matrix(m).frag(1) = F1_ketene + OH_neg; %F1 carboxylate anion
        P_matrix(m).frag(2) = glycerophosphoinositol_neg; %glycerophosphoinositol 
        P_matrix(m).frag(3) = glycerophosphoinositol_neg - water; %glycerophosphoinositol
        P_matrix(m).frag(4) = glycerophosphoinositol_neg -2*water; %glycerophosphoinositol
        P_matrix(m).frag(5) = P_matrix(m).mass - F1_ketene - water - inositol; %loss of inositol and FA1 as carboxylic acid
        P_matrix(m).frag(6) = P_matrix(m).mass - F1_ketene - inositol; % loss of inositol and FA1 as ketene
        P_matrix(m).frag(7) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(8) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(9) = 78.9585; %PO3
        P_matrix(m).frag(10) = 96.9691; %H2PO4
        P_matrix(m).frag(11) = 152.995; %G3P - H2O
        P_matrix(m).frag(12) = 171.006; %G3P 
        P_matrix(m).frag(13) = 223.001; %phosphoinositol - 2H2O
        P_matrix(m).frag(14) = 241.011; %phosphoinositol - H2O
        P_matrix(m).frag(15) = 259.022; %phosphoinositol 
        P_matrix(m).frag(16) = 279.027; %glycerophosphoinositol - 3 water        

        P_matrix(m).frag_name = fragments_add{2};
        m = m+1;      
end

end

lipid_name = ('PI');
save('PI_neg_Matrix.mat','P_matrix','lipid_name','polarity');


