% PA Negative Fragmentation Library
% This script will build PA_neg_Matrix.mat file which will contain
% prediction of PA lipid fragmentation 
%-------------------User settings for library generation 
min_FA = 10; % Minimum number of C in fatty acid chain
max_FA = 30; % Maximum number of C in fatty acid chain
double_bonds= 5; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

% PA fragmentation pathways  
polarity = ('-'); 
%Names of fragments 
% metal adducts
fragments_add{1} = {'FA1 carboxylate anion','FA2 carboxylate anion','loss of FA1 as carboxylic acid',...
    'loss of FA2 as carboxylic acid','loss of FA1 as ketene',...
    'loss of FA2 as ketene','PO3','H2PO4','glycerol phosphate - water',...
    'glycerol phosphate'}; % Adduct
fragments_add{2} = {'FA1 carboxylate anion',...
    'loss of FA1 as carboxylic acid','loss of FA1 as ketene',...
    'PO3','H2PO4','glycerol phosphate - water',...
    'glycerol phosphate'}; % Lyso adduct

% = generate FA list 
FA_list = combvec(min_FA:max_FA,0:double_bonds)'; % first is number of carbons, second is number double bonds 

% fragment library masses
PA_group = 97.97689;

% lipid pieaces
glycerol_residue = 74.03678;
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
        P_matrix(m).mass = F1_ketene + F2_ketene + glycerol_residue + PA_group + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};
        
        P_matrix(m).frag(1) = F1_ketene + OH_neg; %F1 carboxylate anion
        P_matrix(m).frag(3) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(5) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(7) = 78.9585; %PO3
        P_matrix(m).frag(8) = 96.9691; %H2PO4
        P_matrix(m).frag(9) = 152.995; %glycerol phosphate - water
        P_matrix(m).frag(10) = 171.006; %glycerol phosphate

        
        
        if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4]))
        P_matrix(m).frag(2) = NaN; %F2 carboxylate anion
        P_matrix(m).frag(4) = NaN; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(6) = NaN; %loss of FA2 as ketene  
        else
        P_matrix(m).frag(2) = F2_ketene + OH_neg; %F2 carboxylate anion
        P_matrix(m).frag(4) = P_matrix(m).mass -F2_ketene - water; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(6) = P_matrix(m).mass -F2_ketene; %loss of FA2 as ketene       
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
        P_matrix(m).mass = F1_ketene + glycerol_residue + PA_group + adduct_mass(o);
        
        P_matrix(m).adduct = sprintf('Lyso %s',adduct{o});
        
        P_matrix(m).frag(1) = F1_ketene + OH_neg; %F1 carboxylate anion
        P_matrix(m).frag(2) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(3) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(4) = 78.9585; %PO3
        P_matrix(m).frag(5) = 96.9691; %H2PO4
        P_matrix(m).frag(6) = 152.995; %glycerol phosphate - water
        P_matrix(m).frag(7) = 171.006; %glycerol phosphate     
        P_matrix(m).frag_name = fragments_add{2};
        m = m+1;      
end

end

lipid_name = ('PA');
save('PA_neg_Matrix.mat','P_matrix','lipid_name','polarity');


