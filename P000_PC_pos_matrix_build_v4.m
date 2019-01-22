% PC Positive Fragmentation Library
% This script will build PC_pos_Matrix.mat file which will contain
% prediction of PC lipid fragmentation 
%-------------------User settings for library generation 
min_FA = 10; % Minimum number of C in fatty acid chain
max_FA = 30; % Maximum number of C in fatty acid chain
double_bonds= 5; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

% Build PC scanning matrix 
polarity = ('+'); 
non_add = ('+H');
%Names of fragments 
% non metal 
fragments_non{1} = {'choline','phosphatidylcholine head-group','R1CO','R2CO',...
    'loss of phosphocholine and FA1 as ketene','loss of phosphocholine and FA2 as ketene',...
    'loss of phosphocholine and FA1 as carboxylic acid','loss of phosphocholine and FA2 as carboxylic acid',...
    'loss of FA1 as carboxylic acid ', 'loss of FA2 as carboxylic acid ',...
    'loss of FA1 as ketene','loss of FA2 as ketene',...
    'loss of phosphocholine', 'loss of NC3H9'}; % +H 
%lyso non metal 
fragments_non{2} = {'choline','phosphatidylcholine head-group','R1CO',...
    'loss of phosphocholine and FA1 as ketene', 'loss of phosphocholine and FA1 as carboxylic acid',...
    'loss of FA1 as carboxylic acid ',...
    'loss of FA1 as ketene','loss of phosphocholine', 'loss of NC3H9'};% Lyso H
% non metal dduct fragments name
fragments_no_metal_add{1} = {'choline','PO4C2H5 & adduct','R1CO','R2CO',...
    'loss of phosphocholine and FA1 as ketene','loss of phosphocholine and FA2 as ketene',...
    'loss of phosphocholine and FA1 as carboxylic acid ','loss of phosphocholine and FA2 as carboxylic acid',...
    'loss of FA1 as ketene','loss of FA2 as ketene','loss of FA1 as carboxylic acid',...
    'loss of FA2 as carboxylic acid','loss of phosphocholine',' loss of NC3H9'}; 
fragments_no_metal_add{2} = {'choline','PO4C2H5 & adduct','R1CO',...
    'loss of phosphocholine and FA1 as ketene',...
    'loss of phosphocholine and FA1 as carboxylic acid ',...
    'loss of FA1 as ketene','loss of FA1 as carboxylic acid',...
    'loss of phosphocholine',' loss of NC3H9'}; 
% metal 
fragments_add{1} = {'choline','phoshocholine & adduct','R1CO','R2CO',...
    'loss of phosphocholine and FA1 as ketene','loss of phosphocholine and FA2 as ketene',...
    'loss of phosphocholine and FA1 as carboxylic acid','loss of phospcholine and FA2 as carboxylic acid',...
    'loss of FA1 as ketene','loss of FA2 as ketene','loss of FA1 as carboxylic acid','loss of FA2 as carboxylic acid',...
    'loss of NC3H9 and FA1 as carboxylic acid','loss of NC3H9 and FA2 as carboxylic acid',...
    'loss of adduct and FA1 as carboxylic acid','loss of adduct and FA2 as carboxylic acid',...
    'loss of choline and FA1 as carboxylic acid','loss of choline and FA2 as carboxylic acid',...
    'loss of choline, adduct and FA1 as carboxylic acid','loss of choline, adduct and FA2 as carboxylic acid',...
    'loss of phosphocholine; loss of adduct','loss of phosphocholine','loss of NC3H9'}; % Adduct
fragments_add{2} = {'choline','PO4C2H5 & adduct','R1CO',...
    'loss of phosphocholine, adduct and FA1 as ketene',...
    'loss of FA1 as carboxylic acid','loss of FA1 as ketene',...
    'loss of phosphocholine', 'loss of NC3H9','loss of phosphocholine and adduct',...
    'loss of phosphocholine and FA1 as carboxylic acid',...
    'loss of NC3H9 and FA1 as carboxylic acid',...
    'loss of FA1 as carboxylic acid andadduct'}; % Lyso adduct

% = generate FA list 
FA_list = combvec(min_FA:max_FA,0:double_bonds)'; % first is number of carbons, second is number double bonds 

% fragment library masses
PC_group = 183.06604; 
glycerol_residue = 74.03678;
water = 18.01056;	
NC3H9 = 59.07350;
PO4C2H5 = 123.992548;
choline = 85.08915;	
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
H_pos = 1.00728;	
Na_pos = 22.98922;	
K_pos = 38.96316;	
Li_pos = 7.01545;	
NH4_pos = 18.03383;

% Building search matrix for H+ 
m = 1;
k= 0;
for i = 1:(size(FA_list,1))
    for j = 1:(size(FA_list,1)-k)
        P_matrix(m).structure(:,[1 2]) = FA_list(i,:);
        P_matrix(m).structure(:,[3 4]) = FA_list(j+k,:);
        F1_ketene = (C12_atom*P_matrix(m).structure(1) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2)) -2));
        F2_ketene = (C12_atom*P_matrix(m).structure(3) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(3)-P_matrix(m).structure(4)) -2));
        P_matrix(m).mass = F1_ketene + F2_ketene + glycerol_residue + PC_group + H_pos;
        P_matrix(m).adduct = non_add;      
        P_matrix(m).frag(1) = choline+H_pos; %choline
        P_matrix(m).frag(2) = PC_group+H_pos; %phosphatidylcholine head-group
        P_matrix(m).frag(3) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(5) = F2_ketene + glycerol_residue +H_pos ; %loss of phosphocholine; loss of FA1 as ketene
        P_matrix(m).frag(7) = F2_ketene + glycerol_residue - water +H_pos ; %loss of phosphocholine; loss of FA1 as carboxylic acid
        P_matrix(m).frag(9) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid 
        P_matrix(m).frag(11) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(13) = P_matrix(m).mass - PC_group; %loss of phosphocholine
        P_matrix(m).frag(14) = P_matrix(m).mass - NC3H9; % loss of NC3H9'
        if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4])) 
        P_matrix(m).frag(4) = NaN; %R2CO
        P_matrix(m).frag(6) = NaN; %loss of phosphocholine; loss of FA2 as ketene
        P_matrix(m).frag(8) = NaN; %loss of phosphocholine; loss of FA2 as carboxylic acid
        P_matrix(m).frag(10) = NaN; %loss of FA2 as carboxylic acid 
        P_matrix(m).frag(12) = NaN;% loss of FA2 as ketene
        else
        P_matrix(m).frag(4) = F2_ketene + H_pos; %R2CO
        P_matrix(m).frag(6) = F1_ketene + glycerol_residue +H_pos; %loss of phosphocholine; loss of FA2 as ketene
        P_matrix(m).frag(8) = F1_ketene + glycerol_residue - water +H_pos; %loss of phosphocholine; loss of FA2 as carboxylic acid
        P_matrix(m).frag(10) = P_matrix(m).mass - F2_ketene - water	; %loss of FA2 as carboxylic acid 
        P_matrix(m).frag(12) = P_matrix(m).mass - F2_ketene; %loss of FA2 as ketene
        end
        P_matrix(m).frag_name = fragments_non{1};
        m = m +1; 
    end
    k = k +1;
end
clear F1_ketene F2_ketene
%  Continue building search matrix for Lyso H+  Cannot delete m variable as
%  it continues looping 
% Building search matrix for Lyso-H+


for i = 1:(size(FA_list,1))
  
        P_matrix(m).structure(:,[1 2]) = FA_list(i,:);
        P_matrix(m).structure(:,[3 4]) = [0 0];
        F1_ketene = (C12_atom*P_matrix(m).structure(1) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2)) -2));
        P_matrix(m).mass = F1_ketene + glycerol_residue + PC_group + H_pos;
        P_matrix(m).adduct = sprintf('Lyso %s',non_add);
        
        P_matrix(m).frag(1) = choline + H_pos; %choline
        P_matrix(m).frag(2) = PC_group + H_pos; %phoshocholine head-group
        P_matrix(m).frag(3) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(4) = glycerol_residue + H_pos; %loss of phosphocholine; loss of FA1 as ketene
        P_matrix(m).frag(5) = glycerol_residue + H_pos - water; %loss of phosphocholine; loss of FA1 as carboxylic acid
        P_matrix(m).frag(6) = P_matrix(m).mass - F1_ketene - water; %loss of FA1 as carboxylic acid 
        P_matrix(m).frag(7) = P_matrix(m).mass - F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(8) = P_matrix(m).mass - PC_group; %loss of phosphocholine
        P_matrix(m).frag(9) = P_matrix(m).mass - NC3H9; % loss of NC3H9'
        P_matrix(m).frag_name = fragments_non{2};       
        m = m +1;
       
end
clear F1_ketene
% Continue adding non metal adducts 
adduct = {'+TEAC','+NH4'}; 
adduct_mass =[130.159024 18.034374];

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
        
        P_matrix(m).frag(1) = choline + H_pos; %choline
        P_matrix(m).frag(2) = PO4C2H5+adduct_mass(o); %PO4C2H5 & adduct
        P_matrix(m).frag(3) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(5) = F2_ketene + glycerol_residue +H_pos; %loss of phosphocholine; loss of FA1 as ketene
        P_matrix(m).frag(7) = F2_ketene + glycerol_residue -water +H_pos; %loss of phosphocholine;,loss of FA1 as carboxylic acid 
        P_matrix(m).frag(9) = P_matrix(m).mass - F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(11) = P_matrix(m).mass - F1_ketene -water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(13) = P_matrix(m).mass - PC_group; %loss of phosphocholine
        P_matrix(m).frag(14) = P_matrix(m).mass - NC3H9; % loss of NC3H9'
        
        if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4]))
        P_matrix(m).frag(4) = NaN;  %R2CO
        P_matrix(m).frag(6) = NaN;  %loss of phosphocholine; loss of FA2 as ketene
        P_matrix(m).frag(8) = NaN;  %loss of phosphocholine;,loss of FA2 as carboxylic acid 
        P_matrix(m).frag(10) = NaN; %loss of FA2 as ketene
        P_matrix(m).frag(12) = NaN; %loss of FA2 as carboxylic acid
    
  
        else
        P_matrix(m).frag(4) = F2_ketene + H_pos;  %R2CO
        P_matrix(m).frag(6) = F1_ketene + glycerol_residue +H_pos;  %loss of phosphocholine; loss of FA2 as ketene
        P_matrix(m).frag(8) = F1_ketene + glycerol_residue -water +H_pos;  %loss of phosphocholine;,loss of FA2 as carboxylic acid 
        P_matrix(m).frag(10) = P_matrix(m).mass - F2_ketene; %loss of FA2 as ketene
        P_matrix(m).frag(12) = P_matrix(m).mass - F1_ketene -water; %loss of FA2 as carboxylic acid
  
        end
        P_matrix(m).frag_name = fragments_no_metal_add{1};
        m = m +1; 
    end
    k = k +1;
end

% Building search matrix for Lyso-non metal adduct+
for i = 1:(size(FA_list,1))
        P_matrix(m).structure(:,[1 2]) = FA_list(i,:);
        P_matrix(m).structure(:,[3 4]) = [0 0];
        
        F1_ketene = (C12_atom*P_matrix(m).structure(1) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2)) -2));
        P_matrix(m).mass = F1_ketene + glycerol_residue + PC_group + adduct_mass(o);
        P_matrix(m).adduct = sprintf('Lyso %s',adduct{o});
        
        P_matrix(m).frag(1) = choline + H_pos; %choline
        P_matrix(m).frag(2) = PO4C2H5+adduct_mass(o); %PO4C2H5 & adduct
        P_matrix(m).frag(3) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(4) = glycerol_residue + H_pos; %loss of phosphocholine; loss of FA1 as ketene
        P_matrix(m).frag(5) = P_matrix(m).mass - F1_ketene -water; %loss of FA1 as carboxylic acid 
        P_matrix(m).frag(6) = P_matrix(m).mass - F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(7) = P_matrix(m).mass - PC_group; %loss of phosphocholine
        P_matrix(m).frag(8) = P_matrix(m).mass - NC3H9; % loss of NC3H9'
        P_matrix(m).frag_name = fragments_no_metal_add{2};
        m = m+1;      
end

end

% Continue adding non metal adducts 
adduct = {'+Na','+K','+Li'}; 
adduct_mass =[22.98922 38.96316 7.01545];

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
        
        P_matrix(m).frag(1) = choline + H_pos; %choline
        P_matrix(m).frag(2) = PC_group+adduct_mass(o); %phosphocholine
        P_matrix(m).frag(3) = F1_ketene+ H_pos; %R1CO
        P_matrix(m).frag(5) = F2_ketene + glycerol_residue + adduct_mass(o); %loss of phosphocholine ; loss of FA1 as ketene
        P_matrix(m).frag(7) = F2_ketene + glycerol_residue - water + adduct_mass(o); % loss of phosphocholine; loss of FA1 as carboxylic acid 
        P_matrix(m).frag(9) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(11) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(13) = P_matrix(m).mass -F1_ketene - water - NC3H9; %loss of NC3H9 and FA1 as carboxylic acid
        P_matrix(m).frag(15) = P_matrix(m).mass -F1_ketene - water - adduct_mass(o) + H_pos; %loss of adduct and FA1 as carboxylic acid
        P_matrix(m).frag(17) = P_matrix(m).mass -F1_ketene - water - choline; %loss of choline and FA1 as carboxylic acid
        P_matrix(m).frag(19) = P_matrix(m).mass -F1_ketene - water - choline - adduct_mass(o)+ H_pos; %loss of choline, adduct and FA1 as carboxylic acid
        P_matrix(m).frag(21) = P_matrix(m).mass - PC_group - adduct_mass(o) + H_pos; %loss of phosphocholine; loss of adduct
        P_matrix(m).frag(22) = P_matrix(m).mass - PC_group; %loss of phosphocholine
        P_matrix(m).frag(23) = P_matrix(m).mass - NC3H9; % loss of NC3H9'        
        
        if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4]))
        P_matrix(m).frag(4) = NaN;  %R2CO
        P_matrix(m).frag(6) = NaN;  %loss of phosphocholine ; loss of FA2 as ketene
        P_matrix(m).frag(8) = NaN;  % loss of phosphocholine; loss of FA2 as carboxylic acid 
        P_matrix(m).frag(10) = NaN; % loss of FA2 as ketene
        P_matrix(m).frag(12) = NaN; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(14) = NaN; %loss of NC3H9 and FA2 as carboxylic acid
        P_matrix(m).frag(16) = NaN; %loss of adduct and FA2 as carboxylic acid
        P_matrix(m).frag(18) = NaN; %loss of choline and FA2 as carboxylic acid
        P_matrix(m).frag(20) = NaN;%loss of choline, adduct and FA2 as carboxylic acid
  
        else
        P_matrix(m).frag(4) = F2_ketene + H_pos; %R2CO
        P_matrix(m).frag(6) = F1_ketene + glycerol_residue + adduct_mass(o); %loss of phosphocholine ; loss of FA2 as ketene
        P_matrix(m).frag(8) = F1_ketene + glycerol_residue - water + adduct_mass(o); % loss of phosphocholine; loss of FA2 as carboxylic acid 
        P_matrix(m).frag(10) = P_matrix(m).mass -F2_ketene; %loss of FA2 as ketene
        P_matrix(m).frag(12) = P_matrix(m).mass -F2_ketene - water; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(14) = P_matrix(m).mass -F2_ketene - water - NC3H9; %loss of NC3H9 and FA2 as carboxylic acid
        P_matrix(m).frag(16) = P_matrix(m).mass -F2_ketene - water - adduct_mass(o) + H_pos; %loss of adduct and FA2 as carboxylic acid
        P_matrix(m).frag(18) = P_matrix(m).mass -F2_ketene - water - choline; %loss of choline and FA2 as carboxylic acid
        P_matrix(m).frag(20) = P_matrix(m).mass -F2_ketene - water - choline - adduct_mass(o)+ H_pos; %loss of choline, adduct and FA2 as carboxylic acid
       
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
        
        P_matrix(m).frag(1) = choline + H_pos; %choline
        P_matrix(m).frag(2) = PO4C2H5 + adduct_mass(o); %PO4C2H5 & adduct
        P_matrix(m).frag(3) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(4) = glycerol_residue + H_pos; %loss of phosphocholine & adduct; loss of FA2 as ketene
        P_matrix(m).frag(5) = P_matrix(m).mass - F1_ketene - water; %loss of FA1 as carboxylic acid 
        P_matrix(m).frag(6) = P_matrix(m).mass - F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(7) = P_matrix(m).mass - PC_group; %loss of phosphocholine
        P_matrix(m).frag(8) = P_matrix(m).mass - NC3H9; % loss of NC3H9'
        P_matrix(m).frag(9) = P_matrix(m).mass - PC_group - adduct_mass(o) + H_pos; %loss of phosphocholine; loss of adduct
        P_matrix(m).frag(10) = P_matrix(m).mass - PC_group -F1_ketene - water;%       loss of phosphocholine; loss of FA2 as carboxylic acid'
        P_matrix(m).frag(11) =  P_matrix(m).mass - NC3H9 - F1_ketene - water;%     'loss of NC3H9; loss of FA2 as carboxylic acid'
        P_matrix(m).frag(12) = P_matrix(m).mass - F1_ketene - water - adduct_mass(o) + H_pos;%     'loss of FA2 as carboxylic acid; loss of adduct'
        P_matrix(m).frag_name = fragments_add{2};
        m = m+1;      
end

end

lipid_name = ('PC');
save('PC_pos_Matrix.mat','P_matrix','lipid_name','polarity');


