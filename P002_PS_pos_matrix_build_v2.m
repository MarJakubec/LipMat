% PS Positive Fragmentation Library
% This script will build PS_pos_Matrix.mat file which will contain
% prediction of PS lipid fragmentation 
%-------------------User settings for library generation 
min_FA = 10; % Minimum number of C in fatty acid chain
max_FA = 30; % Maximum number of C in fatty acid chain
double_bonds= 5; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

% PS fragmentation pathways  
polarity = ('+'); 
non_add = ('+H');
%Names of fragments 
% non +H  
fragments_non{1} = {'serine phosphate','R1CO','R2CO',...
    'loss of FA1 as ketene and head group','loss of FA2 as ketene and head group',...
    'loss of FA1 as ketene','loss of FA2 as ketene',...
    'loss of FA1 as carboxylic acid ', 'loss of FA2 as carboxylic acid',...
    'loss of phosphoserine'}; % +H 
%lyso + H
fragments_non{2} = {'serine phosphate','R1CO',...
    'loss of FA1 as ketene and head group',...
    'loss of FA1 as ketene',...
    'loss of FA1 as carboxylic acid','loss of phosphoserine'};% Lyso H
% non metal dduct fragments name NH4 and TEAC
fragments_no_metal_add{1} = {'FA1 as ketene','FA2 as ketene','loss of FA1 as ketene and head group',...
    'loss of FA2 as ketene and head group','loss of FA1 as ketene','loss of FA2 as ketene',...
    'loss of FA1 as carboxylic acid','loss of FA2 as carboxylic acid',...
    'phosphoserine head-group','loss of phosphoserine'}; 
fragments_no_metal_add{2} = {'FA1 as ketene','loss of FA1 as ketene and head group',...
    'loss of FA1 as ketene','loss of FA1 as carboxylic acid',...
    'phosphoserine head-group','loss of phosphoserine'}; 
% metal adducts
fragments_add{1} = {'serine phosphate','FA1 as ketene','FA2 as ketene','loss of FA1 as ketene, adduct and head group',...
    'loss of FA2 as ketene, adduct and head group','loss of FA1 as carboxylic acid, adduct and head group',...
    'loss of FA2 as carboxylic acid, adduct and head group','loss of FA1 as ketene','loss of FA2 as ketene',...
    'loss of FA1 as carboxylic acid','loss of FA2 as carboxylic acid','loss of serine and FA1 as carboxylic acid',...
    'loss of serine and FA2 as carboxylic acid','loss of serine and FA1 as ketene',...
    'loss of serine and FA2 as ketene','loss of phosphoserine and adduct',...
    'loss of phosphoserine','loss of phosphoserine, adittion of water',...
    'loss of phosphoric acid','loss of serine'}; % Adduct
fragments_add{2} = {'serine phosphate','FA1 as ketene','loss of FA1 as ketene and head group',...
   'loss of FA1 as carboxylic acid and head group',...
    'loss of FA1 as ketene','loss of FA1 as carboxylic acid','loss of serine and FA1 as carboxylic acid',...
    'loss of serine and FA1 as ketene',...
    'loss of phosphoserine and adduct',...
    'loss of phosphoserine','loss of phosphoserine, adittion of water',...
    'loss of phosphoric acid','loss of serine'}; % Lyso adduct

% = generate FA list 
FA_list = combvec(min_FA:max_FA,0:double_bonds)'; % first is number of carbons, second is number double bonds 

% fragment library masses
PS_group = 185.00892;

% lipid pieaces
glycerol_residue = 74.03678;
serine = 87.03203;		
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
        P_matrix(m).mass = F1_ketene + F2_ketene + glycerol_residue + PS_group + H_pos;
        P_matrix(m).adduct = non_add;      
        
        P_matrix(m).frag(1) = PS_group+H_pos; %serine phosphate
        P_matrix(m).frag(2) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(4) = F2_ketene + glycerol_residue +H_pos ; %loss of FA1 as ketene and head group
        P_matrix(m).frag(6) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(8) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(10) = P_matrix(m).mass -PS_group; %loss of phosphoserine

        if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4])) 
        P_matrix(m).frag(3) = NaN; %R2CO
        P_matrix(m).frag(5) = NaN; %loss of FA2 as ketene and head group
        P_matrix(m).frag(7) = NaN; %loss of FA2 as ketene
        P_matrix(m).frag(9) = NaN; %loss of FA2 as carboxylic acid
        else
        P_matrix(m).frag(3) = F2_ketene + H_pos; %R2CO
        P_matrix(m).frag(5) = F1_ketene + glycerol_residue +H_pos ; %loss of FA2 as ketene and head group
        P_matrix(m).frag(7) = P_matrix(m).mass -F2_ketene; %loss of FA2 as ketene
        P_matrix(m).frag(9) = P_matrix(m).mass -F2_ketene - water; %loss of FA2 as carboxylic acid
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
        P_matrix(m).mass = F1_ketene + glycerol_residue + PS_group + H_pos;
        P_matrix(m).adduct = sprintf('Lyso %s',non_add);
        
        P_matrix(m).frag(1) = PS_group+H_pos; %serine phosphate
        P_matrix(m).frag(2) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(3) = glycerol_residue +H_pos ; %loss of FA1 as ketene and head group
        P_matrix(m).frag(4) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(5) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(6) = P_matrix(m).mass -PS_group; %loss of phosphoserine
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
        P_matrix(m).mass = F1_ketene + F2_ketene + glycerol_residue + PS_group + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};        
        P_matrix(m).frag(1) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(3) = F2_ketene + glycerol_residue +H_pos; %loss of FA1 as ketene and head group
        P_matrix(m).frag(5) = P_matrix(m).mass - F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(7) = P_matrix(m).mass - F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(9) = PS_group + H_pos; % phosphoserine head-group
        P_matrix(m).frag(10) = P_matrix(m).mass - PS_group; % loss of phosphoserine
        if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4]))
        P_matrix(m).frag(2) = NaN;  %R2CO
        P_matrix(m).frag(4) = NaN;  %loss of FA2 as ketene and head group
        P_matrix(m).frag(6) = NaN;  %%loss of FA2 as ketene 
        P_matrix(m).frag(8) = NaN; %loss of FA2 as carboxylic acid
        else
        P_matrix(m).frag(2) = F2_ketene + H_pos; %R2CO
        P_matrix(m).frag(4) = F1_ketene + glycerol_residue +H_pos; %loss of FA2 as ketene and head group
        P_matrix(m).frag(6) = P_matrix(m).mass - F2_ketene; %loss of FA2 as ketene
        P_matrix(m).frag(8) = P_matrix(m).mass - F2_ketene - water; %loss of FA2 as carboxylic acid
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
        P_matrix(m).mass = F1_ketene + glycerol_residue + PS_group + adduct_mass(o);
        P_matrix(m).adduct = sprintf('Lyso %s',adduct{o});
        
        P_matrix(m).frag(1) = F1_ketene + H_pos; %R1CO
        P_matrix(m).frag(2) = glycerol_residue +H_pos ; %loss of FA1 as ketene and head group
        P_matrix(m).frag(3) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(4) = P_matrix(m).mass - F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(5) = PS_group + H_pos; %phosphoserine
        P_matrix(m).frag(6) = P_matrix(m).mass - PS_group; %loss of phosphoserine

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
        P_matrix(m).mass = F1_ketene + F2_ketene + glycerol_residue + PS_group + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};
        
        P_matrix(m).frag(1) = PS_group + adduct_mass(o); %serine phosphate
        P_matrix(m).frag(2) = F1_ketene+H_pos; %R1CO
        P_matrix(m).frag(4) = P_matrix(m).mass -F1_ketene - PS_group- adduct_mass(o) + H_pos; % loss of FA1 as ketene and head group
        P_matrix(m).frag(6) = P_matrix(m).mass -F1_ketene - PS_group - water - adduct_mass(o) + H_pos; % loss of FA1 as carboxylic acid and head group
        P_matrix(m).frag(8) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(10) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(12) = P_matrix(m).mass -F1_ketene - water - serine; % loss of serine and FA1 as carboxylic acid
        P_matrix(m).frag(14) = P_matrix(m).mass -F1_ketene - serine; % loss of serine and FA1 as ketene
        P_matrix(m).frag(16) = P_matrix(m).mass -PS_group - adduct_mass(o) + H_pos; % loss of phosphoserine, loss of adduct 
        P_matrix(m).frag(17) = P_matrix(m).mass - PS_group; % loss of phosphoserine
        P_matrix(m).frag(18) = P_matrix(m).mass - PS_group + water; % loss phosphoserine, addition of water
        P_matrix(m).frag(19) = P_matrix(m).mass -97.977 ; % loss of phosphoric acid
        P_matrix(m).frag(20) = P_matrix(m).mass - serine; % loss of serine 
             
        if all(P_matrix(m).structure(:,[1 2]) == P_matrix(m).structure(:,[3 4]))
        P_matrix(m).frag(3) = NaN; %R2CO
        P_matrix(m).frag(5) = NaN; % loss of FA2 as ketene and head group
        P_matrix(m).frag(7) = NaN; % loss of FA2 as carboxylic acid and head group
        P_matrix(m).frag(9) = NaN; %loss of FA2 as ketene
        P_matrix(m).frag(11) = NaN; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(13) = NaN; % loss of ethylamine and FA2 as carboxylic acid
        P_matrix(m).frag(15) = NaN; %loss serine; loss of FA2 as ketene        
        else
        P_matrix(m).frag(3) = F2_ketene+H_pos; %R2CO
        P_matrix(m).frag(5) = P_matrix(m).mass -F2_ketene - PS_group - adduct_mass(o) + H_pos; % loss of FA2 as ketene and head group
        P_matrix(m).frag(7) = P_matrix(m).mass -F2_ketene - PS_group - water - adduct_mass(o) + H_pos; % loss of FA2 as carboxylic acid and head group
        P_matrix(m).frag(9) = P_matrix(m).mass -F2_ketene; %loss of FA2 as ketene
        P_matrix(m).frag(11) = P_matrix(m).mass -F2_ketene - water; %loss of FA2 as carboxylic acid
        P_matrix(m).frag(13) = P_matrix(m).mass -F2_ketene - water - serine; % loss of serine and FA2 as carboxylic acid
        P_matrix(m).frag(15) = P_matrix(m).mass -F2_ketene - serine; %loss of serine and FA2 as ketene
              
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
        P_matrix(m).mass = F1_ketene + glycerol_residue + PS_group + adduct_mass(o);
        
        P_matrix(m).adduct = sprintf('Lyso %s',adduct{o});
        
        P_matrix(m).frag(1) = PS_group + adduct_mass(o); %serine phosphate
        P_matrix(m).frag(2) = F1_ketene+H_pos; %R1CO
        P_matrix(m).frag(3) = P_matrix(m).mass -F1_ketene - PS_group - adduct_mass(o) + H_pos; % loss of FA1 as ketene and head group
        P_matrix(m).frag(4) = P_matrix(m).mass -F1_ketene - PS_group - water - adduct_mass(o) + H_pos; % loss of FA1 as carboxylic acid and head group
        P_matrix(m).frag(5) = P_matrix(m).mass -F1_ketene; %loss of FA1 as ketene
        P_matrix(m).frag(6) = P_matrix(m).mass -F1_ketene - water; %loss of FA1 as carboxylic acid
        P_matrix(m).frag(7) = P_matrix(m).mass -F1_ketene - water - serine; % loss of serine and FA1 as carboxylic acid
        P_matrix(m).frag(8) = P_matrix(m).mass -F1_ketene - serine ; %loss of serine,  and FA1 as ketene
        P_matrix(m).frag(9) = P_matrix(m).mass -PS_group - adduct_mass(o) + H_pos; % loss of phosphoserine, loss of adduct 
        P_matrix(m).frag(10) = P_matrix(m).mass - PS_group; % loss of phosphoserine
        P_matrix(m).frag(11) = P_matrix(m).mass - PS_group + water; % loss phosphoserine, addition of water
        P_matrix(m).frag(12) = P_matrix(m).mass -97.977 ; % loss of phosphoric acid
        P_matrix(m).frag(13) = P_matrix(m).mass - serine; % loss of serine 

        P_matrix(m).frag_name = fragments_add{2};
        m = m+1;      
end

end

lipid_name = ('PS');
save('PS_pos_Matrix.mat','P_matrix','lipid_name','polarity');


