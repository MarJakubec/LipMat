% SMs Positive Fragmentation Library
% This script will build SM_pos_Matrix.mat file which will contain
% prediction of SMs lipid fragmentation 
%Only 18:1 SM
%-------------------User settings for library generation 
min_FA = 10; % Minimum number of C in fatty acid chain
max_FA = 30; % Maximum number of C in fatty acid chain
double_bonds= 5; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

% Build SMS 
% for 
polarity = ('+'); 
non_add = ('+H');

%Names of fragments 
% non metal adducts including protonation  
fragments_non{1} = {'choline','phosphocholine head-group','fatty acid + NH3',...
    'fatty acid + NC2H3','loss of head-group and loss of FA',...
    'loss of headgroup and FA addition of H2O',...
    'loss of headgroup','loss of headgroup and H2O',...
    'loss of NC3H9'};
% metal adducts
fragments_add{1} = {'choline','phosphocholine head-group',...
    'PO4C2H5 & adduct','fatty acid + NH3','fatty acid + NH3 + adduct'...
    'fatty acid + NC2H3','loss of head-group and FA and adduct',...
    'loss of headgroup and FA; addition of H2O',...
    'loss of headgroup','loss of headgroup and CH2O',...
    'loss of headgroup and adduct and H2O','loss of headgroup and adduct',...
    'loss of headgroup and H2O',...
    'loss of NC3H9'}; % Adduct

% = generate FA list 
FA_list = combvec(min_FA:max_FA,0:double_bonds)'; % first is number of carbons, second is number double bonds 
%generate base list 
min_base = 16;
max_base = 20;
double_bonds = 2;
Base_list = combvec(min_base:max_base,1:double_bonds)';
% fragment library masses
PC_group = 183.06574; 
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

% Building search matrix for nonal  
m = 1;

adduct = {'+H','+NH4','+TEAC'}; 
adduct_mass =[1.00728 18.034374 130.159024];

for o = 1: size(adduct,2)
k= 0;


for i = 1:(size(Base_list,1))
    for j = 1:(size(FA_list,1)-k)
        P_matrix(m).structure(:,[1 2]) = Base_list(i,:);
        P_matrix(m).structure(:,[3 4]) = FA_list(j+k,:);
        Base = (C12_atom*P_matrix(m).structure(1) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2))+1)+ N14_atom); %base as a radical
        FA_chain = (C12_atom*P_matrix(m).structure(3) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(3)-P_matrix(m).structure(4)) -2));
        P_matrix(m).mass = Base + FA_chain + PC_group + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};      
        P_matrix(m).frag(1) = choline+H_pos; %choline
        P_matrix(m).frag(2) = PC_group+H_pos; %phosphocholine head-group
        P_matrix(m).frag(3) = FA_chain + N14_atom +(3*H1_atom) + H_pos; %fatty acid + NH3
        P_matrix(m).frag(4) = FA_chain + N14_atom + (2*C12_atom) +(3*H1_atom) +H_pos ; %fatty acid + NC2H3
        P_matrix(m).frag(5) = Base +H_pos ; %loss of head-group and loss of FA
        P_matrix(m).frag(6) = Base - water + H_pos; %loss of headgroup and FA addition of H2O
        P_matrix(m).frag(7) = P_matrix(m).mass - PC_group; %loss of headgroup
        P_matrix(m).frag(8) = P_matrix(m).mass - PC_group - water; %loss of headgroup and H2O
        P_matrix(m).frag(9) = P_matrix(m).mass - NC3H9; % loss of NC3H9'
        
        P_matrix(m).frag_name = fragments_non{1};
        m = m +1; 
    end
    k = k +1;
end
clear F1_ketene F2_ketene
end
clear adduct adduct_mass
%  Continue building search matrix for metal adducts H+  Cannot delete m variable as
%  it continues looping 


% Continue adding  metal adducts 
adduct = {'+Na','+K','+Li'}; 
adduct_mass =[22.98922 38.96316 7.01545];

for o = 1: size(adduct,2)
k= 0;

for i = 1:(size(Base_list,1))
    for j = 1:(size(FA_list,1)-k)
        P_matrix(m).structure(:,[1 2]) = Base_list(i,:);
        P_matrix(m).structure(:,[3 4]) = FA_list(j+k,:);  
        Base = (C12_atom*P_matrix(m).structure(1) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(1)-P_matrix(m).structure(2))+1)+ N14_atom); %base as a radical
        FA_chain = (C12_atom*P_matrix(m).structure(3) + O16_atom +...
            H1_atom*(2*(P_matrix(m).structure(3)-P_matrix(m).structure(4)) -2));
        P_matrix(m).mass = Base + FA_chain + PC_group + adduct_mass(o);
        P_matrix(m).adduct = adduct{o};      
        P_matrix(m).frag(1) = choline+H_pos; %choline
        P_matrix(m).frag(2) = PC_group+H_pos; %phosphocholine head-group
        P_matrix(m).frag(3) = P31_atom + (4*O16_atom) + (2*C12_atom) + (5*H1_atom) + adduct_mass(o); %PO4C2H5 & adduct
        P_matrix(m).frag(4) = FA_chain + N14_atom +(3*H1_atom) + H_pos; %fatty acid + NH3
        P_matrix(m).frag(5) = FA_chain + N14_atom +(3*H1_atom) + adduct_mass(o); %fatty acid + NH3 + adduct
        P_matrix(m).frag(6) = FA_chain + N14_atom + (2*C12_atom) +(3*H1_atom) +H_pos ; %fatty acid + NC2H3
        P_matrix(m).frag(7) = Base +H_pos ; %loss of head-group and FA and adduct
        P_matrix(m).frag(8) = Base + water + H_pos; %loss of headgroup and FA addition of H2O
        P_matrix(m).frag(9) = P_matrix(m).mass - PC_group; %loss of headgroup
        P_matrix(m).frag(10) = P_matrix(m).mass - PC_group - C12_atom - (2*H1_atom)- O16_atom; %loss of headgroup and CH2O
        P_matrix(m).frag(11) = P_matrix(m).mass - PC_group  - adduct_mass(o)- water + H_pos; %loss of headgroup and adduct and H2O
        P_matrix(m).frag(12) = P_matrix(m).mass - PC_group - adduct_mass(o) + H_pos; %loss of headgroup and adduct
        P_matrix(m).frag(13) = P_matrix(m).mass - PC_group - water; %loss of headgroup and H2O
        P_matrix(m).frag(14) = P_matrix(m).mass - NC3H9; % loss of NC3H9'
        
        P_matrix(m).frag_name = fragments_add{1};
        m = m +1; 
    end
    k = k +1;
end


end

lipid_name = ('SM');
save('SM_pos_Matrix.mat','P_matrix','lipid_name','polarity');


