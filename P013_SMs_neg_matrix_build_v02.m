% SMs Negative Fragmentation Library
% This script will build SM_neg_Matrix.mat file which will contain
% prediction of SMs lipid fragmentation 
%-------------------User settings for library generation 
min_FA = 10; % Minimum number of C in fatty acid chain
max_FA = 30; % Maximum number of C in fatty acid chain
double_bonds= 5; % maximum number of double bonds in one fatty acid chain; 0 = no double bonds
%-------------------End of are for user setting----------

%Further changes only for skilled!!!

% Build SMS 
% for 
polarity = ('-'); 

%Names of fragments 
% non  adducts including protonation  
fragments_non{1} = {'loss of methyl group and adduct',...
    'loss of FA and loss of methyl group','loss of headgroup and FA',...
    'FA fragment','C4H11NO4P fragment','loss of C5H12N and adduct',...
    'loss of C3H10N and adduct'};

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
loss_C5H12N = (5*C12_atom) +(12*H1_atom)+N14_atom;
loss_C3H10N = (3*C12_atom) +(10*H1_atom)+N14_atom;
% % Ions and adducts:
% % 
deprotanation_neg = -1.00728;

% Building search matrix for nonal  
m = 1;

adduct = {'-CH3CO2','-CHO2'}; 
adduct_mass =[59.013864,44.998214];

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
        P_matrix(m).frag(1) = P_matrix(m).mass- adduct_mass(o) - 15.0234738 +deprotanation_neg + H1_atom; %'loss of methyl group and adduct'
        P_matrix(m).frag(2) = P_matrix(m).mass- adduct_mass(o) - FA_chain - 15.0234738 + deprotanation_neg ; %'loss of FA and methyl group
        P_matrix(m).frag(3) = P_matrix(m).mass- adduct_mass(o) - FA_chain - PC_group + deprotanation_neg + H1_atom ; %'loss of headgroup and FA
        P_matrix(m).frag(4) = FA_chain + deprotanation_neg + (2*H1_atom); %FA fragment
        P_matrix(m).frag(5) = 168.043117937 ; %C4H11NO4P fragment'
        P_matrix(m).frag(6) = P_matrix(m).mass- adduct_mass(o)- loss_C5H12N + deprotanation_neg + H1_atom; %loss of C5H12N and adduct
        P_matrix(m).frag(7) = P_matrix(m).mass - adduct_mass(o) - loss_C3H10N + deprotanation_neg +H1_atom; %loss of C3H10N and adduct
        
        
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


lipid_name = ('SM');
save('SM_neg_Matrix.mat','P_matrix','lipid_name','polarity');


