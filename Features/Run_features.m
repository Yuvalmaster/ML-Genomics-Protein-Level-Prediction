function [known_with_features, known_features, unknown_with_features, unknown_features] = Run_features()
%% loading data
variants = readtable('Variants_sequence.xlsx');
known = readtable('known_data_set.xlsx');
unknown = readtable('unknown_data_set.xlsx');

sequences = variants.VariantSequence;               % extracting the variants from the files

[train_rows,train_cols] = size(known);
[test_rows,test_cols] = size(unknown);

%% Features extraction - Known data set

% Codons Count
[final_struct_avg,avg] = codon_count_mat(sequences,train_rows,known);

% Count bases 
[Afinal,Cfinal,Gfinal,Tfinal] = basecounting(sequences,train_rows,known);

% Find bases dominants between two halfs of the sequence
[binar_first,binar_second] = Dominant(sequences,train_rows,known);
binar_first_names  = {'A Dominant 1', 'C Dominant 1' ,...
                      'G Dominant 1', 'T Dominant 1'};
binar_second_names = {'A Dominant 2', 'C Dominant 2' ,...
                      'G Dominant 2', 'T Dominant 2'};

% PSSM for average of 3 bases
% A dominance = 0 ; C dominance = 1 ; G dominance = 2 ; T dominance = 3
bins = known(:,1:27);
[A_rate, T_rate, C_rate, G_rate, seqs, dominance_seqs, dominance_bins] = PSSM(sequences,bins);
PSSM_names = { '1 - dominant' , '2 - dominant' ,...
               '3 - dominant' , '4 - dominant' ,...
               '5 - dominant' , '6 - dominant' ,...
               '7 - dominant' , '8 - dominant'};
% CpG
CpG_feature = CpG(variants,known);

% Hydrofobic & Polar codons
[hydrofobic_feature,polar_feature] = hydrofobic_and_polar(variants,known);

%% Build table - Known data set
clc;
known_with_features = known;

codon_names = final_struct_avg.fields';
codon_count_table = array2table(avg,'variableNames',codon_names);
known_with_features = [known_with_features codon_count_table];

known_with_features.Afinal = Afinal;
known_with_features.Cfinal = Cfinal;
known_with_features.Gfinal = Gfinal;
known_with_features.Tfinal = Tfinal;

binar_first_table = array2table(binar_first,'variableNames', binar_first_names);
binar_second_table = array2table(binar_second,'variableNames', binar_second_names);
known_with_features = [known_with_features binar_first_table binar_second_table];

PSSM_table = array2table(dominance_bins,'variableNames', PSSM_names);
known_with_features = [known_with_features PSSM_table];

known_with_features.CpG_feature = CpG_feature;

known_with_features.hydrofobic_feature = hydrofobic_feature;
known_with_features.polar_feature = polar_feature;

known_features = known_with_features(:,28:end);
writetable(known_with_features, 'Known_with_features.xlsx','Sheet',1);
writetable(known_features, 'Known_features.xlsx','Sheet',1);


%% Features extraction - Unknown data set

% Codons Count
[final_struct_avg,avg] = codon_count_mat(sequences,test_rows,unknown);

% Count bases 
[Afinal,Cfinal,Gfinal,Tfinal] = basecounting(sequences,test_rows,unknown);

% Find bases dominants between two halfs of the sequence
[binar_first,binar_second] = Dominant(sequences,test_rows,unknown);
binar_first_names  = {'A Dominant 1', 'C Dominant 1' ,...
                      'G Dominant 1', 'T Dominant 1'};
binar_second_names = {'A Dominant 2', 'C Dominant 2' ,...
                      'G Dominant 2', 'T Dominant 2'};

% PSSM for average of 3 bases
% A dominance = 0 ; C dominance = 1 ; G dominance = 2 ; T dominance = 3
bins = unknown(:,1:27);
[A_rate, T_rate, C_rate, G_rate, seqs, dominance_seqs, dominance_bins] = PSSM(sequences,bins);
PSSM_names = { '1 - dominant' , '2 - dominant' ,...
               '3 - dominant' , '4 - dominant' ,...
               '5 - dominant' , '6 - dominant' ,...
               '7 - dominant' , '8 - dominant'};
% CpG
CpG_feature = CpG(variants,unknown);

% Hydrofobic & Polar codons
[hydrofobic_feature,polar_feature] = hydrofobic_and_polar(variants,unknown);

%% Build table - Unknown data set
clc;
unknown_with_features = unknown;

codon_names = final_struct_avg.fields';
codon_count_table = array2table(avg,'variableNames',codon_names);
unknown_with_features = [unknown_with_features codon_count_table];

unknown_with_features.Afinal = Afinal;
unknown_with_features.Cfinal = Cfinal;
unknown_with_features.Gfinal = Gfinal;
unknown_with_features.Tfinal = Tfinal;

binar_first_table = array2table(binar_first,'variableNames', binar_first_names);
binar_second_table = array2table(binar_second,'variableNames', binar_second_names);
unknown_with_features = [unknown_with_features binar_first_table binar_second_table];

PSSM_table = array2table(dominance_bins,'variableNames', PSSM_names);
unknown_with_features = [unknown_with_features PSSM_table];

unknown_with_features.CpG_feature = CpG_feature;

unknown_with_features.hydrofobic_feature = hydrofobic_feature;
unknown_with_features.polar_feature = polar_feature;

unknown_features=unknown_with_features(:,28:end);
writetable(unknown_with_features, 'unknown_with_features.xlsx','Sheet',1);
writetable(unknown_features, 'unknown_features.xlsx','Sheet',1);

end
