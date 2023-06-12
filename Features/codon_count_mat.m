function [final_struct_avg,avg] = codon_count_mat(sequences,rows,file)
%
%
%sequences should be the variants
%rows and cols should be the size of the data
%file should be known/ unknown
%loop of summing the codon count accroding to 3 reading frames, inserting
%the final result to the codon count mat
%
%
%initiates the codon count mat in order to call it later (Global)
codon_count = zeros(length(sequences), 64); %rows are variants, colums are codons
%
for i=1:length(sequences)
    sequence = char(sequences(i));
    count_frame_1 = codoncount(sequence(27:50),'frame',1);
    count_frame_2 = codoncount(sequence(27:50),'frame',2);
    count_frame_3 = codoncount(sequence(27:50),'frame',3);
    h = horzcat(count_frame_1,count_frame_2,count_frame_3); %cols are codons, rows are frames of reading 1,2,3
    addition = cell2mat(struct2cell(h(1)))+cell2mat(struct2cell(h(2)))+cell2mat(struct2cell(h(3))); %summing the reading frames
    codon_count(i,:) = addition'; %inserting to the mat
end
%
%for calculating the codon count for each bin, we are taking codon count
%for each variant and avaraging the codon count for each codon on the
%variants belong to the same bin
avg = zeros(rows,64); % codon are cols, bins are rows
for j = 1:rows
    vars = table2array(file(j,2:27)); %INDICES IN BIN
    vars = vars(~isnan(vars));
    relevant_codons = zeros(length(vars),64);
    for x=1:length(vars)
        v = vars(x);
        relevant_codons(x,:) = codon_count(v,:); %picking the relenant roew from codon count
    end
    avg(j,:) = mean(relevant_codons); 
end
%creating a struct in order to have an easier view on the codon count.
sequence = char(sequences(1));
count_frame_1 = codoncount(sequence(27:50),'frame',1);
fields = fieldnames(count_frame_1);
final_struct_avg = table(fields, avg');
end

