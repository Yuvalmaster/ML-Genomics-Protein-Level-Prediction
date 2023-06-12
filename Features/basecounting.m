function [Afinal,Cfinal,Gfinal,Tfinal] = basecounting(sequences,rows, file)

%loop of summing the bases in each sequence
A = zeros(length(sequences), 1);
C = zeros(length(sequences), 1);
G = zeros(length(sequences), 1);
T = zeros(length(sequences), 1);
for i=1:length(sequences)
    sequence = char(sequences(i));
    NTStruct = basecount(sequence(27:50));
    A(i,1) = NTStruct.A /length(sequence(27:50));
    C(i,1) = NTStruct.C /length(sequence(27:50));
    G(i,1) = NTStruct.G /length(sequence(27:50));
    T(i,1) = NTStruct.T /length(sequence(27:50));
end
%%
%each of the following vectors is a cloum vector, his rows are (%) percentage
%of this base in the bin. calculated by averaging the % of this base for all sequences in the bin 
Afinal = zeros(rows, 1);
Cfinal = zeros(rows, 1);
Gfinal = zeros(rows, 1);
Tfinal = zeros(rows, 1);
for j = 1:rows
    vars = table2array(file(j,2:27)); %INDICES IN BIN
        Acurrent = zeros(1,length(vars(1)));
        Ccurrent = zeros(1,length(vars(1)));
        Gcurrent = zeros(1,length(vars(1)));
        Tcurrent = zeros(1,length(vars(1)));
        for y = length(vars(1)) %y is gonna be one index in the bin
        Acurrent(y) = A(vars(y),1);
        Ccurrent(y) = C(vars(y),1);
        Gcurrent(y) = G(vars(y),1);
        Tcurrent(y) = T(vars(y),1);
        end
        Afinal(j) = mean(Acurrent);
        Cfinal(j) = mean(Ccurrent);
        Gfinal(j) = mean(Gcurrent);
        Tfinal(j) = mean(Tcurrent);
end
end

