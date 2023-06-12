function [binar_first,binar_second] = Dominant(sequences,rows,file)
%%
%dominant in half of the sequence
binar_first = zeros(rows, 4);
binar_second = zeros(rows, 4);
%
A1 = zeros(length(sequences), 1);
C1 = zeros(length(sequences), 1);
G1 = zeros(length(sequences), 1);
T1 = zeros(length(sequences), 1);
%
A2 = zeros(length(sequences), 1);
C2 = zeros(length(sequences), 1);
G2 = zeros(length(sequences), 1);
T2 = zeros(length(sequences), 1);
for i=1:length(sequences)
     sequence = char(sequences(i));
     NTStruct1 = basecount(sequence(27:38));
     NTStruct2 = basecount(sequence(39:50));
     A1(i,1) = NTStruct1.A ;
     C1(i,1) = NTStruct1.C ;
     G1(i,1) = NTStruct1.G ;
     T1(i,1) = NTStruct1.T ;
     A2(i,1) = NTStruct2.A ;
     C2(i,1) = NTStruct2.C ;
     G2(i,1) = NTStruct2.G ;
     T2(i,1) = NTStruct2.T ;
end
first= horzcat(A1,C1,G1,T1);
second= horzcat(A2,C2,G2,T2);
for j = 1:rows
    vars = table2array(file(j,2:27)); %INDICES IN BIN
    vars = vars(~isnan(vars));
        Acurrent1 = zeros(1,length(vars(1))); %prectange of A in the bin
        Ccurrent1 = zeros(1,length(vars(1))); %prectange of C in the bin
        Gcurrent1 = zeros(1,length(vars(1)));%prectange of G in the bin
        Tcurrent1 = zeros(1,length(vars(1))); %prectange of T in the bin
        for y = 1:length(vars) %y is gonna be one index in the bin
        Acurrent(y) = first(vars(y),1);%adds the precentage of the relevant bin to the vector
        Ccurrent(y) = first(vars(y),2);
        Gcurrent(y) = first(vars(y),3);
        Tcurrent(y) = first(vars(y),4);
        end
        avg_vec = [mean(Acurrent),mean(Ccurrent),mean(Gcurrent),mean(Tcurrent)];
        maximum = max(avg_vec);
        maximal = find(avg_vec == maximum);
        binar_first(j,:) = zeros(1,4);
        binar_first(j,maximal) = 1;
end     
for j = 1:rows
    vars = table2array(file(j,2:27)); %INDICES IN BIN
    vars = vars(~isnan(vars));
        Acurrent1 = zeros(1,length(vars(1)));
        Ccurrent1 = zeros(1,length(vars(1)));
        Gcurrent1 = zeros(1,length(vars(1)));
        Tcurrent1 = zeros(1,length(vars(1)));
        for y = 1:length(vars) %y is gonna be one index in the bin
        Acurrent(y) = second(vars(y),1);
        Ccurrent(y) = second(vars(y),2);
        Gcurrent(y) = second(vars(y),3);
        Tcurrent(y) = second(vars(y),4);
        end
        avg_vec = [mean(Acurrent),mean(Ccurrent),mean(Gcurrent),mean(Tcurrent)];
        maximum = max(avg_vec);
        maximal = find(avg_vec == maximum);
        binar_second(j,:) = zeros(1,4);
        binar_second(j,maximal) = 1;
end
end

