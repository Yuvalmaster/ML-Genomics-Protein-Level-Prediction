function [A_rate, T_rate, C_rate, G_rate, seqs, dominance_seqs, dominance_bins] = PSSM(Seqs,bins)

    % PSSM for each cell
    for i = 1:length(Seqs)                  
      seqs{i} = Seqs{i}(27:50);       
      for k = 1:24
          bases = basecount(seqs{i}(k));% Count bases
          A(k,i) = bases.A;
          C(k,i) = bases.C;
          T(k,i) = bases.T;
          G(k,i) = bases.G;
      end
    end
   
    A_rate = []; C_rate = []; G_rate = []; T_rate = [];
    for i = 1:24
        A_rate = [A_rate ; sum(A(i,:))/length(Seqs)];
        C_rate = [C_rate ; sum(C(i,:))/length(Seqs)];
        G_rate = [G_rate ; sum(G(i,:))/length(Seqs)];
        T_rate = [T_rate ; sum(T(i,:))/length(Seqs)];
    end

   %% Dominance in sequences
    % A dominance = 0 ; C dominance = 1 ; G dominance = 2 ; T dominance = 3
    dominance_seqs = nan(length(Seqs),8);
    for j = 1:length(Seqs)
        % first three slots
        if sum(A(1:3,j)/3) > mean(A_rate(1:3))
            dominance_seqs(j,1) = 0;
        elseif sum(C(1:3,j)/3) > mean(C_rate(1:3))
           dominance_seqs(j,1) = 1;
        elseif sum(G(1:3,j)/3) > mean(G_rate(1:3))
           dominance_seqs(j,1) = 2;
        elseif sum(T(1:3,j)/3) > mean(T_rate(1:3))
           dominance_seqs(j,1) = 3;
        end
        % second three slots
        if sum(A(4:6,j)/3) > mean(A_rate(4:6))
            dominance_seqs(j,2) = 0;
        elseif sum(C(4:6,j)/3) > mean(C_rate(4:6))
           dominance_seqs(j,2) = 1;
        elseif sum(G(4:6,j)/3) > mean(G_rate(4:6))
           dominance_seqs(j,2) = 2;
        elseif sum(T(4:6,j)/3) > mean(T_rate(4:6))
           dominance_seqs(j,2) = 3;
        end
        % third three slots
        if sum(A(7:9,j)/3) > mean(A_rate(7:9))
            dominance_seqs(j,3) = 0;
        elseif sum(C(7:9,j)/3) > mean(C_rate(7:9))
           dominance_seqs(j,3) = 1;
        elseif sum(G(7:9,j)/3) > mean(G_rate(7:9))
           dominance_seqs(j,3) = 2;
        elseif sum(T(7:9,j)/3) > mean(T_rate(7:9))
           dominance_seqs(j,3) = 3;
        end
        % Fourth three slots
        if sum(A(10:12,j)/3) > mean(A_rate(10:12))
            dominance_seqs(j,4) = 0;
        elseif sum(C(10:12,j)/3) > mean(C_rate(10:12))
           dominance_seqs(j,4) = 1;
        elseif sum(G(10:12,j)/3) > mean(G_rate(10:12))
           dominance_seqs(j,4) = 2;
        elseif sum(T(10:12,j)/3) > mean(T_rate(10:12))
           dominance_seqs(j,4) = 3;
        end
        % fifth three slots
        if sum(A(13:15,j)/3) > mean(A_rate(13:15))
           dominance_seqs(j,5) = 0;
        elseif sum(C(13:15,j)/3) > mean(C_rate(13:15))
           dominance_seqs(j,5) = 1;
        elseif sum(G(13:15,j)/3) > mean(G_rate(13:15))
           dominance_seqs(j,5) = 2;
        elseif sum(T(13:15,j)/3) > mean(T_rate(13:15))
           dominance_seqs(j,5) = 3;
        end
        % sixth three slots
        if sum(A(16:18,j)/3) > mean(A_rate(16:18))
           dominance_seqs(j,6) = 0;
        elseif sum(C(16:18,j)/3) > mean(C_rate(16:18))
           dominance_seqs(j,6) = 1;
        elseif sum(G(16:18,j)/3) > mean(G_rate(16:18))
           dominance_seqs(j,6) = 2;
        elseif sum(T(16:18,j)/3) > mean(T_rate(16:18))
           dominance_seqs(j,6) = 3;
        end
        % seventh three slots
        if sum(A(19:21,j)/3) > mean(A_rate(19:21))
           dominance_seqs(j,7) = 0;
        elseif sum(C(19:21,j)/3) > mean(C_rate(19:21))
           dominance_seqs(j,7) = 1;
        elseif sum(G(19:21,j)/3) > mean(G_rate(19:21))
           dominance_seqs(j,7) = 2;
        elseif sum(T(19:21,j)/3) > mean(T_rate(19:21))
           dominance_seqs(j,7) = 3;
        end
        % final three slots
        if sum(A(22:24,j)/3) > mean(A_rate(22:24))
            dominance_seqs(j,8) = 0;
        elseif sum(C(22:24,j)/3) > mean(C_rate(22:24))
           dominance_seqs(j,8) = 1;
        elseif sum(G(22:24,j)/3) > mean(G_rate(22:24))
           dominance_seqs(j,8) = 2;
        elseif sum(T(22:24,j)/3) > mean(T_rate(22:24))
           dominance_seqs(j,8) = 3;
        end

    end
    %% Dominance in bins
    dominance_bins = zeros(size(bins,1),8);
    for i = 1:size(bins,1)
        for j = 1:8
            bin = table2array(bins(i,2:end));
            bin = bin(~isnan(bin));
            count = dominance_seqs(bin,j);
            sums = [sum(count == 0) sum(count == 1) sum(count == 2)...
                    sum(count == 3)];
            [~,max_val] = max(sums);
            if max_val == 1
                dominance_bins(i,j) = 0;
            elseif max_val == 2
                dominance_bins(i,j) = 1;
            elseif max_val == 3
                dominance_bins(i,j) = 2;
            elseif max_val == 4
                dominance_bins(i,j) = 3;
            end
        end
    end

end