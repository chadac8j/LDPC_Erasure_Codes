% This simulation will use LDPC as Erasure Codes
% Chad Cole
% Sept 7th, 2021

close all

% PER = 0.01;  % Packet Error Rate
% PER = 0.05;  % Packet Error Rate
% PER = 0.1;  % Packet Error Rate
% PER = 0.15;  % Packet Error Rate
% PER = 0.2;  % Packet Error Rate
PER = 0.23;  % Packet Error Rate
% PER = 0.25;  % Packet Error Rate
% PER = 0.3;  % Packet Error Rate
% PER = 0.35;  % Packet Error Rate
% PER = 0.36;  % Packet Error Rate
% PER = 0.4;  % Packet Error Rate

m = 1;
GF_SIZE = 2^m;

% load H1000n_500k_3_6
% load H2223_rate075_3_12
load H_n603_k302_no6cycles

n = length(H_sparse(1,:));
n_k = length(H_sparse(:,1));
k = n - n_k;
rate = k/n;

H = full(H_sparse);
H = rearrange_cols(H);
%swap full rank cols to last part
temp = H(:, 1:n-k);
H(:, 1:n-k) = H(:, k+1:n);
H(:, k+1:n) = temp;
C2=H(:, k+1:n);
%only need to do this N^3 operation once
C2_inv = inv_GF2(C2);
H_system = mod(C2_inv*H, 2);
P_trans = H_system(:, 1:k);
G = cat(2, eye(k), P_trans');
H_sparse = sparse(H);

% Set up convenient structures for holding c-node/v-node graph connections
Vlist=zeros(n-k,max(sum(H_sparse,2))+1); 
Clist=zeros(n,max(sum(H_sparse,1))+1); 

for jj=1:n-k
    Vlist(jj,1)=sum(H_sparse(jj,:));
    icnt=0;
    for ii=1:n
        if H_sparse(jj,ii)==1
            icnt=icnt+1;
            Vlist(jj,icnt+1)=ii;
        end
    end
end

for ii=1:n
    Clist(ii,1)=sum(H_sparse(:,ii));
    jcnt=0;
    for jj=1:n-k
        if H_sparse(jj,ii)==1
            jcnt=jcnt+1;
            Clist(ii,jcnt+1)=jj;
        end
    end
end

num_trials = 1000;
rx_symbol_size_hist = zeros(1, n);
PER_vec = [0.2];
% PER_vec = [0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25];
% PER_vec = [0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5];
P_block = zeros(1, length(PER_vec));
P_block_ML = zeros(1, length(PER_vec));
time_tot_optimized = 0;
time_tot_not_optimized = 0;

for ebno_loop = 1:length(PER_vec)
    PER = PER_vec(ebno_loop)
    num_block_errors = 0;
    num_block_errors_ML = 0;
    iterations_tot = 0;
    iter = 0;
    while ((iter < num_trials) && (num_block_errors_ML < 100))
        iter = iter + 1;
        if mod(iter, 100000) == 0
            iter
        end
        source_vec = round(rand(1, k)); %create source symbols over GF(2)
        source_encode_vec = mod(source_vec*G,2);

%         source_vec = zeros(1, k); % Send all zero's codeword to speed up encoding
%         source_encode_vec = zeros(1, n); %Send all 0's vector

        num_erasures = 0;
        % Send n symbols
        recv_vec = source_encode_vec;
        for sym_ind = 1:n
            if (rand(1) <= PER) %Erasure represented by '-1'
                recv_vec(sym_ind) = -1;
                num_erasures = num_erasures + 1;
            end
        end

        % Test whether ML erasure decoder can decode (full rank received
        % vector set)
        recv_col_ind = find(recv_vec >= 0);
        recv_mat = G(:, recv_col_ind);
        rank_rx = gfrank(recv_mat, 2);
        if (rank_rx < k)
            num_block_errors_ML = num_block_errors_ML + 1
        end
        
        tStart1 = tic;
        [Out_LDPC, iterations] = My_ML_LDPC_Erasure_Decoder(recv_vec, G, k);
        tEnd1 = toc(tStart1);
        time_tot_optimized = time_tot_optimized + tEnd1;
        
        tStart2 = tic;
        [Out_LDPC2, iterations] = My_ML_LDPC_Erasure_Decoder_No_Remove_Zero_Rows(recv_vec, G, k);
        tEnd2 = toc(tStart2);
        time_tot_not_optimized = time_tot_not_optimized + tEnd2;
        
        if (sum(Out_LDPC == Out_LDPC2) ~= k)
            disp("Decoders not matching")
        end
        
%         [Out_LDPC, iterations] = My_LDPC_Erasure_Decoder(recv_vec, Vlist, Clist);    % Collect stats on how many symbols we receive
        iterations_tot = iterations_tot + iterations;
        if (num_erasures > 0)
            rx_symbol_size_hist(num_erasures) = rx_symbol_size_hist(num_erasures) + 1;
        end

        if sum(Out_LDPC2(1:k) == source_vec) ~= k
%         if sum(Out_LDPC == source_encode_vec) ~= n  %for non-systematic codes
            num_block_errors = num_block_errors + 1
        end
    end
    P_block(ebno_loop) = num_block_errors/iter;
    P_block_ML(ebno_loop) = num_block_errors_ML/iter;
    out_str = sprintf('The average number of decode iterations for PER %0.2f is: %0.2f ', PER, iterations_tot/num_trials);
    disp(out_str)
    
end

plot(PER_vec*100, P_block, 'r', PER_vec*100, P_block_ML, 'b')
title_str = sprintf("Block Error Rate vs. Packet Error Rate.  (%d, %d) LDPC Code", n, k);
title(title_str)
legend("Message Passing Decoder", "Maximum Likelihood Decoder")
xlabel("Packet Error Rate")
ylabel("Probability of Block Error")

