% This simulation will use LDPC as Erasure Codes
% Chad Cole
% Sept 20th, 2021

% close all

% PER = 0.01;  % Packet Error Rate
% PER = 0.05;  % Packet Error Rate
% PER = 0.1;  % Packet Error Rate
% PER = 0.15;  % Packet Error Rate
% PER = 0.2;  % Packet Error Rate
% PER = 0.23;  % Packet Error Rate
% PER = 0.25;  % Packet Error Rate
% PER = 0.3;  % Packet Error Rate
% PER = 0.35;  % Packet Error Rate
% PER = 0.36;  % Packet Error Rate
% PER = 0.4;  % Packet Error Rate

m = 1;
GF_SIZE = 2^m;

% load H1000n_500k_3_6
% load n1000_k500_H_no6cycles
% load n700_k400_no6cycles H_sparse
% load n1200_k600_triangle_H H_sparse
% load n400_k200_no6cycles_triangleForm H_sparse
% load n400_k200_no6cycles_triangleForm_2 H_sparse
% load n400_k200_no6cycles_triangleForm_3 H_sparse
% load n960_k720_triangleFormLDPC H_sparse
% load n180_k120_triangleFormLDPC H_sparse
% load n600_k450_no6cycles_triangleForm_GF256 H_sparse
% load n1020_k816_no6cycles_triangleForm H_sparse
% load n2040_k1632_no6cycles_triangleForm H_sparse
% load n1020_k765_no6cycles_triangleForm H_sparse
% load n2040_k1530_no6cycles_triangleForm H_sparse
% load n2040_k1530_irreg_H_no6cycles_triangleForm H_sparse  % used for OpenCL
% load n2040_k1530_no6cycles_triangleForm_2 H_sparse
% load n4000_k2000_no6cycles_triangleForm H_sparse
% load n4000_k2000_no4cycles_triangleForm_2 H_sparse
% load n2000_k1000_no4cycles_triangleForm H_sparse
% load n2000_k1000_no6cycles_triangleForm H_sparse
% load n4080_k3060_irreg_H_no6cycles_triangleForm H_sparse

load n2000_k1000_no6cycles_triangleForm_OpenCL_H H_sparse

% load n224_k196_14_square_grid H_sparse
% load n65_k50_10by5_grid H_sparse
% load H2223_rate075_3_12
% load H_n603_k302_no6cycles

% read_from_file = 1;  % Used for debugging test vectors used in OpenCL
read_from_file = 0;  % Generate random test data

n = length(H_sparse(1,:));
n_k = length(H_sparse(:,1));
k = n - n_k;
rate = k/n;

% H = full(H_sparse);
% H = rearrange_cols(H);
% %swap full rank cols to last part
% temp = H(:, 1:n-k);
% H(:, 1:n-k) = H(:, k+1:n);
% H(:, k+1:n) = temp;
% C2=H(:, k+1:n);
% %only need to do this N^3 operation once
% C2_inv = inv_GF2(C2);
% H_system = mod(C2_inv*H, 2);
% P_trans = H_system(:, 1:k);
% G = cat(2, eye(k), P_trans');
% H_sparse = sparse(H);

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

% num_trials = 1000000;
% num_trials = 200000;
num_trials = 1000000;
rx_symbol_size_hist = zeros(1, n);
iterations_hist = zeros(1, 50);
% PER_vec = [0.1];
% PER_vec = [0.125];
% PER_vec = [0.15];
% PER_vec = [0.1875];
% PER_vec = [0.1719];
% PER_vec = [0.1562];
PER_vec = [0.1406];
% PER_vec = [0.2031];
% PER_vec = [0.4062]; % 26/64
% PER_vec = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14];

% PER_vec = [0.05 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25];
% PER_vec = [0.14, 0.16, 0.18, 0.2, 0.22];

% PER_vec = [0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5];
% PER_vec = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
P_block = zeros(1, length(PER_vec));
P_block_ML = zeros(1, length(PER_vec));
P_block_RS = zeros(1, length(PER_vec));
time_tot_iterative = 0;
time_tot_iterative_ML = 0;
Out_LDPC_Iterative = zeros(1, n);
Out_LDPC_Iterative_ML = zeros(1, n);

% % For RS comparison purposes, we need an (n,k) for the RS code
n_RS = 255;
k_RS = ceil(rate*n_RS);
% n_RS = 250;
% k_RS = round(rate*n_RS);

for ebno_loop = 1:length(PER_vec)
    PER = PER_vec(ebno_loop)
    num_block_errors_Iterative = 0;
    num_block_errors_RS = 0;
    num_block_errors_ML = 0;
    iterations_tot = 0;
    iter = 0;
    while ((iter < num_trials) && (num_block_errors_Iterative < 1000))
        iter = iter + 1;
        if mod(iter, 100000) == 0
            iter
        end
%         source_vec = round(rand(1, k)); %create source symbols over GF(2)
%         source_encode_vec = mod(source_vec*G,2);

% Encoding of systematic triangle or IRA codes
        if (read_from_file == 1)
            source_vec = zeros(1, k);
            file_str = sprintf("LDPC_ErasureEncoder_IN_k%d_Shorts.txt", k);
%             [fid_H,msg] = fopen(file_str,'r');
%             assert(fid_H>=3,msg);
            S = readlines(file_str);
            for ii=1:k %Assume S is at least k rows long
%                 source_vec(ii) = str2double(S(ii));
                source_vec(ii) = mod(str2double(S(ii)), 2); % keep in binary for original decoder
            end
%             fclose(fid_H);
        else
            source_vec = round(rand(1, k)); %create source symbols over GF(2)            
        end
        source_encode_vec = zeros(1, n);
        source_encode_vec(1:k) = source_vec; %Send all 0's vector
        for pp = 1:(n-k)
            source_encode_vec(k+pp) = mod(H_sparse(pp,1:k+pp-1)*source_encode_vec(1:k+pp-1)', 2);
% % The below block is for source symbols > 0 bit wide
%             cum_xor = 0;
%             for qq = 2:Vlist(pp, 1) % Don't need to XOR last col of source_encode_vec because we know it is '0'
%                 cum_xor = bitxor(cum_xor, source_encode_vec(Vlist(pp, qq)));
%             end
%             source_encode_vec(k + pp) = cum_xor;
        end


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

%         % Test whether ML erasure decoder can decode (full rank received
%         % vector set)
%         recv_col_ind = find(recv_vec >= 0);
%         recv_mat = G(:, recv_col_ind);
%         rank_rx = gfrank(recv_mat, 2);
%         if (rank_rx < k)
%             num_block_errors_ML = num_block_errors_ML + 1
%         end
        
        % Assume (n_LDPC/n_RS) equivalent rate RS codes of size n_RS are concatenated
        for ll = 1:floor(n/n_RS)
%         for ll = 1:ceil(n/n_RS)
            if (sum(source_encode_vec((ll-1)*n_RS+1:min(ll*n_RS, n)) == recv_vec((ll-1)*n_RS+1:min(ll*n_RS, n))) < k_RS)
                num_block_errors_RS = num_block_errors_RS + 1
            end
        end

        if (num_erasures > (n-k))
%             num_block_errors_RS = num_block_errors_RS + 1
        else
            tStart1 = tic;
            [Out_LDPC_Iterative, iterations] = My_LDPC_Erasure_Decoder(recv_vec, Vlist, Clist);    % Collect stats on how many symbols we receive
            tEnd1 = toc(tStart1);
            time_tot_iterative = time_tot_iterative + tEnd1;
            iterations_hist(iterations) = iterations_hist(iterations) + 1;

            tStart2 = tic;
            [Out_LDPC_Iterative_ML, iterations] = My_LDPC_HybridML_Erasure_Decoder(recv_vec, Vlist, Clist, H_sparse, n);
            tEnd2 = toc(tStart2);
            time_tot_iterative_ML = time_tot_iterative_ML + tEnd2;
        end
        
%         [Out_LDPC, iterations] = My_LDPC_Erasure_Decoder(recv_vec, Vlist, Clist);    % Collect stats on how many symbols we receive
        iterations_tot = iterations_tot + iterations;

        if (num_erasures > 0)
            rx_symbol_size_hist(num_erasures) = rx_symbol_size_hist(num_erasures) + 1;
        end
        
        if sum(Out_LDPC_Iterative == source_encode_vec) ~= n  %for non-systematic codes
            num_block_errors_Iterative = num_block_errors_Iterative + 1
            sum(Out_LDPC_Iterative ~= source_encode_vec)  % See if this is a small stopping set
        end
        
        if sum(Out_LDPC_Iterative_ML == source_encode_vec) ~= n  %for non-systematic codes
            num_block_errors_ML = num_block_errors_ML + 1
        end
        
    end
    P_block(ebno_loop) = num_block_errors_Iterative/iter;
    P_block_RS(ebno_loop) = num_block_errors_RS/(ceil(n/n_RS)*iter);
    P_block_ML(ebno_loop) = num_block_errors_ML/iter;
    out_str = sprintf('The average number of decode iterations for PER %0.2f is: %0.2f ', PER, iterations_tot/iter);
    disp(out_str)
    
end

% plot(PER_vec*100, P_block, 'r', PER_vec*100, P_block_ML, 'b')
% semilogy(PER_vec*100, P_block, 'r', PER_vec*100, P_block_ML, 'b')
semilogy(PER_vec*100, P_block, 'r:', PER_vec*100, P_block_ML, 'b-.', PER_vec*100, P_block_RS, 'k-','LineWidth',3)
% semilogy(PER_vec*100, P_block_ML, 'c')
title_str = sprintf("Block Error Rate vs. Packet Error Rate.  (%d, %d) LDPC Code", n, k);
title(title_str)
% legend("Message Passing Decoder", "MP/Maximum Likelihood Decoder")
legend("LDPC Message Passing Decoder", "LDPC MP/Maximum Likelihood Decoder", "RS Code")
xlabel("Packet Error Rate")
ylabel("Probability of Block Error")

% % save Perf_n400_k200 P_block PER_vec
% figure
% semilogy(PER_vec*100, P_block, 'r', PER_vec*100, P_block_RS, 'b')
% title_str = sprintf("Rate 1/2 (%d, %d) LDPC Code vs. (256, 128) Reed Solomon Code", n, k);
% title(title_str)
% legend("LDPC (Message Passing Decoding)", "RS (Max Likelihood Decoding)")
% xlabel("Packet Error Rate")
% ylabel("Probability of Block Error")
