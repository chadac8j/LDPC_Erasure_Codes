% This simulation will use Non-binary LDPC codes as Erasure Codes
% Chad Cole
% Oct 4th, 2021

% close all

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

m = 8;
GF_SIZE = 2^m;

% load H1000n_500k_3_6
% load n1000_k500_H_no6cycles
% load n700_k400_no6cycles H_sparse
% load n1200_k600_triangle_H H_sparse
% load n400_k200_no6cycles_triangleForm H_sparse
% load n400_k200_no6cycles_triangleForm_GF256_2ndAttempt H_sparse_nb H_sparse
% load n600_k450_no6cycles_triangleForm H_sparse
% load n600_k450_no6cycles_triangleForm_GF256 H_sparse H_sparse_nb
% load n1020_k816_no6cycles_triangleForm H_sparse
% load n1020_k816_no6cycles_triangleForm_GF256 H_sparse H_sparse_nb
load n2040_k1632_no6cycles_triangleForm H_sparse H_sparse_nb

% load n400_k200_no6cycles_triangleForm_GF256 H_sparse_nb H_sparse
% load n400_k200_no6cycles_triangleForm_2 H_sparse
% load n400_k200_no6cycles_triangleForm_3 H_sparse
% load n960_k720_triangleFormLDPC H_sparse
% load n180_k120_triangleFormLDPC H_sparse

% load n224_k196_14_square_grid H_sparse
% load n65_k50_10by5_grid H_sparse
% load H2223_rate075_3_12
% load H_n603_k302_no6cycles


n = length(H_sparse(1,:));
n_k = length(H_sparse(:,1));
k = n - n_k;
rate = k/n;

% % % % Turn binary H_sparse into GF(Q)
% H_sparse_nb = H_sparse;
% for ii = 1:n-k
%     find_ind = find(H_sparse(ii,:));
%     for jj = 1:length(find_ind)
%         H_sparse_nb(ii,find_ind(jj)) = floor((GF_SIZE-1)*rand(1))+1;
%     end
% end
% 
% % Verify that new matrix has same structure
% for ii = 1:n-k
%     find_ind = find(H_sparse(ii,:));
%     find_ind_nb = find(H_sparse_nb(ii,:));
%     non_zero_len = max(length(find_ind), length(find_ind_nb));
%     if sum(find_ind == find_ind_nb) ~= non_zero_len
%        disp("New non-binary H is not same form as old one") 
%     end
% end

prim_poly_m8 = [1 0 1 1 1 0 0 0 1];

% [GF_add_lookup, GF_mult_lookup, GF_inv_lookup, G] = Build_GF256_Lookup_Tables(m, k, n, prim_poly_m8);
load GF_256_add_mult_inv_tables GF_add_lookup GF_mult_lookup GF_inv_lookup


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
Vlist_val=zeros(n-k,max(sum(H_sparse,2))+1); 
Clist_val=zeros(n,max(sum(H_sparse,1))+1); 

for jj=1:n-k
    Vlist(jj,1)=sum(H_sparse(jj,:));
    Vlist_val(jj,1) = Vlist(jj,1);
    icnt=0;
    for ii=1:n
        if H_sparse(jj,ii)==1
            icnt=icnt+1;
            Vlist(jj,icnt+1)=ii;
            Vlist_val(jj, icnt+1)=H_sparse_nb(jj, ii);
        end
    end
end

for ii=1:n
    Clist(ii,1)=sum(H_sparse(:,ii));
    Clist_val(ii,1)=Clist(ii,1);
    jcnt=0;
    for jj=1:n-k
        if H_sparse(jj,ii)==1
            jcnt=jcnt+1;
            Clist(ii,jcnt+1)=jj;
            Clist_val(ii,jcnt+1)=H_sparse_nb(jj,ii);
        end
    end
end

num_trials = 100000;
rx_symbol_size_hist = zeros(1, n);
% PER_vec = [0.4];
% PER_vec = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14];
% PER_vec = [0.05 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25];
% PER_vec = [0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5];
% PER_vec = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5];


alpha_vec = [0.02, 0.04, 0.06, 0.08, 0.1];
beta = 0.4;
transition = 0.1;
good_transition_bias = 10;
PER_vec = zeros(1, length(alpha_vec));
for ii = 1:length(PER_vec)
     PER_vec(ii) = (1/(1 + 1/good_transition_bias))*alpha_vec(ii) + (1 - 1/(1 + 1/good_transition_bias))*beta;
end
PER_vec_sim = zeros(1, length(alpha_vec));


P_block = zeros(1, length(PER_vec));
P_block_ML = zeros(1, length(PER_vec));
P_block_RS = zeros(1, length(PER_vec));
time_tot_iterative = 0;
time_tot_iterative_ML = 0;
Out_LDPC_Iterative = zeros(1, n);
Out_LDPC_Iterative_ML = zeros(1, n);

% For RS comparison purposes, we need an (n,k) for the RS code
n_RS = 255;
k_RS = round(rate*n_RS);

for ebno_loop = 1:length(PER_vec)
    PER = PER_vec(ebno_loop)
    alpha = alpha_vec(ebno_loop)  % param for bursty error channels
    num_block_errors_Iterative = 0;
    num_block_errors_RS = 0;
    num_block_errors_ML = 0;
    num_erasures_tot = 0;
    iterations_tot = 0;
    iter = 0;
    next_state = 0;
    while ((iter < num_trials) && (num_block_errors_Iterative < 100))
        iter = iter + 1;
        if mod(iter, 100000) == 0
            iter
        end
%         source_vec = round(rand(1, k)); %create source symbols over GF(2)
%         source_encode_vec = mod(source_vec*G,2);

% Encoding of Non-binary systematic triangle or IRA codes
        source_vec = floor(GF_SIZE*rand(1, k)); %create source symbols over GF(2)
        source_encode_vec = zeros(1, n);
        source_encode_vec(1:k) = source_vec; %Send all 0's vector
        for pp = 1:(n-k)
            gf_sum = 0;
            for ll = 1:Vlist_val(pp, 1)-1 % Don't use last nonzero value, which is on the diagonal
                gf_sum = GF_add_lookup(gf_sum + 1, GF_mult_lookup(source_encode_vec(Vlist(pp, ll+1)) + 1, Vlist_val(pp, ll+1) + 1) + 1);
            end
            source_encode_vec(k+pp) = GF_mult_lookup(gf_sum + 1, GF_inv_lookup(Vlist_val(pp, Vlist_val(pp, 1)+1)) + 1);
        end


%         source_vec = zeros(1, k); % Send all zero's codeword to speed up encoding
%         source_encode_vec = zeros(1, n); %Send all 0's vector

        num_erasures = 0;
        % Send n symbols
        recv_vec = source_encode_vec;
        for sym_ind = 1:n
            [cur_error, next_state] = Bursty_Error_Channel_Model_Generator(next_state, alpha, beta, good_transition_bias);
            if (cur_error == 1)
%             if (rand(1) <= PER) %Erasure represented by '-1'
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
        for ll = 1:ceil(n/n_RS)
            if (sum(source_encode_vec((ll-1)*n_RS+1:ll*n_RS) == recv_vec((ll-1)*n_RS+1:ll*n_RS)) < k_RS)
                num_block_errors_RS = num_block_errors_RS + 1
            end
        end

        if (num_erasures < (n-k))
            tStart1 = tic;
            [Out_LDPC_Iterative, iterations] = My_LDPC_HybridML_NonBinary_Erasure_Decoder(recv_vec, Vlist, Clist, H_sparse_nb, n, k, GF_add_lookup, GF_mult_lookup, GF_inv_lookup);
            tEnd1 = toc(tStart1);
            time_tot_iterative = time_tot_iterative + tEnd1;
        end
        
%         [Out_LDPC, iterations] = My_LDPC_Erasure_Decoder(recv_vec, Vlist, Clist);    % Collect stats on how many symbols we receive
        iterations_tot = iterations_tot + iterations;

        if (num_erasures > 0)
            rx_symbol_size_hist(num_erasures) = rx_symbol_size_hist(num_erasures) + 1;
        end
        
%         % For Reed Solomon comparison, assume a rate 1/2 (256, 128) code        

        if sum(Out_LDPC_Iterative == source_encode_vec) ~= n  %for non-systematic codes
            num_block_errors_Iterative = num_block_errors_Iterative + 1
        end
                
        num_erasures_tot = num_erasures_tot + num_erasures;
    end
    P_block(ebno_loop) = num_block_errors_Iterative/iter;
    P_block_RS(ebno_loop) = num_block_errors_RS/(ceil(n/n_RS)*iter);
    out_str = sprintf('The average number of decode iterations for PER %0.2f is: %0.2f ', PER, iterations_tot/iter);
    disp(out_str)
    PER_vec_sim(ebno_loop) = num_erasures_tot/(iter*n);
end

% plot(PER_vec*100, P_block, 'r', PER_vec*100, P_block_ML, 'b')
semilogy(PER_vec*100, P_block, 'k', PER_vec*100, P_block_RS, 'b')
% semilogy(PER_vec*100, P_block, 'r', PER_vec*100, P_block_ML, 'b', PER_vec*100, P_block_RS, 'k')
% semilogy(PER_vec*100, P_block, 'b')
title_str = sprintf("Block Error Rate vs. Packet Error Rate.  (%d, %d) LDPC Code", n, k);
title(title_str)
legend("MP/Maximum Likelihood Decoder", "(255, k_RS) Reed Solomon Code")
% legend("Message Passing Decoder", "MP/Maximum Likelihood Decoder")
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
