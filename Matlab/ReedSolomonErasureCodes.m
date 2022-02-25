% Reed Solomon Erasure Codes
% Chad Cole
% Aug 3rd, 2021

close all

% PER = 0.01;  % Packet Error Rate
PER = 0.05;  % Packet Error Rate
% PER = 0.1;  % Packet Error Rate
% PER = 0.25;  % Packet Error Rate

m = 8;
n = 2^m-1;
% n = 2^m-100;
GF_SIZE = 2^m;
% num_symbols_to_correct = 20;
% num_symbols_to_correct = 40;
% num_symbols_to_correct = 50;
% num_symbols_to_correct = 60;
num_symbols_to_correct = 63;  % This n=255, k=192 is the TIA-5041 standard RS code
k = n - num_symbols_to_correct;

prim_poly_m8 = [1 0 1 1 1 0 0 0 1];

[GF_add_lookup, GF_mult_lookup, GF_inv_lookup, G] = Build_GF256_Lookup_Tables(m, k, n, prim_poly_m8);

% field = gftuple([-1:2^m-2]',prim_poly_m8,2);
% Find G matrix; Use Vandermonde matrix as specified in rfc5510
% G = zeros(k, n);
% Eventually want a systematic version of G
G_k_inv = inv(G(1:k, 1:k));
G_sys = G_k_inv*G;

num_block_errors = 0;
num_symbol_errors = 0;
num_trials = 10000;
rx_symbol_size_hist = zeros(1, n);
rx_symbol_size_errors_hist = zeros(1, n);

G_sys_real = zeros(size(G_sys));
for ll = 1:prod(size(G_sys))
    G_sys_real(ll) = double(G_sys.x(ll));
end
    
num_sys_symbols_hist = zeros(1, num_trials);
dec_time_hist = zeros(1, num_trials);

for iter = 1:num_trials
    if mod(iter, 1000) == 0
        iter
    end
    source_vec = round(255*rand(1, k)); %create source symbols over GF(2^8)
    source_encode_vec = source_vec*G_sys;
    source_encode_vec_real = zeros(1, n);
    for ll = 1:(length(source_encode_vec))
        source_encode_vec_real(ll) = double(source_encode_vec.x(ll));
    end
    source_encode_vec = source_encode_vec_real;
    
%     source_vec = zeros(1, k); % Send all zero's codeword to speed up encoding
%     source_encode_vec =  gf(zeros(1, n), m, prim_poly_m8_number);
%     source_encode_vec =  zeros(1, n);
    
    recv_vec_ind = zeros(1, n);
    num_sym_rx = 0;
    num_systematic_rx = 0;
    % Send n symbols
    for sym_ind = 1:n
        if (rand(1) > PER)
       % add to rx indices list
            num_sym_rx = num_sym_rx + 1;
            recv_vec_ind(num_sym_rx) = sym_ind;
            if (sym_ind <= k)
                num_systematic_rx = num_systematic_rx + 1;
            end
        end
    end
    Out_RS = zeros(1, k);
    
    if num_sym_rx >= k
        rx_vec = source_encode_vec(recv_vec_ind(1:k)); 
%         tic
        tStart = tic; 
%         Out_RS = My_RS_Decode(recv_vec_ind(1:k), rx_vec, m, n, k, prim_poly_m8_number, G_sys, log_lookup);
        Out_RS = My_RS_Decode_Optimize_With_GFTables(recv_vec_ind(1:k), rx_vec, m, n, k, G_sys_real, GF_add_lookup, GF_mult_lookup, GF_inv_lookup);
%         toc
        tEnd = toc(tStart);
%         disp(["Time: ", tEnd, "Number of Systematic symbols: ", num_systematic_rx]);
        num_sys_symbols_hist(iter) = num_systematic_rx;
        dec_time_hist(iter) = tEnd;
    end
    % Collect stats on how many symbols we receive
    rx_symbol_size_hist(num_sym_rx) = rx_symbol_size_hist(num_sym_rx) + 1;
    % Does rx_mat have at least k linearly independent rows?
    if sum(Out_RS == source_vec) ~= k
        num_block_errors = num_block_errors + 1
    end
end

% percent_err = (rx_symbol_size_errors_hist./(eps + rx_symbol_size_hist))*100;
% 
% figure
% subplot(2, 1, 1)
% b = bar(percent_err);
% xtips1 = b(1).XEndPoints;
% ytips1 = b(1).YEndPoints;
% labels1 = string(rx_symbol_size_hist);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% title("Percent error when received X symbols.  k = 10.  Random Code")
% ylim([0, 110]);
% subplot(2, 1, 2)
% 

stem(num_sys_symbols_hist, dec_time_hist)
title_str = [string(PER*100), "% Packet Error Rate.  (156, 96) RS Code"];
title(title_str(1) + title_str(2))
xlabel("Source Symbols Received")
ylabel("Decode Time (seconds)")

