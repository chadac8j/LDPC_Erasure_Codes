% Chad Cole
% Aug 2nd, 2021

close all

% PER = 0.01;  % Packet Error Rate
PER = 0.05;  % Packet Error Rate

k = 10;
%n = 12;
n = 15;

gen_new_G = 0;

if (gen_new_G == 1)
% Try random codes first
    G = [eye(k), round(rand(k, n-k))];
% G_parity = [eye(k), ones(k, 1), round(rand(k, n-k-1))];  % parity + random code.  The parity column guarantees can correct one source symbol
    G_parity = [eye(k), ones(k, 1), (zeros(k, n-k-1))];
    for jj = k+2:n  %This block will give k/2 ones in each G column
    % Assume n-k-1 is even
        temp_perm = randperm(k);
        for ii = 1:(k)/2
            G_parity(temp_perm(ii), jj)=1;
        end
    end
end

figure
subplot(2, 1, 1)
spy(G)
subplot(2, 1, 2)
spy(G_parity)

num_block_errors = 0;
num_block_errors_2 = 0;

num_symbol_errors = 0;
num_symbol_errors_when_rx_at_least_k_symbols = 0;
num_symbol_errors_2 = 0;
num_symbol_errors_when_rx_at_least_k_symbols_2 = 0;
num_trials = 100000;
rx_symbol_size_hist = zeros(1, n);
rx_symbol_size_errors_hist = zeros(1, n);
rx_symbol_size_errors_hist_2 = zeros(1, n);

for iter = 1:num_trials
    rx_sym_vec = zeros(1, n);
    num_sym_rx = 0;
    % Send n symbols
    for sym_ind = 1:n
        if (rand(1) > PER)
       %received symbol ok
            num_sym_rx = num_sym_rx + 1;
            rx_sym_vec(num_sym_rx) = sym_ind;  %keep track of received symbol index
        end
    end
    rx_mat = zeros(num_sym_rx, k);
    rx_mat_2 = zeros(num_sym_rx, k);
    for ii = 1:num_sym_rx
        rx_mat(ii, :) = G(:, rx_sym_vec(ii));
        rx_mat_2(ii, :) = G_parity(:, rx_sym_vec(ii));
    end
    % Collect stats on how many symbols we receive
    rx_symbol_size_hist(num_sym_rx) = rx_symbol_size_hist(num_sym_rx) + 1;
    % Does rx_mat have at least k linearly independent rows?
    rank_1 = gfrank(rx_mat, 2);
    if (rank_1 < k)
        num_block_errors = num_block_errors + 1;
        num_symbol_errors = num_symbol_errors + (k - rank_1);
        rx_symbol_size_errors_hist(num_sym_rx) = rx_symbol_size_errors_hist(num_sym_rx) + 1;
        if (num_sym_rx >= k)
            num_symbol_errors_when_rx_at_least_k_symbols = num_symbol_errors_when_rx_at_least_k_symbols + (k - rank_1);
        end
    end
    % compare performance to second code
    rank_2 = gfrank(rx_mat_2, 2);
    if (rank_2 < k)
        num_block_errors_2 = num_block_errors_2 + 1;
        num_symbol_errors_2 = num_symbol_errors_2 + (k - rank_2);
        rx_symbol_size_errors_hist_2(num_sym_rx) = rx_symbol_size_errors_hist_2(num_sym_rx) + 1;
        if (num_sym_rx >= k)
            num_symbol_errors_when_rx_at_least_k_symbols_2 = num_symbol_errors_when_rx_at_least_k_symbols_2 + (k - rank_2);
        end
    end
end

num_symbol_errors
rx_symbol_size_hist;
num_symbol_errors_when_rx_at_least_k_symbols

percent_err = (rx_symbol_size_errors_hist./(eps + rx_symbol_size_hist))*100;
percent_err_2 = (rx_symbol_size_errors_hist_2./(eps + rx_symbol_size_hist))*100;

figure
subplot(2, 1, 1)
b = bar(percent_err);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(rx_symbol_size_hist);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title("Percent error when received X symbols.  k = 10.  Random Code")
ylim([0, 110]);
subplot(2, 1, 2)
b2 = bar(percent_err_2);
xtips2 = b2(1).XEndPoints;
ytips2 = b2(1).YEndPoints;
labels2 = string(rx_symbol_size_hist);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title("Percent error when received X symbols.  k = 10.  Random+Parity Code")
ylim([0, 110]);

