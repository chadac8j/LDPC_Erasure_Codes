% This is a script to test My_RS_Decode

% Generate test matrix

m = 8;
% n = 2^m-1;
n = 7;
num_symbols_to_correct = 2;
k = n - num_symbols_to_correct;
log_lookup = []; %Maybe use later?

prim_poly_m8 = [1 0 1 1 1 0 0 0 1];
prim_poly_m8_number = 0;
for ii = 1:length(prim_poly_m8)
    prim_poly_m8_number = prim_poly_m8_number + prim_poly_m8(length(prim_poly_m8) - ii + 1)*2^(ii-1);
end
% Form log table
gftable(m, prim_poly_m8_number);
x_gf = gf(0:2^m-1, m, prim_poly_m8_number);

x_gf = gf(0:2^m-1, m, prim_poly_m8_number);
alpha = x_gf(3);

G = gf(zeros(k, n), m, prim_poly_m8_number);
% G(1:k,1:k) = eye(5);
% G(:,k+1) = [1 0 1 1 0]';
% G(:,k+2) = [0 1 0 0 1]';
% G(:,k+3) = [0 1 1 0 1]';

for row=1:k
    for col=1:n
        G(row, col) = alpha^(col*row);
    end
end
% Eventually want a systematic version of G
G_k_inv = inv(G(1:k, 1:k));
G = G_k_inv*G; % Find systematic G

% G_rx = [G(:, 1) G(:, 3) G(:, 5) G(:, 6) G(:, 7)];
% G_rx = [G(:, 1) G(:, 3) G(:, 5) G(:, 6) G(:, 8)];
% recv_vec_ind = [1 3 5 6 7];
num_trials = 100000;
num_errors = 0;

for iter = 1:num_trials
    source_vec = round(255*rand(1, k)); %create source symbols over GF(2^8)
%     source_vec = round(rand(1, k));
    source_encode_vec = source_vec*G;
    recv_vec_ind = randperm(n);
    recv_vec_ind = sort(recv_vec_ind(1:k));
    rx_vec = source_encode_vec(recv_vec_ind);

    Out_RS = My_RS_Decode(recv_vec_ind, rx_vec, m, n, k, prim_poly_m8_number, G, log_lookup);
    if sum(Out_RS == source_vec) ~= k
        num_errors = num_errors + 1;
    end
       
end


