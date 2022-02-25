row_weights = sum(H_sparse, 2);
row_hist = zeros(1,20);
for ii = 1:length(row_weights)
    row_hist(row_weights(ii)) = row_hist(row_weights(ii)) + 1;
end

col_weights = sum(H_sparse, 1);
col_hist = zeros(1,20);
for ii = 1:length(col_weights)
    col_hist(col_weights(ii)) = col_hist(col_weights(ii)) + 1;
end


%%
num_6_cycle = 0;
for ii = 1:n
    if Cycle_Finder_length6(Vlist, Clist, ii)
%     if Cycle_Finder_length4_fromroot(Vlist, Clist, ii)
        str_temp = strcat(num2str(ii), "  ", num2str(Clist(ii, 2:Clist(ii, 1)+1)));
        disp(str_temp)
        num_6_cycle = num_6_cycle + 1;
    end
end
num_6_cycle


%% RS Exact FER Probability - runs into numerical stability issues for many parameters
% n_RS = 255;
% k_RS = 204; % rate = 0.8
n_RS = 150;
k_RS = 142; % rate = 0.8
P_sym_err = 0.001;

P_C = (1-P_sym_err)^(n_RS);  % initialize with all symbols correct case 
for ii = 1:(n_RS - k_RS)
    P_C = P_C + nchoosek(n_RS, ii)*P_sym_err^ii*(1-P_sym_err)^(n_RS - ii);
end
PER = 1 - P_C
