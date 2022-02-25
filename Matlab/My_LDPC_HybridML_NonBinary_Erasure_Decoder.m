% LDPC Erasure decoder function
% This version works over non-binary GF(Q)

function [Msg, iterations] = My_LDPC_HybridML_NonBinary_Erasure_Decoder(recv_vec_val, Vlist, Clist, H_sparse, n, k, GF_add_lookup, GF_mult_lookup, GF_inv_lookup)

do_ML_decode = 1;
% do_ML_decode = 0;

% Let's represent an erasure by '-1'
y_current = recv_vec_val;


itenum = 10;
% itenum = 50;

stopsig=0;  % If a valid codeword has been found, set this to `1'
itestep=0;  % Iteration number

erasure_hist = zeros(1, itenum);

while (stopsig==0) && (itestep<itenum)

    itestep=itestep+1;

%%%%%% Algorithm will first look for any c-nodes connected to exactly one
%%%%%% erasure
    for ii = 1:length(Vlist(:, 1))
        num_erasures = 0;
        erasure_ind = 0;
        for jj = 1:Vlist(ii, 1)
            if (y_current(Vlist(ii, jj+1)) == -1)
                num_erasures = num_erasures + 1;
                erasure_ind = Vlist(ii, jj+1);
            end
        end
            % If we can correct the erasure
        if (num_erasures == 1)
%             y_current(erasure_ind) = mod(sum(y_current(setxor(Vlist(ii, 2:Vlist(ii,1)+1), erasure_ind))), 2);
            check_indices = setxor(Vlist(ii, 2:Vlist(ii,1)+1), erasure_ind);
            gf_sum = 0;
            for kk = 1:length(check_indices)
                if (y_current(check_indices(kk)) > 255) || (y_current(check_indices(kk)) < 0)
                   poop = 1; 
                end
                gf_sum = GF_add_lookup(gf_sum + 1, GF_mult_lookup(y_current(check_indices(kk)) + 1, H_sparse(ii, check_indices(kk)) + 1) + 1);
            end
            y_current(erasure_ind) = GF_mult_lookup(gf_sum + 1, GF_inv_lookup(H_sparse(ii, erasure_ind)) + 1);
        end
    end

    num_cur_erasures = length(find(y_current == -1));
    if(num_cur_erasures == 0)
        stopsig = 1;
    end

%Keep a history of erasures for the current decoding iteration.
    erasure_hist(itestep) = num_cur_erasures;

end

if ((num_cur_erasures > 0) && (do_ML_decode == 1)) % now do ML decoding on residual erasures
    % Find H_erasure = H_known columns
    erasure_ind = find(y_current == -1);
    num_erasures = length(erasure_ind);
    find_inv = H_sparse(:, erasure_ind);
    
%     rank_rx = rank(full(find_inv));
%     if (rank_rx < num_erasures)
%         not_full_rank = 1
%     end

    non_erasure_ind = setdiff(1:n, erasure_ind);
%     rhs = mod(H_sparse(:, non_erasure_ind)*y_current(non_erasure_ind)', 2);
    rhs = zeros(n-k, 1);
    for kk = 1:(n-k)
        non_erasure_ind_kk = intersect(Vlist(kk, 2:Vlist(kk,1)+1), non_erasure_ind);
        gf_sum = 0;
        for ll = 1:length(non_erasure_ind_kk)
            gf_sum = GF_add_lookup(gf_sum + 1, GF_mult_lookup(y_current(non_erasure_ind_kk(ll)) + 1, H_sparse(kk, non_erasure_ind_kk(ll)) + 1) + 1);
        end
        rhs(kk) = gf_sum;
    end
    dont_do_jordan = 0;
    % Need to solve linear equation find_inv*x = rhs
    for col=1:num_erasures
        non_zero_ind = find(find_inv(col:end, col)) + col - 1;
        if length(non_zero_ind) == 0
            dont_do_jordan = 1;
            break; % linearly dependent columns within erasure space
        end
        %swap first non-zero row
        temp_val = rhs(col);
        rhs(col) = rhs(non_zero_ind(1));
        rhs(non_zero_ind(1)) = temp_val;
        temp_row = find_inv(col, :);
        find_inv(col, :) = find_inv(non_zero_ind(1), :);
        find_inv(non_zero_ind(1), :) = temp_row;
        % multiply by inverse to make 1's in diagonal
        non_zero_indices = find(find_inv(col, :));
        multiplier = GF_inv_lookup(find_inv(col, non_zero_indices(1)));
        find_inv(col, col) = GF_mult_lookup(find_inv(col, col) + 1, multiplier + 1);  % this should equal 1
        rhs(col) = GF_mult_lookup(rhs(col) + 1, multiplier + 1);
        for kk = 2:length(non_zero_indices)
            find_inv(col, non_zero_indices(kk)) = GF_mult_lookup(find_inv(col, non_zero_indices(kk))+1, multiplier + 1);
        end
        % zero out other non-zero rows below diagonal
        for ii = 2:length(non_zero_ind)
            non_zero_row_indices = union(find(find_inv(non_zero_ind(ii), :)), find(find_inv(col, :)));
            multiplier = find_inv(non_zero_ind(ii), col);
            for ll = 1:length(non_zero_row_indices)
                find_inv(non_zero_ind(ii), non_zero_row_indices(ll)) = GF_add_lookup(find_inv(non_zero_ind(ii), non_zero_row_indices(ll)) + 1, GF_mult_lookup(multiplier + 1,find_inv(col, non_zero_row_indices(ll))  + 1) + 1);
            end
            rhs(non_zero_ind(ii)) = GF_add_lookup(rhs(non_zero_ind(ii)) + 1, GF_mult_lookup(multiplier + 1, rhs(col) + 1) + 1);
        end
    end
    % Now do Jordan elimination on upper triangle
    if ~dont_do_jordan
        for col=num_erasures:-1:2
            non_zero_ind = find(find_inv(1:col-1, col));
            % zero out other non-zero rows above diagonal
            for ii = 1:length(non_zero_ind)
                rhs(non_zero_ind(ii)) = GF_add_lookup(rhs(non_zero_ind(ii)) + 1, GF_mult_lookup(find_inv(non_zero_ind(ii), col) + 1, rhs(col) + 1) + 1);
                find_inv(non_zero_ind(ii), col) = 0;
            end
        end
    end
    y_current(erasure_ind) = rhs(1:num_erasures);
end
Msg = y_current;
iterations = itestep;
