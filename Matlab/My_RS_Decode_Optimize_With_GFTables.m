% Reed Solomon Erasure Codes Decoder
% Chad Cole
% Aug 9th, 2021
%
% For erasure channel, the goal is to solve a system of linear equations in 
%  GF(2^m) to reconstruct original message.
% Need to take advantage of the fact that up to n-k erasures can be
% corrected and that the code is in systematic form, so that received
% source symbols do not need any processing, only the repair symbols.
%
% recv_vec_ind contains an ordered list of received symbol #'s.  Their
% corresponding value in GF(2^8) is in recv_vec_gf256_val
% This version takes in pre-computed GF tables for +, inv, *

function [Msg] = My_RS_Decode_Optimize_With_GFTables(recv_vec_ind, recv_vec_gf256_val, m, n, k, G, GF_add_lookup, GF_mult_lookup, GF_inv_lookup)

% build initial matrix that we will put into Gauss Jordan form. 
% It is made up of the k columns of the systematic G that we received.
    GJ_mat = zeros(k, k);
    %build initial decode matrix
    for ii = 1:k
        GJ_mat(ii, :) = G(:, recv_vec_ind(ii));
    end
    
    % Now create a num_sys_symbols X num_sys_symbols identity matrix in
    % upper left corner.  Since we assume received symbols arrive in order,
    % the first num_sys_symbols rows are systematic rows.  We need to swap
    % columns and their corresponding variable positions.
    num_sys_symbols = length(find(recv_vec_ind <= k));
    bit_order_vec = 1:k;
    sys_ind = 1;
    non_sys_ind = k;
    for ii = 1:num_sys_symbols
        col_ind = 0;
        ind = 1;
        while(col_ind == 0)
            if (GJ_mat(ii, ind) ~= 0)
                col_ind = ind;
            end
            ind = ind + 1;
        end
        temp_col = GJ_mat(:, ii);
        GJ_mat(:, ii) = GJ_mat(:, col_ind);
        GJ_mat(:, col_ind) = temp_col;
        temp_col_ind = bit_order_vec(ii);
        bit_order_vec(ii) = col_ind;
        bit_order_vec(col_ind) = temp_col_ind;
    end
    % Create zeros below diagonal
    % Only last k - num_sys_symbols of this array should be used
    repair_multiply_accumulator = recv_vec_gf256_val;
    row_index = num_sys_symbols+1;
    swap_ind = row_index + 1; %reset swap index
    NotDone = 1;
    while (row_index <= k) && (NotDone == 1)
        % First num_sys_symbols have coefficients of '1' along diagonal
        for jj = 1:num_sys_symbols
            repair_multiply_accumulator(row_index) = GF_add_lookup(repair_multiply_accumulator(row_index) + 1, GF_mult_lookup(GJ_mat(row_index, jj) + 1, repair_multiply_accumulator(jj) + 1) + 1);
            GJ_mat(row_index, jj) = 0;
        end
        for jj = num_sys_symbols+1:(row_index-1)
            repair_multiply_accumulator(row_index) = GF_add_lookup(repair_multiply_accumulator(row_index) + 1, GF_mult_lookup(GJ_mat(row_index, jj) + 1, repair_multiply_accumulator(jj) + 1) + 1);
            row_multiplier = GJ_mat(row_index, jj);  % Put this in temp variable so we don't overwrite to zero in first column and then product in all following columns as zero
            for ll = jj:k
                GJ_mat(row_index, ll) = GF_add_lookup(GJ_mat(row_index, ll) + 1, GF_mult_lookup(row_multiplier + 1, GJ_mat(jj, ll) + 1) + 1);
            end
        end
        % Check to see if diagonal element is non-zero.  If so, Make diagonal element '1' in GF(2^m)
        %  If not, need to swap with a below row
        if (GJ_mat(row_index, row_index) ~= 0)
            GF_mult = GF_inv_lookup(GJ_mat(row_index, row_index));
            for ll = row_index:k
                GJ_mat(row_index, ll) = GF_mult_lookup(GF_mult + 1, GJ_mat(row_index, ll) + 1);
            end
            repair_multiply_accumulator(row_index) = GF_mult_lookup(GF_mult + 1, repair_multiply_accumulator(row_index) + 1);
            row_index = row_index + 1;
            swap_ind = row_index + 1; %reset swap index
        else % swap rows from below
            if (swap_ind > k)
                NotDone = 0;
            else
                row_temp = GJ_mat(row_index, :);
                GJ_mat(row_index, :) = GJ_mat(swap_ind, :);
                GJ_mat(swap_ind, :) = row_temp;
                repair_multiply_accumulator_temp = repair_multiply_accumulator(row_index);
                repair_multiply_accumulator(row_index) = repair_multiply_accumulator(swap_ind);
                repair_multiply_accumulator(swap_ind) = repair_multiply_accumulator_temp;
                swap_ind = swap_ind + 1;
            end
        end
    end
    % At this point, the repair_multiply_accumulator for the final row
    % should have calculated the value for the bit_order_vec(k)th source
    % symbol
    if (row_index <= k)
        % decoding matrix is not full rank and at least one erasure remains
    end
    % Need to do Jordan operation to get to full diagonal matrix.  Start at
    % bottom and go back up to first repair row.  
    for ii = k-1:-1:num_sys_symbols+1
        for jj = (ii + 1):k
            repair_multiply_accumulator(ii) = GF_add_lookup(repair_multiply_accumulator(ii) + 1, GF_mult_lookup(repair_multiply_accumulator(jj) + 1, GJ_mat(ii, jj) + 1) + 1);
            GJ_mat(ii, jj) = 0;
        end
    end
            
    % At this point repair_multiply_accumulator should have the values of
    % the symbols that were not received as systematic, but only in a repair symbol.  
    %  These will be in the permuted order given by bit_order_vec
    final_output = zeros(1, k);
    for ii = 1:num_sys_symbols
        final_output(bit_order_vec(ii)) = recv_vec_gf256_val(ii);
    end
    for ii = (num_sys_symbols+1):k
        final_output(bit_order_vec(ii)) = repair_multiply_accumulator(ii);
    end
 
    Msg = final_output;
    