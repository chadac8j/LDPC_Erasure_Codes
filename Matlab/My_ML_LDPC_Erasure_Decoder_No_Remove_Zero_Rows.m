% LDPC Erasure decoder function

function [Msg, iterations] = My_ML_LDPC_Erasure_Decoder_No_Remove_Zero_Rows(recv_vec_val, G, k)

    iterations = 1;  % In future, can keep track of how many row operations are necessary?
% Let's represent an erasure by '-1'
    recv_vec_ind = find(recv_vec_val >= 0);
    num_recv_symbols = length(recv_vec_ind);
    GJ_mat = zeros(num_recv_symbols, k);
    recv_vec_val_no_erasures = zeros(1, num_recv_symbols);
    %build initial decode matrix
    for ii = 1:num_recv_symbols
        GJ_mat(ii, :) = G(:, recv_vec_ind(ii));
        recv_vec_val_no_erasures(ii) = recv_vec_val(recv_vec_ind(ii));
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
    % So we have num_recv_symbols - num_sys_symbols rows to try to
    % reconstruct an identity matrix.  Many of these will probably be
    % linearly dependent
    
    repair_multiply_accumulator = recv_vec_val_no_erasures;
    row_index = num_sys_symbols+1;
    swap_ind = row_index + 1; %reset swap index
    NotDone = 1;
    while (row_index <= k) && (NotDone == 1)
        % First num_sys_symbols have coefficients of '1' along diagonal
        for jj = 1:num_sys_symbols
            if GJ_mat(row_index, jj) == 1
                repair_multiply_accumulator(row_index) = mod(repair_multiply_accumulator(row_index) + repair_multiply_accumulator(jj), 2);                
                GJ_mat(row_index, jj) = 0;
            end
        end
        for jj = num_sys_symbols+1:(row_index-1)
            if (GJ_mat(row_index, jj) == 1)
                repair_multiply_accumulator(row_index) = mod(repair_multiply_accumulator(row_index) + repair_multiply_accumulator(jj), 2);
                GJ_mat(row_index, jj:k) = mod(GJ_mat(row_index, jj:k) + GJ_mat(jj, jj:k), 2);
            end
        end
        % Check to see if diagonal element is non-zero (since binary, has to be '1').
        %  If not, need to swap with a below row
        if (GJ_mat(row_index, row_index) ~= 0) % since binary, must be '1'
            row_index = row_index + 1;
            swap_ind = row_index + 1; %reset swap index
        else % swap rows from below
            if (swap_ind > num_recv_symbols) % we have run out of rows to choose from
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
    % At this point, the new matrix should be full rank in first k rows (if
    % possible)
    if (row_index <= k)
        % decoding matrix is not full rank and at least one erasure remains
        Msg = recv_vec_val;  % Don't do any decoding if we know we can't succeed?
        return 
    end
    % Need to do Jordan operation to get to full diagonal matrix.  Start at
    % bottom and go back up to first repair row.  
    for ii = k-1:-1:num_sys_symbols+1
        for jj = (ii + 1):k
            repair_multiply_accumulator(ii) = mod(repair_multiply_accumulator(ii) + repair_multiply_accumulator(jj)*GJ_mat(ii, jj), 2);
            GJ_mat(ii, jj) = 0;
        end
    end
            
    % At this point repair_multiply_accumulator should have the values of
    % the symbols that were not received as systematic, but only in a repair symbol.  
    %  These will be in the permuted order given by bit_order_vec
    final_output = zeros(1, k);
    for ii = 1:num_sys_symbols
        final_output(bit_order_vec(ii)) = recv_vec_val_no_erasures(ii);
    end
    for ii = (num_sys_symbols+1):k
        final_output(bit_order_vec(ii)) = repair_multiply_accumulator(ii);
    end


    Msg = final_output;

