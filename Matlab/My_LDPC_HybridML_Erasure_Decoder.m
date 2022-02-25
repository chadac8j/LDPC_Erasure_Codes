% LDPC Erasure decoder function

function [Msg, iterations] = My_LDPC_HybridML_Erasure_Decoder(recv_vec_val, Vlist, Clist, H_sparse, n)

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
            y_current(erasure_ind) = mod(sum(y_current(setxor(Vlist(ii, 2:Vlist(ii,1)+1), erasure_ind))), 2);
        end
    end

    num_cur_erasures = length(find(y_current == -1));
    if(num_cur_erasures == 0)
        stopsig = 1;
    end

%Keep a history of erasures for the current decoding iteration.
    erasure_hist(itestep) = num_cur_erasures;

end

if (num_cur_erasures > 0) % now do ML decoding on residual erasures
    % Find H_erasure = H_known columns
    erasure_ind = find(y_current == -1);
    num_erasures = length(erasure_ind);
    find_inv = H_sparse(:, erasure_ind);
    non_erasure_ind = setdiff(1:n, erasure_ind);
    rhs = mod(H_sparse(:, non_erasure_ind)*y_current(non_erasure_ind)', 2);
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
        % zero out other non-zero rows below diagonal
        for ii = 2:length(non_zero_ind)
            find_inv(non_zero_ind(ii), :) = mod(find_inv(non_zero_ind(ii), :) + find_inv(col, :), 2);
            rhs(non_zero_ind(ii)) = mod(rhs(non_zero_ind(ii)) + rhs(col), 2);
        end
    end
    % Now do Jordan elimination on upper triangle
    if ~dont_do_jordan
        for col=num_erasures:-1:2
            non_zero_ind = find(find_inv(1:col-1, col));
            % zero out other non-zero rows above diagonal
            for ii = 1:length(non_zero_ind)
                find_inv(non_zero_ind(ii), :) = mod(find_inv(non_zero_ind(ii), :) + find_inv(col, :), 2);
                rhs(non_zero_ind(ii)) = mod(rhs(non_zero_ind(ii)) + rhs(col), 2);
            end
        end
    end
    y_current(erasure_ind) = rhs(1:num_erasures);
end
Msg = y_current;
iterations = itestep;
