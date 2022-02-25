% LDPC Erasure decoder function

function [Msg, iterations] = My_LDPC_Erasure_Decoder(recv_vec_val, Vlist, Clist)

% Let's represent an erasure by '-1'
y_current = recv_vec_val;


% itenum = 100;
itenum = 50;
% itenum = 20;

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

Msg = y_current;
iterations = itestep;
