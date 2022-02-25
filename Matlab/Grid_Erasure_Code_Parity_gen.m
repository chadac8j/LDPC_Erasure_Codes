% Grid code parity matrix

% Take grid of variable nodes row by row to create source code vec, i.e.
% x1 x2 x3
% x4 x5 x6
% x7 x8 x9
% becomes
% x1 x2 x3 x4 x5 x6 x7 x8 x9


% Say 14 X 14 code
% column_height = 14;
% row_width = 14;
column_height = 10;
row_width = 5;
k = column_height*row_width;
n_k = column_height + row_width; % one parity check per row and column
n = n_k + k;
rate = k/n

H_sparse = spalloc(n-k,n,2*column_height*row_width);
% create row parity first
for ii = 1:column_height
        H_sparse(ii, (ii-1)*row_width+1:ii*row_width) = 1;
        H_sparse(ii, k+ii) = 1; %parity bit
end

% now create column parity
for ii = (column_height+1):(column_height + row_width)
    for jj = 1:column_height
        H_sparse(ii, (ii - column_height) + (jj-1)*row_width) = 1;
    end
    H_sparse(ii, k+ii) = 1; %parity bit
end
