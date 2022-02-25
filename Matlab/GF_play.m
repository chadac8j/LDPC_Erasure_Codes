% M = [1, 1; 0, 1];
% dim = size(M);
% rows = dim(1);
% I_2 = eye(2);
% 
% for i=1:rows
%    M_inv(1:rows,i)=gflineq(M,I_2(1:rows,i),2);
% end
% mod(M_inv*M, 2)

% row_weights = sum(G_parity(:, k+1:n), 2);
% for ii = 1:length(row_weights)
%     if(row_weights(ii) < 2)
%         indices = find(G_parity(ii, k+1:n)==0);
%         %exchange with row that has most ones already
%         swap_col = k + indices(1);
%         col_ind = find(G_parity(:, swap_col));
%         [val, ind] = max(row_weights(col_ind));
%         %swap rows
%         G_parity(col_ind(ind), swap_col) = 0;
%         G_parity(ii, swap_col) = 1;        
%     end
% end

m = 8;
GF_SIZE = 2^m;

% prim_poly_m8 = [1 0 1 1 1 0 0 0 1]; % prim poly given in rfc 5510
prim_poly_m8 = fliplr([1 1 1 0 0 0 0 1 1]); % prim poly given in TIA 5041
prim_poly_m8_number = 0;
for ii = 1:length(prim_poly_m8)
    prim_poly_m8_number = prim_poly_m8_number + prim_poly_m8(length(prim_poly_m8) - ii + 1)*2^(ii-1);
end
% Form log table
gftable(m, prim_poly_m8_number);
x_gf = gf(0:2^m-1, m, prim_poly_m8_number);


% Form log table
gf_log_inv = zeros(1, 2^m-1);
gf_log_inv(1) = 0;
gf_log_inv(2) = 1;
alpha = x_gf(3);
gf_temp = alpha;
for ii = 3:2^m
    gf_log_inv(ii) = double(gf_temp.x);
    gf_temp = alpha*gf_temp;
end
%Appears to be a bug in the log() func for GF
Log_array = [-inf 0:n-1];
log_lookup(gf_log_inv+1) = Log_array;
% gf_log_inv = gf(gf_log_inv, m, prim_poly_m8_number);

gf_inv = zeros(1, GF_SIZE);
gf_inv(1) = 0;
gf_inv(2) = 1;
for ii = 3:GF_SIZE
	gf_inv(ii) = gf_log_inv(GF_SIZE-log_lookup(ii)+1);
end
gf_inv = gf_inv(2:end);

% Say we want X * Y, so take log_inv(log(X) + log(Y))

X = alpha;
Y = alpha^2;
log_eq_ind = gf_log_inv(mod(log_lookup(X.x + 1) + log_lookup(Y.x + 1), 255) + 2);

% % The following code is used in Build_GF256_Lookup_Tables()
% % Unit test for my GF(256) lookup tables vs Matlab built in implementation
% my_GF_mult_table = zeros(GF_SIZE, GF_SIZE);
% for row = 2:length(x_gf)
%     for col = 2:length(x_gf)
%         matlab_imp = x_gf(row)*x_gf(col);
%         my_imp = gf_log_inv(mod(log_lookup(double(x_gf.x(row)) + 1) + log_lookup(double(x_gf.x(col)) + 1), 255) + 2);
%         my_GF_mult_table(row, col) = my_imp;
%         if (matlab_imp.x ~= my_imp)
%             my_imp
%         end
%     end
% end
% 
% % Unit test against Matlab built-in function for inverse of a GF()
% for ii = 2:GF_SIZE
%     matlab_imp = inv(x_gf(ii));
%     my_imp = gf_inv(x_gf.x(ii) );
%     if (matlab_imp.x ~= my_imp)
%             my_imp
%     end
% end
% 
% % Unit test for bitxor addition/subtraction logic vs Matlab built in GF implementation
% my_GF_add_table = zeros(GF_SIZE, GF_SIZE);
% for row = 1:length(x_gf)
%     for col = 1:length(x_gf)
%         matlab_imp = x_gf(row) + x_gf(col);
%         my_imp = bitxor( double(x_gf.x(row)), double(x_gf.x(col)) );
%         my_GF_add_table(row, col) = my_imp;
%         if (matlab_imp.x ~= my_imp)
%             my_imp
%         end
%     end
% end

% Testing some properties of expanded RS code in TIA-5041 standard
alpha_4 = [1 0 1 0 1 1 0 1;
           1 0 0 1 0 1 0 1;
           1 0 0 0 1 0 0 1;
           1 0 0 0 0 1 1 1;
           1 0 0 0 0 0 0 0;
           0 1 0 0 0 0 0 0;
           0 0 1 0 0 0 0 0;
           0 0 0 1 0 0 0 0];
       
alpha_12 = zeros(m, m);
for ii = 1:m
    alpha_12(m-ii+1, :) = mydec2binvec(gf_log_inv(12+ii), 8);
end

alpha_24 = mod(alpha_12*alpha_12, 2)