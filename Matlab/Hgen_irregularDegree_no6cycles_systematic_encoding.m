%%%%%%%%%%%%%
%% Chad Cole
%% Envistacom
%% Sept 2021
%%
%%  This program builds LDPC codes with a mostly regular degree profile and a
%%  girth of eight (no 4 or 6-cycles).  They are also lower triangular in the
%%  right part of the parity matrix, thus facilitating systematic encoding.
%%
%%  deg_v/c_prof data structure:  dimension - num_v_deg X 2
%%       first col is number of variable nodes w/ degree = second col value
%%
%%  IMPORTANT:  Always order node degrees in DESCENDING order, from 
%%               the largest deg in first row down to smallest in last.
%%%%%%%%%%%%%%%%

clear all

%%% Specify desired code variable/check profile
% deg_c_prof = [600 8];
% deg_v_prof = [1200 4];
% deg_c_prof = [102 6];
% deg_v_prof = [204 3];
% % (2000, 1500) code
% deg_c_prof = [500 12];
% deg_v_prof = [2000 3];
% % % (2000, 1000) code
% deg_c_prof = [1000 6];
% deg_v_prof = [2000 3];
% (4000, 2000) code
% deg_c_prof = [2000 6];
% deg_v_prof = [4000 3];
% (2040, 1530) code
% % deg_c_prof = [499 12; 11 11];
% deg_c_prof = [479 13; 31 12];
% deg_v_prof = [145 10; 1389 3; 476 2; 30 1];
% 4080, 3060) code
% deg_c_prof = [499 12; 11 11];
% deg_c_prof = [682 15; 338 14];
% deg_v_prof = [422 12; 2638 3; 964 2; 56 1];
deg_c_prof = [2000 6];
deg_v_prof = [4000 3];

%% sum(deg_v_prof(:,1).*deg_v_prof(:,2))
%% sum(deg_v_prof(:,1))

row_w = deg_c_prof(1, 2); % assumes regular row weight
col_w = deg_v_prof(1, 2);

num_v_deg = length(deg_v_prof(:,1));
num_c_deg = length(deg_c_prof(:,1));
n = sum(deg_v_prof(:,1));
n_k = sum(deg_c_prof(:,1));
k = n - n_k;

num_edges_v = 0;
for ii = 1:num_v_deg
    num_edges_v = num_edges_v + deg_v_prof(ii, 1)*deg_v_prof(ii, 2);
end

num_edges_c = 0;
for ii = 1:num_c_deg
    num_edges_c = num_edges_c + deg_c_prof(ii, 1)*deg_c_prof(ii, 2);
end

if(num_edges_v ~= num_edges_c)
    text = ['bad degree profile']
end
num_edges = num_edges_v;

%%%Initialize structure containing number of edges that each check node
%%%contains in the partially constructed H matrix - dc_current.  This
%%%variable gets updated for each added edges.  The dc structure, on the
%%%other hand, is fixed here and doesn't change.  It is a vector with
%%%length equal to the number of total check nodes and it contains the
%%%intended degree of each of the n-k check nodes.
dv_current = zeros(1, n);
dv = zeros(1,n);
count = 0;
for ii = 1:num_v_deg
    dv(count+1:count+deg_v_prof(ii, 1)) = deg_v_prof(ii, 2);
    count = count + deg_v_prof(ii, 1);
end

not_done = 1;
num_tries = 0;

%%% The num_count variable keeps track of how many variable nodes were able to have
%%% all of their edges connected in the graph without creating 6-cycles.
%%% If the number is close to the block length of the code, then after a
%%% practical number of tries, the desired H matrix should be created.
num_count = zeros(1,100);

while (not_done)
%%%num_tries counts how many attempts were made to create desired H matrix
num_tries = num_tries + 1;

%%% Initialize variables for each new try at creating the H matrix
temp_dv = dv;
dv_current = zeros(1, n);
count = 0;
%%% In V/Clist, the number of columns is the maximum variable/check degree
%%% plus 1 - remember deg_c/v_prof should have degrees in descending order.
Vlist=zeros(n-k,deg_c_prof(1,2)+1);  % list of V-nodes
Clist=zeros(n,deg_v_prof(1,2)+1);    % list of C-nodes

ii = 0;
ok = 1;
while (ii < (n-k-1)) && ok
    ii = ii + 1;
    count = count + 1;
    if (ii/(n-k) > 0.997)
        temp_dv = dv + 1;
    end
  %%% for each of row_w - 1 v-nodes, must find one with no 4-cycles
    v_count = 0;
           %%% Stay in this loop until all edges at current check node
           %%% have been assigned to available var nodes or else the 
           %%% constraint can't be met because there are no var node
           %%% candidates left.

           %%% Find vars (cols) which still need '1''s to satisfy
           %%% deg_v_prof constraint.
    avail_vars = find(temp_dv - dv_current);
    avail_vars = setdiff(avail_vars, (k+ii:n)); % to ensure triangle property
    already_tried_vlist = [];
    while ( length(avail_vars) > 0) && (v_count < (row_w - 1))
        
           c_nodes = [];
           v_nodes = [];
           
           unirv = rand(1,1);
            %%% Calculate the number of edges left to connect for available nodes.
%             num_edges_left = sum(temp_dv(avail_vars) - dv_current(avail_vars));
            num_edges_left = sum((temp_dv(avail_vars) - dv_current(avail_vars)).^3); % Use a square law to put more probability on low-weight vnodes
            value = ceil(num_edges_left*unirv);
            cum_sum = 0;
            d_v_index = 0;
            %%% This is the probabilistic edge assignment, assigning a
            %%% higher probability of edge connection to those var nodes
            %%% which need a larger number of connections to meet their
            %%% degree constraint.
            while (cum_sum < value)
               d_v_index = d_v_index + 1;
%                cum_sum = cum_sum + (temp_dv(avail_vars(d_v_index)) - dv_current(avail_vars(d_v_index)));
               cum_sum = cum_sum + (temp_dv(avail_vars(d_v_index)) - dv_current(avail_vars(d_v_index)))^3;
            end
            current_var = avail_vars(max(d_v_index, 1));
            no_4_cycles = 1;
            v_set = [current_var];
            c_nodes = Clist(current_var, 2:Clist(current_var, 1)+1);
            for cnode_ind = 1:length(c_nodes)
                v_set = [v_set Vlist(c_nodes(cnode_ind), 2:Vlist(c_nodes(cnode_ind), 1)+1)];
            end
            v_set = unique(v_set);
            % Want to create a new graph WITH current vnode candidate
            VlistTemp = Vlist;
            ClistTemp = Clist;
            VlistTemp(ii, 1) =  VlistTemp(ii, 1) + 1;
            VlistTemp(ii, VlistTemp(ii, 1) + 1) = current_var;
            ClistTemp(current_var, 1) = ClistTemp(current_var, 1) + 1;
            ClistTemp(current_var, ClistTemp(current_var, 1) + 1) = ii;
            
            no_4_cycles = max(0, 1-Cycle_Finder_length4_fromroot(VlistTemp, ClistTemp, current_var  ));
            no_6_cycles = max(0, 1-Cycle_Finder_length6(VlistTemp, ClistTemp, current_var));
%             for vnode_ind = 1:length(v_set)
%                 if (Cycle_Finder_length4_fromroot(VlistTemp, ClistTemp, v_set(vnode_ind) ))
%                     no_4_cycles = 0;
%                 end
%             end
            if (no_4_cycles && no_6_cycles)
                v_count = v_count + 1;
                %%% Update the code (C/Vlist) structures.
                dv_current(current_var) = dv_current(current_var) + 1;
                Vlist(ii, 1) =  Vlist(ii, 1) + 1;
                Vlist(ii, Vlist(ii, 1) + 1) = current_var;
                Clist(current_var, 1) = Clist(current_var, 1) + 1;
                Clist(current_var, Clist(current_var, 1) + 1) = ii;
            end
            %%% Once a candidate var node for edge connection has been
            %%% made, delete this var node from the available list
            already_tried_vlist = [already_tried_vlist current_var];
            avail_vars = find(temp_dv - dv_current);
            avail_vars = setdiff(avail_vars, union((k+ii:n), already_tried_vlist) ); % to ensure triangle property
    end
    if (v_count < (row_w - 1) )  % were not able to assign all row_w - 1 edges without making a 4-cycle
        ok = 0;
    end
    % Add triangle edge
    dv_current(k + ii) = dv_current(k + ii) + 1;
    Vlist(ii, 1) =  Vlist(ii, 1) + 1;
    Vlist(ii, Vlist(ii, 1) + 1) = k + ii;
    Clist(k + ii, 1) = Clist(k + ii, 1) + 1;
    Clist(k + ii, Clist(k + ii, 1) + 1) = ii;
end

%% This is to keep track of how close we got to n in each try.
    num_count(num_tries) = count

    if (ii == (n-k-1)) && (ok == 1) 
        not_done = 0;
    end
end %outer done while


%% So, the code is now contained in the Clist Vlist structures and to
%% convert this to H_sparse, the following lines are called.

H_sparse = spalloc(n-k,n,sum(Clist(:,1)) + 1);
for ii = 1:n
    for jj = 1:Clist(ii,1)
        H_sparse(Clist(ii, jj+1), ii) = 1;
    end
end
% Add final triangle edge in bottom right of H
H_sparse(n-k, n) = 1;
% Clean up any deg 1 variable nodes
for ii = k+1:n-1
    % Check for column weight 1
    if (sum(H_sparse(:, ii))==1)
        % add staircase 1
        H_sparse(ii+1-k, ii) = 1;
    end
end



%% Save your new H matrix!
%save test_2640H_irreg_1 H_sparse
% save n2040_k1530_irreg_H_no6cycles_triangleForm H_sparse
