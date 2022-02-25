%%%%%%%%%%%%%
%% Hgen_no6cycles.m - Requires no input files.
%%
%%  This program builds LDPC codes with arbitrary degree profile and a
%%  girth of eight (no 6-cycles) 
%%
%%  deg_v/c_prof data structure:  dimension - num_v_deg X 2
%%       first col is number of variable nodes w/ degree = second col value
%%
%%  IMPORTANT:  Always order node degrees in DESCENDING order, from 
%%               the largest deg in first row down to smallest in last.
%%%%%%%%%%%%%%%%

clear all
%format compact
% function Hgen_rand_no6cycles_fun ()

%%% Specify desired code variable/check profile
% deg_c_prof = [600 8];
% deg_v_prof = [1200 4];
% deg_c_prof = [102 6];
% deg_v_prof = [204 3];
deg_c_prof = [500 6];
deg_v_prof = [1000 3];
%deg_c_prof = [302 9];
%deg_v_prof = [54 20; 21 7; 7 6; 92 5; 133 3; 295 2];
% deg_v_prof = [60 8;241 4; 302 2];
% deg_c_prof = [242 7; 59 6];
% deg_c_prof = [252 8];
% deg_v_prof = [504 4];
% deg_c_prof = [16 16];
% deg_v_prof = [128 2];
% deg_c_prof = [408 6];
% deg_v_prof = [816 3];
% deg_c_prof = [5000 6];
% deg_v_prof = [10000 3];
% deg_c_prof = [820 15];
% deg_v_prof = [4100 3];
% deg_c_prof = [1320 6];
% deg_v_prof = [2640 3];
% deg_c_prof = [1000 15];
% deg_v_prof = [5000 3];
num_v_deg = length(deg_v_prof(:,1));
num_c_deg = length(deg_c_prof(:,1));
n = sum(deg_v_prof(:,1));
k = n - sum(deg_c_prof(:,1));

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
dc_current = zeros(1, n-k);
dc = zeros(1,n-k);
count = 0;
for ii = 1:num_c_deg
    dc(count+1:count+deg_c_prof(ii, 1)) = deg_c_prof(ii, 2);
    count = count + deg_c_prof(ii, 1);
end

%%% This variable specifys how many total edges still need to be connected
%%% in the partially constructed H matrix.
num_edges_left = num_edges - sum(dc_current);


not_done = 1;
num_tries = 0;
no_repeats = 0;
num_picked = 0;

%%% The num_count variable keeps track of how many variable nodes were able to have
%%% all of their edges connected in the graph without creating 6-cycles.
%%% If the number is close to the block length of the code, then after a
%%% practical number of tries, the desired H matrix should be created.
num_count = zeros(1,100);

while (not_done)
%%%num_tries counts how many attempts were made to create desired H matrix
num_tries = num_tries + 1;

%%% Initialize variables for each new try at creating the H matrix
dc_current = zeros(1, n-k);
num_edges_left = num_edges;
count = 0;
%%% In V/Clist, the number of columns is the maximum variable/check degree
%%% plus 1 - remember deg_c/v_prof should have degrees in descending order.
Vlist=zeros(n-k,deg_c_prof(1,2)+1);  % list of V-nodes
Clist=zeros(n,deg_v_prof(1,2)+1);    % list of C-nodes

ii = 0;
ok = 1;
while (ii < num_v_deg) && ok
    ii = ii + 1;
    jj = 0;
    while (jj < deg_v_prof(ii, 1)) && ok
           jj = jj + 1;
           count = count + 1;
           
           %%% c_nodes will contain all of the check nodes at the 2nd check
           %%% tier in the tree with current variable node as root.  As
           %%% each edge connection is specified for a given variable node
           %%% as root, this list will grow.
           c_nodes = [];
           
           v_nodes = [];
           %%% Find checks (rows) which still need '1''s to satisfy
           %%% deg_c_prof constraint.
           avail_checks = find(dc - dc_current);
           kk = 0;
           %%% Stay in this loop until all edges at current variable node
           %%% have been assigned to available check nodes or else the 
           %%% constraint can't be met because there are no check node
           %%% candidates left.
           while (length(avail_checks) > 0) && (kk < deg_v_prof(ii, 2))
            unirv = rand(1,1);
            %%% Calculate the number of edges left to connect in graph.
            num_edges_left = sum(dc(avail_checks) - dc_current(avail_checks));
            value = ceil(num_edges_left*unirv);
            cum_sum = 0;
            d_c_index = 0;
            %%% This is the probabilistic edge assignment, assigning a
            %%% higher probability of edge connection to those check nodes
            %%% which need a larger number of connections to meet their
            %%% deg_c_prof constraint.
            while (cum_sum < value)
               d_c_index = d_c_index + 1;
               cum_sum = cum_sum + (dc(avail_checks(d_c_index)) - dc_current(avail_checks(d_c_index)));
            end
            current_check = avail_checks(d_c_index);
            %%% Once a candidate check node for edge connection as been
            %%% made, delete this check node from the available list for
            %%% the next edge connections of this variable node.
            avail_checks = setxor(avail_checks,current_check);

            %%% This is the list of current v-nodes in first tier, not
            %%% including the root.
            cycle_nodes_tent = Vlist(current_check, 2:Vlist(current_check,1)+1);
            c_nodes_tent = [current_check];
            for qq = 1:Vlist(current_check,1)
                c_nodes_tent = [c_nodes_tent, setxor(current_check, Clist(Vlist(current_check, qq+1), 2:Clist(Vlist(current_check, qq+1), 1)+1))];
            end

            c_nodes_tent = [c_nodes_tent, c_nodes];
        %%% Ensuring the first tier of v-nodes is unique will give
        %%% girth 6.  For a girth of 8 we must also make sure no check
        %%% nodes in check tier 2 are duplicates.  If the following if
        %%% statement is true, then the current edge does not create any
        %%% 6-cycles.
            if (length(unique(c_nodes_tent)) == length(c_nodes_tent))
                kk = kk + 1;  %move on to adding next edge
                %%% Update the code (C/Vlist) structures.
                c_nodes = c_nodes_tent;
                dc_current(current_check) = dc_current(current_check) + 1;
                Vlist(current_check, 1) =  Vlist(current_check, 1) + 1;
                Vlist(current_check, Vlist(current_check, 1) + 1) = count;
                Clist(count, 1) = Clist(count, 1) + 1;
                Clist(count, Clist(count, 1) + 1) = current_check;
            end
           end
           if (kk < deg_v_prof(ii, 2))
               ok = 0;
           end
           
    end
end

%% The following checks whether the check node (deg_c_prof) constraint is
%% valid for this code.  If it is not imperative that the exact check
%% profile is as specified, then this extra constraint can be commented
%% out.
  all_rows_good = 1;
  count_i = 0;
  for iii = 1:num_c_deg
     for jjj = 1:deg_c_prof(iii, 1)
        count_i = count_i + 1;
        if(Vlist(count_i,1) ~= deg_c_prof(iii, 2))
          all_rows_good = 0;
        end
     end
  end
  if (all_rows_good == 1)
       not_done = 0;
  end
  %% If check node degree profile does not need to be explicitly followed,
  %% then uncomment the following three lines.
    if (ok == 1)
        not_done = 0;
    end

%% This is to keep track of how close we got to n in each try.
    num_count(num_tries) = count
end %outer done while

%% So, the code is now contained in the Clist Vlist structures and to
%% convert this to H_sparse, the following lines are called.

H_sparse = spalloc(n-k,n,sum(Clist(:,1)));
for ii = 1:n
    for jj = 1:Clist(ii,1)
        H_sparse(Clist(ii, jj+1), ii) = 1;
    end
end

%% Save your new H matrix!
%save test_2640H_irreg_1 H_sparse

