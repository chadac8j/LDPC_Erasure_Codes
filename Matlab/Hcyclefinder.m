%%%%%%%%%%%%%
%% Hcyclefinder.m - Requires an H_sparse as input.
%%
%%  This utility searches LDPC codes for cycles of up to length 8.  It
%%  counts the number of each type of cycle.  It could be useful to compare
%%  the cycle distributions of two equal-parameter codes and select the one
%%  with the least number of cycles.
%%
%%%%%%%%%%%%%%%%


% clear all
% function sumout = Hcyclefinder ()

% load H1008_nox2
% load H1200_4_8_no6cycle
% load H2048_8023an

n = length(H_sparse(1,:));
n_k = length(H_sparse(:,1));
% n_k = gfrank(full(H_sparse));
k = n - n_k;
rate = k/n;

%% The num_x_cycles length-n vector contains the number of cycles
%% encountered for the tree rooted at each variable nodes.
no_4_cycle = 1;
num_4_cycles = zeros(1,n);
no_6_cycle = 1;
num_6_cycles = zeros(1,n);
no_8_cycle = 1;
num_8_cycles = zeros(1,n);

Vlist=zeros(n-k,max(sum(H_sparse,2))+1);  % list of V-nodes
Clist=zeros(n,max(sum(H_sparse,1))+1);    % list of C-nodes

%Generate C/Vlist from H_sparse data.
for jj=1:n-k
    Vlist(:,1)=sum(H_sparse,2);
    icnt=0;
    for ii=1:n
        if H_sparse(jj,ii)==1
            icnt=icnt+1;
            Vlist(jj,icnt+1)=ii;
        end
    end
end

for ii=1:n
    Clist(:,1)=(sum(H_sparse,1))';
    jcnt=0;
    for jj=1:n-k
        if H_sparse(jj,ii)==1
            jcnt=jcnt+1;
            Clist(ii,jcnt+1)=jj;
        end
    end
end

t = 0;
while (t < n)
    t = t + 1;
    c_tier_1 = Clist(t, 2:Clist(t,1)+1);
    v_tier_1 = [];
    v_tier_1_ind = [];
    vtemp = [];
    %%% Find all variable nodes in tier 1 and their corresponding
    %%% connections to check nodes in tier 1.
    for ii = 1:length(c_tier_1)
        v_tier_1 = [v_tier_1, Vlist(c_tier_1(ii), 2:Vlist(c_tier_1(ii),1)+1)];
        v_tier_1_ind = [v_tier_1_ind, c_tier_1(ii)*ones(1,Vlist(c_tier_1(ii),1))];
    end
        index = 0;
    %Get rid of the duplicates of the root v-node in tier 1
    for jj = 1:length(v_tier_1)
        if (v_tier_1(jj) ~= t)
            index = index + 1;
            vtemp(index) = v_tier_1(jj);
            vtemp_ind(index) = v_tier_1_ind(jj);
        end
    end
    v_tier_1 = vtemp;
    v_tier_1_ind = vtemp_ind;

    %%% Search for 4-cycles.  Order variable nodes and find duplicates.
    [ordered, ind] = sort(v_tier_1);
    for ii = 1:length(ordered)-1
        if ordered(ii) == ordered(ii+1)
            no_4_cycle = 0;  % There are 4-cycles.
            num_4_cycles(t) = num_4_cycles(t) + 1;  %% Count number of 4-cycles.
            %build a list of all 4-cycle participants
            var_temp_1 = v_tier_1(ind(ii));
% Add the cycle bits to the list, it's understood that the root is present
            dont_add = 0;
        end
    end
    %now start to look at 6-cycles
    cindex = 0;
    c_tot = 0;
    for ii = 1:length(v_tier_1)
        %Second row of c_tier_2 gives info about connection to upper level of tree
        c_tier_2(2,c_tot+1:c_tot+Clist(v_tier_1(ii),1)-1) = v_tier_1(ii);
        %Find check nodes in check tier two
        throwaway_checks = intersect(Clist(v_tier_1(ii), 2:Clist(v_tier_1(ii),1)+1),v_tier_1_ind(ii));
        add_checks = setxor(throwaway_checks, Clist(v_tier_1(ii), 2:Clist(v_tier_1(ii),1)+1));
        c_tier_2(1,cindex+1:cindex+length(add_checks)) = add_checks;
        cindex = cindex + length(add_checks);
        c_tot = c_tot + Clist(v_tier_1(ii), 1)-1;
    end
    [ordered, ind] = sort(c_tier_2(1,:));
    for ii = 1:length(ordered)-1
        %%% This loop looks for duplicate check nodes at check tier two.
        %%% If any are found, the following if statement is true.
        if ordered(ii) == ordered(ii+1)
            no_6_cycle = 0;
            num_6_cycles(t) = num_6_cycles(t) + 1;
            c_temp_1 = c_tier_2(1,ind(ii));
            c_temp_2 = c_tier_2(1,ind(ii+1));
            v_parent_1 = c_tier_2(2,ind(ii));
            v_parent_2 = c_tier_2(2,ind(ii+1));
        end
    end
    vindex = 0;
    c_tot = 0;
    for ii = 1:length(c_tier_2(1,:))
        v_tier_2(2,c_tot+1:c_tot+Vlist(c_tier_2(1, ii),1)-1) = c_tier_2(1,ii);
        for jj = 1:Vlist(c_tier_2(1, ii),1)
            if Vlist(c_tier_2(1,ii), jj+1) ~= c_tier_2(2,ii)
                vindex = vindex + 1;
                v_tier_2(1,vindex) = Vlist(c_tier_2(1,ii), jj+1);
            end
        end
        c_tot = c_tot + Vlist(c_tier_2(1, ii), 1)-1;
    end
    [ordered, ind] = sort(v_tier_2(1,:));
    for ii = 1:length(ordered)-1
        %%% This loop looks for duplicate variable nodes at variable tier two.
        %%% If any are found, the following if statement is true.
        if ordered(ii) == ordered(ii+1)
            no_8_cycle = 0;
            num_8_cycles(t) = num_8_cycles(t) + 1;
        end
    end
end %while

%sumout = sum(num_6_cycles)

