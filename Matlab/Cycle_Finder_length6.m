
function [Num_6_cycles] = Cycle_Finder_length6(Vlist, Clist, vroot)
            Num_6_cycles = 0;
            vnodes_tent = [vroot];
            for qq = 1:Clist(vroot, 1)
                vnodes_tent = [vnodes_tent, setxor(vroot, Vlist(Clist(vroot, qq+1), 2:Vlist(Clist(vroot, qq+1), 1)+1))];
            end

        %%% Ensuring the first tier of v-nodes is unique will give
        %%% girth 6.  For a girth of 8 we must also make sure no check
        %%% nodes in check tier 2 are duplicates.  If the following if
        %%% statement is true, then the current edge does not create any
        %%% 6-cycles.
            if (length(unique(vnodes_tent)) == length(vnodes_tent))
%                 disp("No 4 cycles")
                % We have established that there are no 4-cycles.
                %  Now look for 6-cycles
                %           v
                %         / | \
                %        c  c  c   tier 1
                %      / | \
                %     v  v  v      tier 2
                %         / | \
                %        c  c  c   tier 3

                    %now start to look at 6-cycles
                c_tier_1 = Clist(vroot, 2:Clist(vroot,1)+1);
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
                vtemp = [];
                vtemp_ind = [];
    %Get rid of the duplicates of the root v-node in tier 1
                for jj = 1:length(v_tier_1)
                    if (v_tier_1(jj) ~= vroot)
                        index = index + 1;
                        vtemp(index) = v_tier_1(jj);
                        vtemp_ind(index) = v_tier_1_ind(jj);
                    end
                end
                v_tier_1 = vtemp;
                v_tier_1_ind = vtemp_ind;
                cindex = 0;
                c_tot = 0;
                c_tier_2 = zeros(2,0);
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
                        Num_6_cycles = 1;
                    end
                end

            
            else  % There is a 4-cycle
                Num_6_cycles = 1;
            end

