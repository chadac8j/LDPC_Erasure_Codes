%  855        1391        1584        1642        1676

function [Is_4_cycle] = Cycle_Finder_length4_fromroot(Vlist, Clist, vroot)
            Is_4_cycle = 0;
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
            else
                Is_4_cycle = 1;
            end
