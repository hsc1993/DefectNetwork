function [rn,links,links_newloop,links_start_end_list] = addnewloop(nlist,b,numSegs,shift,rn,links,links_start_end_list,existing_loops,seglength_multiply)
    % links0 is the set of links associated with new loop

    [rn1,links1] = generatePDL(nlist,b,numSegs,shift,existing_loops,seglength_multiply);
    [rn,links,links0] = addPDLtolist(rn,links,rn1,links1);

    links_newloop = links0;
    if size(links_start_end_list,1) == 0
        links_start_end_list = [links_start_end_list;1,6*numSegs];
    else
        links_start_end_list = [links_start_end_list;links_start_end_list(end,2)+1,links_start_end_list(end,2)+6*numSegs];
    end

end

