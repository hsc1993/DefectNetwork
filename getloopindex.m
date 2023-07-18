function [m1,m2] = getloopindex(idx1,idx2,links_start_end_list)
    % if currently considered segment belongs to a dislocation that has
    % reacted, freeze it
    
    for ii = 1:size(links_start_end_list,1)
        if idx1>=links_start_end_list(ii,1) && idx1<=links_start_end_list(ii,2)
            m1 = ii; % find out which loop (index m1) does the current node belong to
            continue
        end
    end
    
    for jj = 1:size(links_start_end_list,1)
        if idx2>=links_start_end_list(jj,1) && idx2<=links_start_end_list(jj,2)
            m2 = jj;
            continue
        end
    end

end