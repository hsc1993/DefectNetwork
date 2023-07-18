function reacted_node = detect_reaction(links)
    
    node_id = [];
    for i = 1:size(links,1)
        node_id = [node_id;links(i,1);links(i,2)];
    end
    unique_node = unique(node_id);
    
    flag_reaction = 2<histc(node_id,unique(node_id));
    
    if all(flag_reaction==0)
        reacted_node = NaN;
    else

        reacted_node = unique_node(2<histc(node_id,unique(node_id)));
        
    end


end



