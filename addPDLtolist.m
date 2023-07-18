function [rn,links,links0] = addPDLtolist(rn,links,rn0,links0)

    links0(:,1:2) = links0(:,1:2)+size(rn,1);
    links = [links;links0;];
    rn = [rn;rn0];

end