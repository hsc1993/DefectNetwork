function [rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,numloops,range_numSegs_list,numSegs_list,SurfacePlane,nlist,b,rn,links,links_start_end_list,shift_z,existing_loops,fraction_plim,seglength_multiply)

% use 1.5 because full range is 2*plim, and we only put full loop in
% 1.5*plim range to ensure full loops do not get wasted by diffusing out


xi_list1 = randperm(2*fraction_plim*plim,numloops)-fraction_plim*plim;
yi_list1 = randperm(2*fraction_plim*plim,numloops)-fraction_plim*plim;

xy_list1 = [];
for i = 1:size(xi_list1,2)
    xy_list1 = [xy_list1;xi_list1(i) yi_list1(i)];
end

shift = [];


for i = 1:numloops
    idx1 = randi(size(range_numSegs_list,2),1);
    current_numSegs = range_numSegs_list(idx1);
    numSegs_list = [numSegs_list;current_numSegs];

    xi = xy_list1(i,1);
    yi = xy_list1(i,2);
    zi = -SurfacePlane(2)/SurfacePlane(3)*yi-SurfacePlane(4)/SurfacePlane(3)+SurfacePlane(4)-5*rand(1,1)+shift_z;
    shift1i = [xi yi zi];
    shift = [shift;shift1i];


    existing_loops = existing_loops+1;
    [rn,links,links_newloop,links_start_end_list] = addnewloop(nlist,b,current_numSegs,shift1i,rn,links,links_start_end_list,existing_loops,seglength_multiply);
    
    
end

end


