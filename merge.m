clear;clc;close all;
format short

f.Name = 'Final Configuration';
f.Position = [100 100 1000 800];


fileID1 = fopen(strcat('dislocation.txt'),'w');
fileID2 = fopen(strcat('dislocation_burgers.txt'),'w');

rn_all = [];
for idx = 1:9
    load(strcat('relaxround_',num2str(idx),'.mat'))


    for idx_link = 1:size(links,1)
            
        idx_rn1 = links(idx_link,1);
        idx_rn2 = links(idx_link,2);
        rn1 = rn(idx_rn1,1:3);
        rn2 = rn(idx_rn2,1:3);

        if idx == 2
            rn1 = rn1+[2*plim,0,0];
            rn2 = rn2+[2*plim,0,0];
        elseif idx == 3
            rn1 = rn1+[4*plim,0,0];
            rn2 = rn2+[4*plim,0,0];
        elseif idx == 4
            rn1 = rn1+[0,2*plim,plim];
            rn2 = rn2+[0,2*plim,plim];
        elseif idx == 5
            rn1 = rn1+[2*plim,2*plim,plim];
            rn2 = rn2+[2*plim,2*plim,plim];
        elseif idx == 6
            rn1 = rn1+[4*plim,2*plim,plim];
            rn2 = rn2+[4*plim,2*plim,plim];
        elseif idx == 7
            rn1 = rn1+[0,4*plim,2*plim];
            rn2 = rn2+[0,4*plim,2*plim];
        elseif idx == 8
            rn1 = rn1+[2*plim,4*plim,2*plim];
            rn2 = rn2+[2*plim,4*plim,2*plim];
        elseif idx == 9
            rn1 = rn1+[4*plim,4*plim,2*plim];
            rn2 = rn2+[4*plim,4*plim,2*plim];
        end
        rn_all = [rn_all;rn1];
        rn_all = [rn_all;rn2];

        b = links(idx_link,3:5);
        fprintf(fileID1,'[[%6.6f,%6.6f,%6.6f],[%6.6f,%6.6f,%6.6f]],',rn1,rn2);
        fprintf(fileID2,'[%6.6f,%6.6f,%6.6f],',b);

    end


end


rn_all_xmin = min(rn_all(:,1))
rn_all_xmax = max(rn_all(:,1))
rn_all_ymin = min(rn_all(:,2))
rn_all_ymax = max(rn_all(:,2))
rn_all_zmin = min(rn_all(:,3))
rn_all_zmax = max(rn_all(:,3))

fclose(fileID1);


f.Name = 'Initial Configuration';
f.Position = [100 100 1000 800];

fg1 = figure(1);
fg2 = figure(2);
set(0,'CurrentFigure',fg1);
set(gcf, 'Position', [10 400 600 600]);
plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane)
% plotnodes(rn,links_init0,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)

set(0,'CurrentFigure',fg2);
set(gcf, 'Position', [700 400 600 600]);
% plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)
plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane)








