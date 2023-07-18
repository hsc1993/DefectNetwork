% This DEMO calculates a Voronoi diagram with arbitrary points in arbitrary
% polytope/polyheron in 2D/3D

% 20	50	200	500	1000	2000
clear all;close all;clc
%% generate random samples
n = 10;        % number of points
m = 0;         % number of boundary point-candidates
d = 3;          % dimension of the space
tol = 1e-07;            % tolerance value used in "inhull.m" (larger value high precision, possible numerical error)
% pos0 = rand(n,d);       % generate random points
% bnd0 = rand(m,d)       % generate boundary point-candidates
% % bnd0 = [1 0 0; 0 0 1;0 1 0; 1 1 0; 0 1 1; 1 0 1; 1 1 1; 0 0 0];       % generate boundary point-candidates
num_bndpoints_inclusion = 8;
case_box = 1;
% case_box = 2;

% 
if case_box == 1
    xmax = 500;
    xmin = 100;
    ymax = 500;
    ymin = 100;
    zmax = 500;
    zmin = 100;
end

bnd0 = [0.5*(xmin+xmax) 0.5*(ymin+ymax) 0.5*(zmin+zmax); xmin ymin zmax; xmax ymin zmin;  xmin ymax zmin; xmax ymax zmin; xmin ymax zmax; xmax ymin zmax; xmax ymax zmax; xmin ymin zmin];       % generate boundary point-candidates

rng('shuffle')
xplane_bndpnts = 10;
yplane_bndpnts = 10;
zplane_bndpnts = 10;

% x plane
for i = 1:xplane_bndpnts
    yrand1 = (ymax-ymin)*rand(1) + ymin;
    yrand2 = (ymax-ymin)*rand(1) + ymin;
    
    zrand1 = (zmax-zmin)*rand(1) + zmin;
    zrand2 = (zmax-zmin)*rand(1) + zmin;    
    bnd0 = [bnd0;xmin yrand1 zrand1];
    bnd0 = [bnd0;xmax yrand2 zrand2];
end



% y plane
for i = 1:yplane_bndpnts
    xrand1 = (xmax-xmin)*rand(1) + xmin;
    xrand2 = (xmax-xmin)*rand(1) + xmin;
    
    zrand1 = (zmax-zmin)*rand(1) + zmin;
    zrand2 = (zmax-zmin)*rand(1) + zmin;    
    bnd0 = [bnd0;xrand1 ymin zrand1];
    bnd0 = [bnd0;xrand2 ymax zrand2];
end

% z plane
for i = 1:zplane_bndpnts
    xrand1 = (xmax-xmin)*rand(1) + xmin;
    xrand2 = (xmax-xmin)*rand(1) + xmin;   
    
    yrand1 = (ymax-ymin)*rand(1) + ymin;
    yrand2 = (ymax-ymin)*rand(1) + ymin;
    
    bnd0 = [bnd0;xrand1 yrand1 zmin];
    bnd0 = [bnd0;xrand2 yrand2 zmax];
end


K = convhull(bnd0);
bnd_pnts = bnd0(K,:);   % take boundary points from vertices of convex polytope formed with the boundary point-candidates
% bnd_pnts = unique(bnd_pnts,'rows');

hold on
scatter3(bnd0(:,1),bnd0(:,2),bnd0(:,3))
scatter3(bnd_pnts(:,1),bnd_pnts(:,2),bnd_pnts(:,3),'*')
hold off


pos0 = [];
for i=1:n
    pos_x = randi([xmin xmax],1,1);
    pos_y = randi([ymin ymax],1,1);
    pos_z = randi([zmin zmax],1,1);
    pos01_temp = [pos_x,pos_y,pos_z];
    pos0 = [pos0;pos01_temp];
end


% indice_bndpoints_inclusion = randperm(8,num_bndpoints_inclusion);
% 
% K = convhull(bnd0);
% bnd_pnts = bnd0(K,:);   % take boundary points from vertices of convex polytope formed with the boundary point-candidates

%% take points that are in the boundary convex polytope
in = inhull(pos0,bnd0,[],tol); 
% inhull.m is written by John D'Errico that efficiently check if points are
% inside a convex hull in n dimensions
% We use the function to choose points that are inside the defined boundary
u1 = 0;
for i = 1:size(pos0,1)
    if in(i) ==1
        u1 = u1 + 1;
        pos(u1,:) = pos0(i,:);
    end
end


%% 
% =========================================================================
% INPUTS:
% pos       points that are in the boundary      n x d matrix (n: number of points d: dimension) 
% bnd_pnts  points that defines the boundary     m x d matrix (m: number of vertices for the convex polytope
% boundary d: dimension)
% -------------------------------------------------------------------------
% OUTPUTS:
% vornb     Voronoi neighbors for each generator point:     n x 1 cells
% vorvx     Voronoi vertices for each generator point:      n x 1 cells
% =========================================================================

[vornb,vorvx] = polybnd_voronoi(pos,bnd_pnts);


%% PLOT

for i = 1:size(vorvx,2)
    col(i,:)= rand(1,3);
end


switch d
    case 2
        figure('position',[0 0 600 600],'Color',[1 1 1]);
        for i = 1:size(pos,1)
        plot(vorvx{i}(:,1),vorvx{i}(:,2),'-r')
        hold on;
        end
        plot(bnd_pnts(:,1),bnd_pnts(:,2),'-');
        hold on;
        plot(pos(:,1),pos(:,2),'Marker','o','MarkerFaceColor',[0 .75 .75],'MarkerEdgeColor','k','LineStyle','none');
        axis('equal')
        axis([0 1 0 1]);
        set(gca,'xtick',[0 1]);
        set(gca,'ytick',[0 1]);        
    case 3
        figure('position',[0 0 600 600],'Color',[1 1 1]);
        for i = 1:size(pos,1)
            flag_outside = 0;
            
            if max(vorvx{i}(:,1))>xmax+0.1
                flag_outside = 1;
            end
            if min(vorvx{i}(:,1))<xmin-0.1
                flag_outside = 1;
            end
            if max(vorvx{i}(:,2))>ymax+0.1
                flag_outside = 1;
            end
            if min(vorvx{i}(:,2))<ymin-0.1
                flag_outside = 1;
            end
            if max(vorvx{i}(:,3))>zmax+0.1
                flag_outside = 1;
            end
            if min(vorvx{i}(:,3))<zmin-0.1
                flag_outside = 1;
            end
            
            if flag_outside == 1
                continue
            end



        K = convhulln(vorvx{i});
%         scatter3(vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');

        trisurf(K,vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'FaceColor',col(i,:),'FaceAlpha',0.5,'EdgeAlpha',0)
%         nb = cell2mat(vornb(1))
%         nb_idx = nb(1)
%         trisurf(K,vorvx{1}(:,1),vorvx{1}(:,2),vorvx{1}(:,3),'FaceColor','r','FaceAlpha',0.1,'EdgeAlpha',1)
%         trisurf(K,vorvx{nb_idx}(:,1),vorvx{nb_idx}(:,2),vorvx{nb_idx}(:,3),'FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1)
        hold on;
        end
%         scatter3(pos(:,1),pos(:,2),pos(:,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
%         scatter3(pos(1,1),pos(1,2),pos(1,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
%         scatter3(pos(nb_idx,1),pos(nb_idx,2),pos(nb_idx,3),'Marker','x','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
        axis('equal')
%         axis([0 1 0 1 0 1]);
%         set(gca,'xtick',[0 1]);
%         set(gca,'ytick',[0 1]);
%         set(gca,'ztick',[0 1]);
%         set(gca,'FontSize',16);
        xlabel('X');ylabel('Y');zlabel('Z');
        
        axis off
        
% Convert cell to a table and use first row as variable names
size(vorvx,2)
vorvx_namestring = ['n',num2str(n),'_vorvx.txt'];
vorvx_fileID = fopen(vorvx_namestring,'a');

for i = 1:size(vorvx,2)
    M = cell2mat(vorvx(i));
    A1 = M(:,1);
    A2 = M(:,2);
    A3 = M(:,3);
    
    fprintf(vorvx_fileID, '[');
    for j = 1:size(A1)
        formatSpec1 = '[%8.4f,';
        fprintf(vorvx_fileID,formatSpec1,A1(j));

        formatSpec2 = '%8.4f,';
        fprintf(vorvx_fileID,formatSpec2,A2(j));

        formatSpec3 = '%8.4f],\n';
        fprintf(vorvx_fileID,formatSpec3,A3(j));
    end
    fprintf(vorvx_fileID, ']');
    fprintf(vorvx_fileID, ',');
end
fclose(vorvx_fileID);


vornb_namestring = ['n',num2str(n),'_vornb.txt'];
vornb_fileID = fopen(vornb_namestring,'a');

for i = 1:size(vornb,2)
    M = cell2mat(vornb(i));
    A = M;
    
    fprintf(vornb_fileID, '[');
    for j = 1:size(A,2)
        formatSpec = '%i,';
        fprintf(vornb_fileID,formatSpec,A(j));        
    end
    fprintf(vornb_fileID, ']');
    fprintf(vornb_fileID, ',');

end
fclose(vornb_fileID);

end
