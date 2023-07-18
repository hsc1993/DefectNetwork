%Dislocation Dynamics simulation in 3-dimension
%
%the best way to run it is by "run dd3d"
%
%Features:
%infinite boundary condition (no pbc)
%linear mobility law (mobfcc0,mobfcc1)
%no remesh, no collision detection, no reaction
%N^2 interaction (no neighbor list, no fast multipole)

%Data structure:
%NMAX:    maximum number of nodes (including disabled ones)
%LINKMAX: maximum number of links (including disabled ones)
%rn: (NMAX,4) array of nodal positions (last column is flag: -1 means disabled)
%vn: (NMAX,3) array of nodal velocities
%fn: (NMAX,3) array of nodal forces
%links: (LINKMAX,8) array of links (id1,id2,bx,by,bz,nx,ny,nz)

clc;close all;clear;
load data_61
doremesh = 1;
rn_beforerelax = rn;
integrator='int_eulerbackward_relax_noprismatic';


dt_dd_limit=5.5e-16;

format short

mex -O SegSegForces.c
mex -O SegSegForcesVector.c

if ~exist('restrictSurfaceNodes') %If ~exist, create a bool that allows surface ndoes
    restrictSurfaceNodes=false;
end

if ~exist('writeMovie') %Create the writeMovie variable if not predefined
    writeMovie = false;
end

if writeMovie %Record a movie if specified to do so   
    moviename='video.avi'; %REMEMBER TO CHANGE THE NAME OF THE FILE
    aviobj=VideoWriter(moviename);

    aviobj.Quality=100;
    %aviobj.Compression='Indeo5';
    %open(aviobj);
    
    frame =1;
end

%default value if run by itself (e.g. not through "rundd3d")
if(~exist('rn'))
    initparams;  %default parameter settings (can be override by restart files below)
    load frsource; %read in dislocation configuration and setting files
end
if(~exist('dt'))
    dt=dt0;
end
if(~exist('doshielding'))
    doshielding=0;
end
if(~exist('dopartials'))
    dopartials=0;
    stackforcevec=[0  0  0;
                   0  0  0;
                   0  0  0;
                   0  0  0];
end
if(~exist('doSFT'))
    doSFT=0;
end

if(~exist('dobox'))
    dobox=0;
    boxD = [-202.2 202.2; -202.2 202.2; -202.2 202.2];
end
if(~exist('SFT_plane') )
    SFT_plane=[];

    
else
    SFTtol=1e-5;
    for i=0:doSFT-1
        SFT_index=4*i+1;
        SFT_plane(SFT_index,4:6)=SFT_plane(SFT_index,4:6)-SFTtol.*SFT_plane(SFT_index,1:3);
        SFT_plane(SFT_index+1,4:6)=SFT_plane(SFT_index+1,4:6)-SFTtol.*SFT_plane(SFT_index+1,1:3);
        SFT_plane(SFT_index+2,4:6)=SFT_plane(SFT_index+2,4:6)-SFTtol.*SFT_plane(SFT_index+2,1:3);
        SFT_plane(SFT_index+3,4:6)=SFT_plane(SFT_index+3,4:6)-SFTtol.*SFT_plane(SFT_index+3,1:3);
    end
end
if(~exist('doinclusion'))
    doinclusion=0;
    inclusion_pos_rad=0;
end

% cleanup the empty node and link entries at the end of the initial data structures
[rn,links]=cleanupnodes(rn,links);

%add the partials to the code. A flag is needed to have the right sense of the stacking fault force(Enrique Sep 2005)


if((findstr(mobility,'fcc'))&(dopartials))
    numlinks=size(links,1);
    if(~exist('stackforcevec'))
        stackforcevec=zeros(numlinks,3);
    end
    %[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
    %consistencycheck(rn,links,connectivity,linksinconnect,rntol);
    %[rn,links,stackforcevec]=stackfaultsplit(rn,links,connectivity,rann,gamma,stackforcevec,maxconnections);
end

% cleanup the empty node and link entries at the end of the initial data structures
[rn,links]=cleanupnodes(rn,links);

% genererate the connectivity list from the list of links
[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
consistencycheck(rn,links,connectivity,linksinconnect,rntol);

%plot dislocation structure
% set(0,'CurrentFigure',fg1);
% set(gcf, 'Position', [10 400 600 600]);
% plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)
% drawnow
% 
% set(0,'CurrentFigure',fg2);
% set(gcf, 'Position', [700 400 600 600]);
% plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)
% 
% drawnow

totaltime=0;
dt=min(dt,dt0);
mdold=10;
counter=0;
ppt=0;

frame_counter1=1;
frame_counter2=1;

dt_dd = 0;
timecounter_list = zeros(1,numprepartial);
reacted_list = [];
loop_list = [1:numprepartial];
% nextloop_idx = numprepartial+1;

% collided_loop contains all loops that have collided and thus won't
% diffuse
collided_loop = [];
absorbed_loop = [];
num_deletedloops = 0;

relax_round = 0;
for curstep=1:totalsteps
%     while dt_dd<dt_dd_limit && relax_round<400
    while relax_round<800
        disp(strcat('dt_dd=',num2str(dt_dd)))
        idx_loopinsert = floor(dt_dd/dt_loopinsert)+1+numprepartial;
        
        relax_round = relax_round+1;
        flag_loopinsert = 0;
        if (dt_dd-loop_time_start)>dt_loopinsert
            flag_loopinsert = 1;
        end
        % use nodes_newloop and 'links_old vs links' to determine if reaction
        % has taken place for the new loop

        nodes_w_3links = detect_reaction(links);
        reacted_node = [nodes_w_3links];
        for i = 1:size(links,1)
            for j = 1:size(nodes_w_3links,1)
                if links(i,1) == nodes_w_3links(j)
                    reacted_node = [reacted_node;links(i,2)];
                end
                if links(i,2) == nodes_w_3links(j)
                    reacted_node = [reacted_node;links(i,1)];
                end
            end
        end
        reacted_node = unique(reacted_node);

    if writeMovie && curstep==1 %Open the movie file for writing
        open(aviobj);
    end

%% check if [100] segment has been created
% note that some reaction does not generate [100] but a 'collision' does
% happen
idx_reactedloop = [];
for i = 1:size(links,1)
    b_diffuse = links(i,3:5);

    if abs(b_diffuse(1))+abs(b_diffuse(2))+abs(b_diffuse(3)) == b_norm
        idx_reactedloop = [idx_reactedloop;rn(i,end)];
        continue
    end 
end

%% check if loop absorption has taken place
% detect absorbed loop
absorption_ratio = 0.1;
num_loops = size(loop_list,2);
for idx_loop = 1:num_loops
    rn_absorbedloop = rn(rn(:,5)==idx_loop,:); % get the rn for the current absorbed loop
    num_absorbednodes = sum(rn_absorbedloop(:,4)==4);
    num_allnodes = size(rn_absorbedloop,1);
    current_absorption_ratio = num_absorbednodes/num_allnodes;
    if current_absorption_ratio > absorption_ratio
        absorbed_loop = [absorbed_loop;idx_loop];
    end
end
absorbed_loop = unique(absorbed_loop);


    %% Check if nodes reaches a surface
    [lrn,lrn2]=size(rn);
    rnnew_checksurf=rn;
    for i=1:lrn
       if dot(SurfacePlane,[rnnew_checksurf(i,1:3) 1])< 0
           rnnew_checksurf(i,4) = 4; % change the forth column to represent attaching to surface
       end
    end

    %% DD
    %integrating equation of motion
    rn = rnnew_checksurf;
    
    for i = 1:1
        [rn,vn,dt,fn,fseg,totaltime]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,totaltime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curstep,printfreq,doshielding,SurfacePlane,dobox,boxD);
        % there is force due to surface plane in segforcevec.m
        dt_dd = dt_dd+dt;
        timecounter_list = timecounter_list+dt;
        
%         set(0,'CurrentFigure',fg1);
%         set(gcf, 'Position', [10 400 600 600]);
%         plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)
%         drawnow
%    
%         if writeMovie
%             theframe = getframe(gcf);
%             framelist1(frame_counter1) = theframe;
%             frame_counter1=frame_counter1+1;
%         end
%         pause(0.01);
%     
%         set(0,'CurrentFigure',fg2);
%         set(gcf, 'Position', [700 400 600 600]);
%         plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)
%         drawnow
%     
%         if writeMovie
%             theframe = getframe(gcf);
%             framelist2(frame_counter2) = theframe;
%             frame_counter2=frame_counter2+1;
%         end
%         pause(0.01);
    
    end

    name_file = strcat('relaxround_',num2str(relax_round));
    save(name_file)

    rnnew_collide=[rn(:,1:3) vn rn(:,4) rn(:,end)];
    rnnew_genbelong=[rn(:,1:3) vn rn(:,4) rn(:,5)];
    
    linksnew=links;
    connectivitynew=connectivity;
    linksinconnectnew=linksinconnect;
    fsegnew=fseg;
    stackforcevecnew=stackforcevec;
    if(docrossslip)
        %Fleischer and Escaig Cross-slip
        crossslip=0;
        [rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,crossslip,printnode]=crossslipfuncnewFrank(rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,crossslip,integrator,dt,dt0,rmax,totaltime,areamin,lmin,printnode,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);
%         view(viewangle);
%         drawnow
        %pause;
    end
    if(doseparation)
        %spliting of nodes with 4 or more connections
        [rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=separationnew(rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,mobility,MU,NU,a,Ec,rann,appliedstress,dopartials,stackforcevecnew,rann,gamma,rntol,maxconnections,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);
    end
    if(docollision)
        %collision detection and handling
%         [rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=collision_test(rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);
        [rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,node1_index,node2_index]=collision_nodereturn(rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);
    end
    
    if(doremesh)
        %remesh
%         size(rnnew_collide)
        [rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=remesh(rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);      
        size(rnnew_collide)
    end

    if ~isnan(node1_index) && ~isnan(node2_index)
        collided_loop = [collided_loop;node1_index];
        collided_loop = [collided_loop;node2_index];
    end
    collided_loop = unique(collided_loop);
%     make sure that both collided loops are added to the 'collided_loop'

    rn_belonglist = [];
    for i = 1:size(rnnew_collide,1)
        success = 0;
        for j = 1:size(rnnew_genbelong,1)
            if isequal(rnnew_collide(i,:),rnnew_genbelong(j,1:end-1))
                rn_belonglist = [rn_belonglist;rnnew_genbelong(j,end)];
                success = 1;
                continue
            end
        end
        if success == 0
%             if no matching is found, assign '0' to the belonglist,
%             indicating the node no longer belongs to a single loop
            rn_belonglist = [rn_belonglist;0];
        end
    end

    rn=[rnnew_collide(:,1:3) rnnew_collide(:,7) rn_belonglist];
%     rn=[rnnew_collide(:,1:3) rnnew_collide(:,7) rnnew_genbelong(:,end)];
    vn=rnnew_collide(:,4:6);

    links=linksnew;
    connectivity=connectivitynew;
    linksinconnect=linksinconnectnew;
    fseg=fsegnew;
    stackforcevec=stackforcevecnew;
    %store run time information
    %time step
    data(curstep,1)=dt;
    %save restart
    consistencycheck(rn,links,connectivity,linksinconnect,rntol);
%     conservstakforce(stackforcevec,links,rn,connectivity);
   
    end % end of dd
end

% %% output video
% if writeMovie %Record a movie if specified to do so   
%     newVid = VideoWriter('video1', 'MPEG-4'); % New
%     newVid.FrameRate = 10;
%     newVid.Quality = 100;
%     open(newVid);
%     for frame = 1:length(framelist1)
%         writeVideo(newVid,framelist1(frame));
%     end
%     close(newVid);
% end
% if writeMovie %Record a movie if specified to do so   
%     newVid = VideoWriter('video2', 'MPEG-4'); % New
%     newVid.FrameRate = 10;
%     newVid.Quality = 100;
%     open(newVid);
%     for frame = 1:length(framelist2)
%         writeVideo(newVid,framelist2(frame));
%     end
%     close(newVid);
% end

% set(0,'CurrentFigure',fg1);
% set(gcf, 'Position', [10 400 600 600]);
% plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)
% 
% set(0,'CurrentFigure',fg2);
% set(gcf, 'Position', [700 400 600 600]);
% plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)

save finalconfig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%