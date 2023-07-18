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

input_bes()


format longEng
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

if flag_plot == 1
    %plot dislocation structure
    set(0,'CurrentFigure',fg1);
    set(gcf, 'Position', [10 400 600 600]);
    plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)
    drawnow
    
    set(0,'CurrentFigure',fg2);
    set(gcf, 'Position', [700 400 600 600]);
    plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)
    
    drawnow
end

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


name_file = strcat('data_',num2str(0));
save(name_file)


%% insert the first loop
xi = 0;
yi = xi;
zi = xi;
shift2 = [xi yi zi 0]-[20*b1 0];
idx2 = randi(size(range_numSegs_list,2),1);

% before inserting loop, save the config incase the loop does not react or
% collide
rn_beforeinsertion = rn;
links_beforeinsertion = links;
[connectivity_beforeinsertion,linksinconnect_beforeinsertion]=genconnectivity(rn,links,maxconnections);

[rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,1,range_numSegs_list,numSegs_list,SurfacePlane,n1list,b1,rn,links,links_start_end_list,fullloop_shift_z,existing_loops,fraction_plim,seglength_multiply);

if seglength_multiply == 1
    num_nodes_perloop = [num_nodes_perloop;6];
    current_loop_numnodes = 6;
else
    num_nodes_perloop = [num_nodes_perloop;numSegs_list(end)*6];
    current_loop_numnodes = numSegs_list(end)*6;
end

loop_time_start = dt_dd;
rn_init = rn;
stuck_indicator = [];

nodes_newloop = [];
for i = 1:size(links_newloop)
    nodes_newloop = [nodes_newloop;links_newloop(i,1);links_newloop(i,1);];
end
nodes_newloop = unique(nodes_newloop);

% add timer for the additional loop for kMC
timecounter_list = [timecounter_list 0];

% add index of current additionalloop into additionalloop_list
loop_list = [loop_list,size(loop_list,2)+1];

[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
consistencycheck(rn,links,connectivity,linksinconnect,rntol);
numlinks=size(links,1);
stackforcevec=zeros(numlinks,3);

% calculate diffusion time for every loop
dt_diffusion_list = [];
for i = 1:size(numSegs_list,1)
    dt_diffusion_i = caldtDiff(T,numSegs_list(i),b_norm,a);
    dt_diffusion_list = [dt_diffusion_list;dt_diffusion_i];
end

% collided_loop contains all loops that have collided and thus won't
% diffuse
stepnum = 0;
collided_loop = [];
absorbed_loop = [];
num_deletedloops = 0;
for curstep=1:totalsteps
    while dt_dd<dt_dd_limit && existing_loops<numprepartial+num_success_loops
        stepnum = stepnum+1;
        disp(strcat('dt_dd=',num2str(dt_dd)))
        idx_loopinsert = floor(dt_dd/dt_loopinsert)+1+numprepartial;
        
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


        %% no reaction happening, continue running until reaction happens, then add new loop in
        flag_newloopreaction = 0; 
        if ~isnan(reacted_node) 
            % reaction happened
            % check if reacted_node contains the node within last inserted loop
            % pre-existing loop reaction does not matter
            for i = 1:size(reacted_node)
                if any(nodes_newloop==reacted_node(i))
                    flag_newloopreaction = 1;
                    break
                end
            end
        end

        flag_newloopabsorption = 0;
        if any(absorbed_loop==  existing_loops)
           
            flag_newloopabsorption = 1;
        end

        flag_newloopcollision = 0;
        if any(collided_loop== existing_loops)
            flag_newloopcollision = 1;
        end

        if flag_newloopreaction==1 || flag_newloopabsorption==1 || flag_newloopcollision == 1 ...
            || flag_loopinsert==1  %idx_loopinsert increases with time, and if time has come for next loop, add it
%           if the last loop did not react nor collide when new loop is
%           inserted, discard the last loop
        
%             if existing_loops+num_deletedloops==idx_loopinsert && flag_newloopreaction==0 && flag_newloopabsorption==0 && flag_newloopcollision == 0
            if flag_loopinsert==1 && flag_newloopreaction==0 && flag_newloopabsorption==0 && flag_newloopcollision == 0
                rn = rn_beforeinsertion;
                links = links_beforeinsertion;
                connectivity = connectivity_beforeinsertion;
                linksinconnect = linksinconnect_beforeinsertion;
                loop_list = loop_list(1:end-1);
                numSegs_list = numSegs_list(1:end-1);
                num_deletedloops = num_deletedloops+1;
                existing_loops = existing_loops-1;
                num_nodes_perloop = num_nodes_perloop(1:end-1);
            end


            name_file = strcat('data_',num2str(existing_loops));
            save(name_file)

            idx2 = randi(size(range_numSegs_list,2),1);
            numSegs2 = range_numSegs_list(idx2);
            
            rn_beforeinsertion = rn;
            links_beforeinsertion = links;
            [connectivity_beforeinsertion,linksinconnect_beforeinsertion]=genconnectivity(rn,links,maxconnections);

            % insert loop  
            if rand(1,1)<0.25 % [11-1]  red
                [rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,1,range_numSegs_list,numSegs_list,SurfacePlane,n1list,b1,rn,links,links_start_end_list,partialloop_initialz(1)+fullloop_shift_z,existing_loops,fraction_plim,seglength_multiply);

            elseif rand(1,1)<0.5 && rand(1,1)>=0.25  % [-1-1-1]  blue
                [rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,1,range_numSegs_list,numSegs_list,SurfacePlane,n2list,b2,rn,links,links_start_end_list,partialloop_initialz(1)+fullloop_shift_z*0.5,existing_loops,fraction_plim,seglength_multiply);
            
            elseif rand(1,1)<0.75 && rand(1,1)>=0.5  % [1-11]  green
                [rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,1,range_numSegs_list,numSegs_list,SurfacePlane,n3list,b3,rn,links,links_start_end_list,partialloop_initialz(1)+fullloop_shift_z,existing_loops,fraction_plim,seglength_multiply);

            elseif rand(1,1)<1 && rand(1,1)>=0.75  % [-111]  yellow
                [rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,1,range_numSegs_list,numSegs_list,SurfacePlane,n4list,b4,rn,links,links_start_end_list,partialloop_initialz(1)+fullloop_shift_z,existing_loops,fraction_plim,seglength_multiply);
            end

            if seglength_multiply == 1
                num_nodes_perloop = [num_nodes_perloop;6];
                current_loop_numnodes = 6;
            else
                num_nodes_perloop = [num_nodes_perloop;numSegs_list(end)*6];
                current_loop_numnodes = numSegs_list(end)*6;
            end
            
            if size(num_nodes_perloop,1) ~= existing_loops
                disp('wrong node number')
            end

            loop_time_start = dt_dd;
            stuck_indicator = [];  % for a new loop, initiate an empty list to see if the loop is stuck

            nodes_newloop = [];
            for i = 1:size(links_newloop)
                nodes_newloop = [nodes_newloop;links_newloop(i,1);links_newloop(i,1);];
            end
            nodes_newloop = unique(nodes_newloop);
            
            dt_diffusion_new = caldtDiff(T,numSegs2,b_norm,a);
            dt_diffusion_list = [dt_diffusion_list;dt_diffusion_new];

            % add timer for the additional loop for kMC
            timecounter_list = [timecounter_list 0];
            % add index of current additional loop into additionalloop_list
            loop_list = [loop_list,size(loop_list,2)+1];
    
            [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
            consistencycheck(rn,links,connectivity,linksinconnect,rntol);
            numlinks=size(links,1);
            stackforcevec=zeros(numlinks,3);
            flag_newloopreaction = 0;
        end
%         end        
    if writeMovie && curstep==1 %Open the movie file for writing
        open(aviobj);
    end

%% use DD to calculate force on each loop

% [~,~,~,fn_diffusion,fseg_diffusion,~]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,totaltime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curstep,printfreq,doshielding,SurfacePlane,dobox,boxD);
rn_allmobile = rn;
rn_allmobile(:,4) = 0;
[~,~,~,fn_diffusion,fseg_diffusion,~]=feval(integrator,rn_allmobile,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,totaltime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curstep,printfreq,doshielding,SurfacePlane,dobox,boxD);
fn_norm = [];
for i = 1:size(fn_diffusion,1)
    fn_norm = [fn_norm;norm(fn_diffusion(i,:))];
end

loopforce_list = zeros(1,size(numSegs_list,1));

sequence_rn = [1:size(rn,1)]';
for ii = 1:size(numSegs_list,1)  % for each loop
    idnodes_currentloop = rn(:,5)==ii;
    sequence_currentloop = sequence_rn(idnodes_currentloop);

    % find all links that connect within the current loop
    idx_links_withinloop = [];
    for jj = 1:size(links,1)
        if any(sequence_currentloop==links(jj,1)) && any(sequence_currentloop==links(jj,2))
            idx_links_withinloop = [idx_links_withinloop;jj];
        end
    end
    b_withinloop = links(idx_links_withinloop,3:5);
    
    for jj = 2:size(b_withinloop,1)
        if isequal(b_withinloop(1,:),b_withinloop(jj,:))
            b = b_withinloop(jj:end,:);
        end
    end
    fn_list = fn_diffusion(idnodes_currentloop,:);
    loopforce_list(ii) = sum(fn_list*b'/norm(b));
%     loopforce_tensor(1:size(fn_list,1),ii) = fn_list*b'/norm(b)
end
%%

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

%% diffusion of nodes 
    for i = 1:size(links,1)
        b_diffuse = links(i,3:5);
        idx1 = links(i,1); % index for node 1
        idx2 = links(i,2); % index for node 2 
%         nodes on the GB do not diffuse
        if rn(idx1,4) == 4
            continue
        end
        % <100> dislocation does not diffuse
        if abs(b_diffuse(1))+abs(b_diffuse(2))+abs(b_diffuse(3)) == b_norm
%             idx_reactedloop = rn(i,end);
            continue
        end
        % collision nodes do not diffuse
        if any(collided_loop==rn(idx1,end))
            continue
        end
        % absorbed loops do not diffuse
        if any(absorbed_loop==rn(idx1,end))
            continue
        end
        % reacted nodes become '0' labeled, and they do not move
        if rn(idx1,end) == 0
            continue
        end
%       apply different diffusion speed for loops 
%       for each loop, consider its net force resolved along burgers vector
%       due to all forces, diffuse along direction of the force
        
        idx_loop = rn(idx1,end);
%       use force on entire loop
%       if three last moves were back and forth, move the loop in the
%       direction of GB
        if size(stuck_indicator,1)<3  % first three jumps, let it jump freely
            if loopforce_list(idx_loop) > 0  % b tau b
                rn(idx1,1:3) = rn(idx1,1:3)+b_diffuse;
                stuck_indicator = [stuck_indicator;1];
            else
                rn(idx1,1:3) = rn(idx1,1:3)-b_diffuse;
                stuck_indicator = [stuck_indicator;-1];
            end
        else
%             when three jumps are in the same direction, let it jump
%             freely
            if stuck_indicator(end) == stuck_indicator(end-1) && stuck_indicator(end)==stuck_indicator(end-2)
                if loopforce_list(idx_loop) > 0  % b tau b
                    rn(idx1,1:3) = rn(idx1,1:3)+b_diffuse;
                    stuck_indicator = [stuck_indicator;1];
                else
                    rn(idx1,1:3) = rn(idx1,1:3)-b_diffuse;
                    stuck_indicator = [stuck_indicator;-1];
                end
%                 when the loop is stuck, force it to jump towards GB
            else
                if b_diffuse(end)>0
                    rn(idx1,1:3) = rn(idx1,1:3)-b_diffuse;
                else
                    rn(idx1,1:3) = rn(idx1,1:3)+b_diffuse;
                end
            end
        end

    end
    
    % subtract timecounter by their diffusion time
    for i = 1:size(timecounter_list,2)
        if timecounter_list(i)>dt_diffusion_list(i)
            timecounter_list(i) = timecounter_list(i)-dt_diffusion_list(i);
        end
    end

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
    
        if flag_plot == 1
            set(0,'CurrentFigure',fg1);
            set(gcf, 'Position', [10 400 600 600]);
            plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)
            drawnow
    
            if writeMovie
                theframe = getframe(gcf);
                framelist1(frame_counter1) = theframe;
                frame_counter1=frame_counter1+1;
            end
            pause(0.01);
    
            set(0,'CurrentFigure',fg2);
            set(gcf, 'Position', [700 400 600 600]);
            plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)
            drawnow
        
            if writeMovie
                theframe = getframe(gcf);
                framelist2(frame_counter2) = theframe;
                frame_counter2=frame_counter2+1;
            end
            pause(0.01);
        end
    
    end

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
%         [rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,node1_index,node2_index]=collision_nodereturn(rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);
        [rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=collision_test(rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);
    end
    
    if(doremesh)
        %remesh
        [rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=remesh(rnnew_collide,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);      
    end

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
   
    
    % count number of nodes for each loop, if its less than 6, its collided
    for i = 1:existing_loops
        num = sum(rn(:,end)==i);
        if num<current_loop_numnodes
            collided_loop = [collided_loop;i];
        end
    end
    collided_loop = unique(collided_loop);
    
    %%  relax when collision happens
    if any(collided_loop==existing_loops) && 1==0
        disp('relaxing')
        for i = 1:20
        [rn,vn,dt,fn,fseg,totaltime]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,totaltime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curstep,printfreq,doshielding,SurfacePlane,dobox,boxD);
        % there is force due to surface plane in segforcevec.m
        dt_dd = dt_dd+dt;
        timecounter_list = timecounter_list+dt;
        set(0,'CurrentFigure',fg1);
        set(gcf, 'Position', [10 400 600 600]);
        plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)
        drawnow
        if writeMovie
            theframe = getframe(gcf);
            framelist1(frame_counter1) = theframe;
            frame_counter1=frame_counter1+1;
        end
        pause(0.01);
        set(0,'CurrentFigure',fg2);
        set(gcf, 'Position', [700 400 600 600]);
        plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)
        drawnow
        if writeMovie
            theframe = getframe(gcf);
            framelist2(frame_counter2) = theframe;
            frame_counter2=frame_counter2+1;
        end
        pause(0.01);
        end
    end

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

%% output video
if writeMovie %Record a movie if specified to do so   
    newVid = VideoWriter('video1', 'MPEG-4'); % New
    newVid.FrameRate = 10;
    newVid.Quality = 100;
    open(newVid);
    for frame = 1:length(framelist1)
        writeVideo(newVid,framelist1(frame));
    end
    close(newVid);
end
if writeMovie %Record a movie if specified to do so   
    newVid = VideoWriter('video2', 'MPEG-4'); % New
    newVid.FrameRate = 10;
    newVid.Quality = 100;
    open(newVid);
    for frame = 1:length(framelist2)
        writeVideo(newVid,framelist2(frame));
    end
    close(newVid);
end

if flag_plot == 1
    set(0,'CurrentFigure',fg1);
    set(gcf, 'Position', [10 400 600 600]);
    plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)

    set(0,'CurrentFigure',fg2);
    set(gcf, 'Position', [700 400 600 600]);
    plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)
end

save relaxround.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%