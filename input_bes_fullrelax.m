%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input file for dislocation - precipitate interaction   %  
% Eshelby inclustions                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;

format short
% From Ninive et al: 
% - pure Al GGA-PBE
% mu = 4.5628e-09 N/A2*a2, E = 6.38793e-09 N/A2*a2, B = -1.16291e+08 (N/A2*a2)^-1
writeMovie = true; %Bool to save a movie or not
restrictSurfaceNodes = false; %Bool to keep nodes off of the particles surface or not
sim_params = simulation_parameters; %Define all the simulation parameters needed for emission/vacancies

%% simulation constants
a0 = 2.86;
b_globalframe = a0/2*[1 1 1];
Va = a0^3/2;   % BCC atomic volume  A^3
w = a0*sqrt(6)/3;  %  1/6[112]
b_norm = a0/2*norm([1 1 1]);
kb = 8.617333262145e-5;	% eV * K^-1
gamma=6.612e-11; %Lu et al 1999    related to SF
MU = 4.5628e-09; % Ninive 2014    
MU = 94*1e-11; % shear modulus   94 GPa=94e-11 N/A^2
NU = 0.33;% Ninive 2014     poisson ratio
a=1/sqrt(2);%b %1.75*(1/sqrt(2));% 3b  lmin/sqrt(3)*0.5;

Ec = [MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1)];% [Perfect Shockley StairRod Hirth Frank] [J/l3]MU/(4*pi)*log(a/0.1);
Ri=5; % Radius inclusion - to be varied 

%% dislocation line configuration
numSeg_smallest = 5;
numSeg_biggest = 20;
numSeg_biggest = 7;

range_numSegs_list = [numSeg_smallest:numSeg_biggest];

numSegs = 5; % 1->4.6A, 5-> 23A, 10->46A, 12->55A 
Diameter = 2*(sqrt(3)+2)/4*numSegs*b_norm;
% 20b diameter 25b
numloops2 = 0;

segmentlength = sqrt(6)*a0/3;

plim = 80;  % plot limit
% SurfacePlane=[0 0 0 0];


boxlength = 100;

dt0=3e-11; % 1ps = 1e-12s
dt=3e-10;  % this dt was used in original DD, but when KMC is introduced, dt is replaced with poisson variate
stress0 = [0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]*1e-11; %in units of N/A^2, before exponential is in units of GPa
stress = 0*1e-11; %in units of N/A^2, before exponential is in units of GPa
appliedstress = [stress 0 0;  % shear   shear
                0 stress 0
                0 0 stress];

%% prismatic loop parameters
T = 2000;  % k
totalsteps = 1;

rn = [];
links = [];
links_start_end_list = [];

%% simulation time
dt_dd_limit = 10e-11;
dt_loopinsert = 5e-12;
% dt_loopinsert = 1e-11;

%% two sets of burgers vectors can produce [100] type:
%  [11-1] [-1-1-1] [1-11] [-111]  we choose this set
%  [111] [1-1-1] [-11-1] [-1-11]

%% loop 1 
b1 = 1/sqrt(3)*[1 1 -1]*b_norm;
n11 = 1/sqrt(2)*[ 0 1 1];
n12 = 1/sqrt(2)*[1 0 1];
n13 = 1/sqrt(2)*[1 -1 0];
n14 = 1/sqrt(2)*[0 -1 -1];
n15 = 1/sqrt(2)*[-1 0 -1];
n16 = 1/sqrt(2)*[-1 1 0];
n1list = [n11;n12;n13;n14;n15;n16];

%% loop 2 
b2 = 1/sqrt(3)*[-1 -1 -1 ]*b_norm;
n21 = 1/sqrt(2)*[ 0 -1 1];
n22 = 1/sqrt(2)*[1 -1 0];
n23 = 1/sqrt(2)*[1 0 -1];
n24 = 1/sqrt(2)*[0 1 -1];
n25 = 1/sqrt(2)*[-1 1 0];
n26 = 1/sqrt(2)*[-1 0 1];
n2list = [n21;n22;n23;n24;n25;n26];

%% loop 3
b3 = 1/sqrt(3)*[1 -1 1]*b_norm;
n31 = 1/sqrt(2)*[ 1 1 0];
n32 = 1/sqrt(2)*[0 1 1];
n33 = 1/sqrt(2)*[-1 0 1];
n34 = 1/sqrt(2)*[-1 -1 0];
n35 = 1/sqrt(2)*[0 -1 -1];
n36 = 1/sqrt(2)*[1 0 -1];
n3list = [n31;n32;n33;n34;n35;n36];

%% loop 4
b4 = 1/sqrt(3)*[-1 1 1]*b_norm;
n41 = 1/sqrt(2)*[1 1 0];
n42 = 1/sqrt(2)*[0 1 -1];
n43 = 1/sqrt(2)*[-1 0 -1];
n44 = 1/sqrt(2)*[-1 -1 0];
n45 = 1/sqrt(2)*[0 -1 1];
n46 = 1/sqrt(2)*[1 0 1];
n4list = [n41;n42;n43;n44;n45;n46];



SurfacePlane=[0 -1 2 0];
% numSegs_list stores num of segments for each loop
numSegs_list = [];

existing_loops = 0;

numloops1 = 5;
numloops2 = numloops1;
numloops3 = numloops1;
numloops4 = numloops1;

partialloop_initialz = [0 -15 -10 0];
fullloop_shift_z = 40;

% [11-1]  red
[rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,numloops1,range_numSegs_list,numSegs_list,SurfacePlane,n1list,b1,rn,links,links_start_end_list,partialloop_initialz(1),existing_loops,0.7);

% [-1-1-1]  blue
[rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,numloops2,range_numSegs_list,numSegs_list,SurfacePlane,n2list,b2,rn,links,links_start_end_list,partialloop_initialz(2),existing_loops,0.7);

% [1-11]  green
[rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,numloops3,range_numSegs_list,numSegs_list,SurfacePlane,n3list,b3,rn,links,links_start_end_list,partialloop_initialz(3),existing_loops,0.7);

% [-111]  yellow
[rn,links,links_newloop,links_start_end_list,numSegs_list,existing_loops] = addnewloop_main(plim,numloops4,range_numSegs_list,numSegs_list,SurfacePlane,n4list,b4,rn,links,links_start_end_list,partialloop_initialz(4),existing_loops,0.7);

% fraction of area for full loop introduction
fraction_plim=0.65;


numprepartial = numloops1+numloops2+numloops3+numloops4;

links_init0 = links;
rn_init0 = rn;


[lrn,lrn2]=size(rn);
rnnew_checksurf=rn;
for i=1:lrn
   if dot(SurfacePlane,[rnnew_checksurf(i,1:3) 1]) < 0
       rnnew_checksurf(i,4) = 4; % change the forth column to represent attaching to surface
   end
end
rn = rnnew_checksurf;

% select only segments above GB plane to evaluate force

% find first and last node of each loop
idxlist_firstlast = [];
for i = 1:size(rn,1)-1
    if rn(i,5)~=rn(i+1,5)
        idxlist_firstlast = [idxlist_firstlast;i;i+1];
    end
end
% rn(idxlist_firstlast,4) = 4;
% generate selection for rn's with two endpoints being '4'
endpoints_partialloop = [];
for i = 1:size(rn,1)-1
    % for the first link of the current loop, if its '0', make it '4'
    if rn(i,4) == 0 && rn(i+1,4) == 4 
        endpoints_partialloop = [endpoints_partialloop;i+1];
    end
    if rn(i,4) == 4 && rn(i+1,4) == 0
        endpoints_partialloop = [endpoints_partialloop;i];
    end
end

if rn(1,4) == 0
    endpoints_partialloop = [1;endpoints_partialloop];
end
if rn(end,4) == 0
    endpoints_partialloop = [endpoints_partialloop;size(rn,1)];
end


% % use endpoints generate links for parital loop segments (notice this process will also link different loop, need to delete those later)
% idx_partialloop_rn = [];
% idx_partialloop_links = [];
% links_id12_list = [];
% for i = 1:size(endpoints_partialloop,1)-1
%     if mod(i,2) == 1
%         idx_list_onepartial = [endpoints_partialloop(i):endpoints_partialloop(i+1)]';
%         idx_partialloop_rn = [idx_partialloop_rn;idx_list_onepartial];
%         idx_partialloop_links = [idx_partialloop_links;idx_list_onepartial(1:end-1)];
% %         
%         current_linksize = size(links_id12_list,1);
%         if size(links_id12_list,1) == 0
%             links_id12_list_start = -1;
%         else
%             links_id12_list_start = links_id12_list(end,1);
%         end
% 
%         for j = 1:size(idx_list_onepartial,1)-1
% %             links_id12_list = [links_id12_list;[current_linksize+j+floor(i/2) current_linksize+j+1+floor(i/2)]];
% %             links_id12_list = [links_id12_list;[current_linksize+j+i-1 current_linksize+j+1+i-1]];
%             links_id12_list = [links_id12_list;[links_id12_list_start+j+1 links_id12_list_start+j+2]];
%         end
%     end
% end
% rn_partialloop = rn(idx_partialloop_rn,:);
% links_partialloop = links(idx_partialloop_links,:);
% % links_partialloop(:,1:2) = links_id12_list;
% 
% rn_init1 = rn_partialloop;
% links_init1 = links_partialloop;

% for i = 1:size(links_partialloop,1)
%     links_partialloop(i,1) = i;
%     links_partialloop(i,2) = i+1;
% end
% links_partialloop = links_partialloop(1:end-1,:);


% % correct the links that connect different loops
% idx_links_partialloop_removeinterloop = [];
% idx_links_partialloop_removeinterloop_inverse = [];
% for i = 1:size(links_partialloop,1)
%     if rn_partialloop(links_partialloop(i,1),5) ~= rn_partialloop(links_partialloop(i,2),5)
%         idx_links_partialloop_removeinterloop_inverse = [idx_links_partialloop_removeinterloop_inverse;i];
%         continue
%     end
%     idx_links_partialloop_removeinterloop = [idx_links_partialloop_removeinterloop;i];
% end
% 
% idx_reconstruct_links_partialloop = [1;idx_links_partialloop_removeinterloop_inverse;size(rn_partialloop,1)];
% links_partialloop = links_partialloop(idx_links_partialloop_removeinterloop,:);


links_partialloop = links;
rn_partialloop = rn;

% remove links that connect different loops
links_onsurface = [];
idx_links_onsurface = [];
idx_links_onsurface_inverse = [];
for i = 1:size(links_partialloop,1)
    if rn_partialloop(links_partialloop(i,1),4) == 4 && rn_partialloop(links_partialloop(i,2),4) == 4
        idx_links_onsurface_inverse = [idx_links_onsurface_inverse;i];
        continue
    end
    idx_links_onsurface = [idx_links_onsurface;i];
end
links_onsurface = links_partialloop(idx_links_onsurface,:);


% find intersects of links with GB, and use the intersects in partial loops
for i = 1:size(links_onsurface)
    idx1 = links_onsurface(i,1);
    idx2 = links_onsurface(i,2);
    % links that cross the GB
    if rn_partialloop(idx1,4) ~= rn_partialloop(idx2,4)
        rn1 = rn_partialloop(idx1,1:3);
        rn2 = rn_partialloop(idx2,1:3);
        rn1 = [rn1 1];
        rn2 = [rn2 1];
        delta_r = rn2-rn1;
        k = (-rn1*SurfacePlane')/((delta_r*SurfacePlane'));
        r0 = rn1+k*delta_r;
        idx_segment = rn_partialloop(idx1,5);
        r0 = [r0(1:3) 4 idx_segment];

        if rn_partialloop(idx1,4) == 4
            rn_partialloop(idx1,:) = r0;
        end
        if rn_partialloop(idx2,4) == 4
            rn_partialloop(idx2,:) = r0;
        end

    end
end


rn_partialloop(1,4) = 4;
rn_partialloop(end,4) = 4;
rn = rn_partialloop;
rn(:,4) = 4;
links = links_onsurface;

rn_init = rn;
links_init = links;


% 1 eV/Å3 = 160.2176621 GPa
multiplier_GPaToeV = 1/160.2176621;


% Inclusion
doinclusion=0;
inclusion_pos_rad=[0 8 0 Ri 0]; % 11 -11
maxconnections=500;
lmax = 4;
lmin = 1.25;
rann = a;
rntol = rann/2;
areaminmag(1)=2*lmin*rann*cos(asin(2*rann/(2*pi*lmin))); 
areaminmag(2)=2*lmax*rntol*cos(asin(2*rntol/(2*pi*lmax)));
areamin=min(areaminmag);
areamax=(1/4)*lmax*lmax;

%% simulation parameters
mobility='mobbcc_climb_with_Bt_varying';
% mobility='mobbcc0';
integrator='int_eulerbackward';
integrator='int_eulerbackward_prismatic';
dopartials=0;
stackforcevec=zeros(size(rn,1),3);
doSFT=0;
SFT_plane=[];
doshielding=0;
dobox=0;
boxD = [-boxlength boxlength; -boxlength boxlength; -boxlength boxlength];
docrossslip=0;
doremesh=0;
docollision=1;
doseparation=0;
% SurfacePlane=[0 0 0 0];
     

%rotation matrix
% Q = [ dot(e1,e1p) dot(e2,e1p) dot(e3,e1p)
%       dot(e1,e2p) dot(e2,e2p) dot(e3,e2p)
%       dot(e1,e3p) dot(e2,e3p) dot(e3,e3p) ];
%   
% %Transform sigma into coordinate system 1
% appliedstress = Q*sigma*Q';
%16.257e-14. = conversion factor for Aluminium
%appliedstress = -1500*13.0676e-14.*(1/(2*sqrt(6))).*([-2 0 -1; 0 2 1; -1 1 0]+(1/3).*[2 0 -1; 0 -2 1; -1 1 0]);

% viewangle = [1,3,1]; %Upward tilted - able to see top and climb
% viewangle=[-1,1,-0.5];  %Looking down z-axis - top view

%viewangle=[0,90];  %Looking down z-axis - top view
%viewangle = [-5,10,10]; %From behind the dislocation
%viewangle = [0,10,0]; %Looking straight at the percipitate - climb angle
% viewangle = [-8,0,0]; %profile view
% viewangle = [1,3,1]; %Upward tilted - able to see top and climb
 viewangle = [1,1,1]; %Upward tilted - able to see top and climb
 viewangle = [-1,-2,0.5]; %
%  viewangle = [0,2,1]; %
 viewangle = [1,0,0.1]; %
viewangle1 = [1,0,0.1]; %

viewangle2 = [-1,-2,1]; %


% x: glide direction y: superjog height z: dislocation length 
% XYZ: for plot purpose, x y z -> Z X Y
printfreq=10;      
printnode=10;
plotfreq=1;       
rmax=0.1;
f.Name = 'Initial Configuration';
f.Position = [100 100 1000 800];

fg1 = figure(1);
fg2 = figure(2);
set(0,'CurrentFigure',fg1);
set(gcf, 'Position', [10 400 600 600]);
plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)
% plotnodes(rn,links_init0,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle1)

set(0,'CurrentFigure',fg2);
set(gcf, 'Position', [700 400 600 600]);
% plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)
plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,viewangle2)


