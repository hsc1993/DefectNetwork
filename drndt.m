function [vnvec,fn,fseg,flag] = drndt(t,rnvec,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt)

%unscramble rn
rn=[reshape(rnvec,length(rnvec)/3,3),flag];

%rn(:,1:3)

%nodal driving force
fseg=segforcevec(MU,NU,a,Ec,rn,links,appliedstress,[],mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);


% %% calculate field stress influence on each node
% b_list = [];
% seg_start_list = [];
% seg_end_list = [];
% mu=MU;
% nu=NU;
% 
% for j = 1:size(links,1)
%     b_list = [b_list;links(j,3:5)];
%     id_rnstart = links(j,1);
%     id_rnend = links(j,2);
%     seg_start_list = [seg_start_list;rn(id_rnstart,1:3)];
%     seg_end_list = [seg_end_list;rn(id_rnend,1:3)];
% end
% 
% rn_aditional_force = [];
% for i = 1:size(links,1)   
%     idx_r1 = links(i,1);
%     idx_r2 = links(i,2);
%     fieldPointStress = FieldPointStress(rn(idx_r1,1:3),seg_start_list,seg_end_list,b_list,a,mu,nu); % same as [mu] -> N/A^2
%     stress =  [fieldPointStress(1) fieldPointStress(4) fieldPointStress(6);
%             fieldPointStress(4) fieldPointStress(2) fieldPointStress(5);
%             fieldPointStress(6) fieldPointStress(5) fieldPointStress(3)];
%     r1= rn(idx_r1,1:3);
%     r2= rn(idx_r2,1:3);
%     l = norm(r2-r1);
%     t = (r2-r1)/l; % line direction unit vector
%     b = links(i,3:5);
%     g = b/norm(b);
%     n = cross(t,g);
% 
%     stress_b = [stress(1,1)*b(1)+stress(1,2)*b(1)+stress(1,3)*b(1)
%         stress(2,1)*b(2)+stress(2,2)*b(2)+stress(2,3)*b(2)
%         stress(3,1)*b(3)+stress(3,2)*b(3)+stress(3,3)*b(3)];
%                     
%     line_vector = r2-r1;
%     fpk = cross(stress_b',line_vector);
%     fpk_half = fpk/2;
% 
%     rn_aditional_force = [rn_aditional_force;fpk_half];
%     fseg(i,1:3) = fseg(i,1:3)+fpk_half;
%     fseg(i,4:6) = fseg(i,4:6)+fpk_half;
% end


%%
% rn_aditional_force
% 
% fseg
% 
% size(rn_aditional_force)
% size(fseg)
% fseg = fseg+rn_aditional_force


%mobility function
% [vn,fn]=feval(mobility,fseg,rn,links,connectivity,[],[],mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt);
[vn,fn]=feval(mobility,fseg,rn,links,connectivity,[],[],mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt);

%fixed nodes (flag~=0) are not allowed to move
vn=vn.*((rn(:,4)==0 | rn(:,4)==3)*[1 1 1]);

%nodes on precipitate are not allowed to move
% num_inclusions=size(inclusion_pos_rad,1);
% for i=1:num_inclusions
%     vn=vn.*((sqrt(sum((rn(:,1:3)-repmat(inclusion_pos_rad(i,1:3),size(rn,1),1)).^2,2))>repmat(inclusion_pos_rad(i,4),size(rn,1),1))*[1 1 1]);
% end
    
%if flag==4 then is a surface node and can move only in the surface plane.
%The normal to the surface plane is estimated with the links to the surface
%node
vn = image_vel(rn,vn,SurfacePlane);
flag=rn(:,4);
%make up the last
vnvec=reshape(vn,length(vn(:,1))*3,1);

%vn
%disp(sprintf('t=%20.10e vn(1,:)=(%20.10e,%20.10e,%20.10e)',t,vn(1,1),vn(1,2),vn(1,3)));
%pause(0.2)
%if(t>0)
%   pause
%end