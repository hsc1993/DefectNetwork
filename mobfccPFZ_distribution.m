function [vn,fn,rnnew,ppt] = mobfccPFZ(fseg,rn,links,connectivity,nodelist,conlist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,ppt,dt)
%mobility law function (model: FCC0)
global Bscrew Bedge Bclimb Bline

%Drag (Mobility) parameters
%Bscrew=7.5e-25;%IF THE UNIT LENGTH IS THE LATTICE PARAMETER (Enrique Sep 2005)
%Bedge=3.9e-25;
%Beclimb=1e-12;%Almost no climb at all
%Bline=1.0e-2*min(Bscrew,Bedge);
%numerical tolerance
eps=1e-12;
rnnew=rn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables precipitate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vf = 4.9*10^(-8); % Number density [1/Å^3]
r_pfz = 1500; % Å, with pfz MOVE
mean = 18.127; % mean valu [nm]   
mean = mean*10; % Convert to Å
std = 7.583; %   
mu = 2.817; %   
sigma2 = 0.161;  
theta = 0.6155; % Angle between dislocation glide plane and precipitate (35.3 deg)
corr_factor = 0.5; % Correlation factor between precipitate area and length
ppt_n = [1 0 0]; % Orientation precipitates
ppt_new = []; % Initialize list of pinned precipitates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% length of the nodelist for which the velocity will be calculated
L1=size(nodelist,1);
% if no nodelist is given then the nodelist becomes the whole node population
% this portion of the code sets up that nodelist along with the connlist
% that contains all of the nodal connections
if L1==0
    L1=size(rn,1);
    nodelist=linspace(1,L1,L1)';
    [L2,L3]=size(connectivity);
    conlist=zeros(L2,(L3-1)/2+1);
    conlist(:,1)=connectivity(:,1);
    for i=1:L2
        connumb=conlist(i,1);
        conlist(i,2:connumb+1)=linspace(1,connumb,connumb);
    end
end

% now cycle through all of the nodes for which the velocity must be calculated
for n=1:L1
    n0=nodelist(n);                 %n0 is the nodeid of the nth node in nodelist
    numNbrs=conlist(n,1);           %numNbrs is the number of connections for node n0 in conlist
    fn(n,:)=zeros(1,3);             % initialize the total force and the total drag matrix
    Btotal=zeros(3,3);
    
    % Loop through neighboring nodes and calculate drag matrix
    for i=1:numNbrs
        ii=conlist(n,i+1);        % connectionid for this connection
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);% calculate the length of the link and its tangent line direction. Always outward the node we are watching at
        L=norm(rt);
        if L>0.0
            fsegn0=fseg(linkid,3*(posinlink-1)+[1:3]);
            fn(n,:)=fn(n,:)+fsegn0; % nodeid for the node that n0 is connected to
            burgv=links(linkid,3:5); %(Enri modify)% burgers vector of the link                                                           
            linedir=rt./L;
            
            glideplane=cross(burgv,linedir);% Enrique Sep 2005
            nmag=norm(glideplane);
            if nmag<eps
                nplane=links(linkid,6:8);
                nmag=norm(nplane);
                if nmag>eps
                    nplane=nplane./nmag;
                end
            else
                nplane=glideplane/nmag;
            end
            % if abs(abs(nplane(1)*nplane(2))+abs(nplane(2)*nplane(3))+abs(nplane(1)*nplane(3))-1)>0.1 % this can be the first condition that checks for good glide planes
            if abs(nplane(1)*nplane(2)*nplane(3))<0.1 % this is a second condition that checks for good glide planes
                % this is a special glide plane that is not of 111 type signifying a junction dislocation
                Btotal=Btotal+(2.0*L).*((Bclimb).*eye(3)+(Bline-Bclimb).*(linedir'*linedir));
            else
                cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
                mdir=cross(nplane,linedir);
                Bglide=2.081*(1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2));
                Btotal=Btotal+(2.0*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Bclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );%( Bclimb ) .* ( nplane' * nplane ) +
            end   
        end
    end
    
    % If already pinned by precipitate
    if rn(n0,4)==3 && L1==1 %Added node in remesh next to pinned node do not inherit flag
        rn(n0,4)=0;
    end
    
    % Calculate force for pinned node
    if rn(n0,4)==3 && L1>1
        r_ppt=ppt(find(ppt(:,1)==n),2);
        s_ppt = 0.05*r_ppt; % Stress to break presipitate, need to be taken from table or function depending on presipitate size! UPDATE!
        sn = s_ppt*burgv'; % Need to be a tensor with s_ppt in glide direction
        f_ppt = cross(sn',linedir);
        fn(n,:)=fn(n,:)-f_ppt;
    end
    
    % Calculate nodal velocities
    if rcond(Btotal)<eps
        [evec,eval]=eig(Btotal);                    % find eigenvalues and eigen vectors of drag matrix
        evalmax=eval(1,1);
        eval=eval./evalmax;
        fvec=fn(n,:)'./evalmax;
        for i=2:3                                   % invert drag matrix and keep zero eigen values as zero
            if eval(i,i)>eps
                eval(i,i)=1/eval(i,i);
            else
                eval(i,i)=0.0d0;
            end
        end
        vn(n,:)=(evec*eval*evec'*fvec)';  % calculate the velocity
    else
        vn(n,:)=(Btotal\fn(n,:)')';                 % Btotal was wellconditioned so just take the inverse
    end
    
    if rn(n0,4)==3 && L1>1
        disp(sprintf('Ppt pinned node %i: r=%d f_ppt=%d fn=%d',n0,r_ppt,norm(f_ppt),norm(fn(n,:))));
        
        % Unpin pinned node if driving force is high enough
        if norm(fn(n,:))>norm(f_ppt)
            disp(sprintf('Unpinned node %i!',n0));
            rnnew(n0,4)=0;
        else
            % Save pinned precipitate if not
            ppt_new=cat(1,ppt_new,[n0,r_ppt]);
        end
        
    end
%    if numNbrs==2
%        ii=conlist(n,2);                                                                      
%        n1=links(connectivity(n0,2*ii),3-connectivity(n0,2*ii+1));
%        ii=conlist(n,3);                                                                      
%        n2=links(connectivity(n0,2*ii),3-connectivity(n0,2*ii+1));
%        rt=rn(n1,1:3)-rn(n2,1:3);
%        L=norm(rt);
%        linedir=rt./L;
%        vn(n,:)=((eye(3)-linedir'*linedir)*vn(n,:)')';
%    end
    lrn2=size(rn,2);
    if(rn(n0,lrn2)~=0)
        vn(n,:)=[0 0 0];
    end
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for new presipitates %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % TO DO
    % - Make sure units are correct!
    % - Lognormal distribution radius instead of length
    % - Update stress nessescary to break precipitate depending on radius
    % - What to do with velocities conserning shearable/non-shearable?
    % - Improve pfz definition

% Loop through nodes and check if any hits precipitates
for n=1:L1
    
    n0=nodelist(n);                 %n0 is the nodeid of the nth node in nodelist
    numNbrs=conlist(n,1);           %numNbrs is the number of connections for node n0 in conlist

    Btotal=zeros(3,3);
    
    % Skip if outside pfz or if already pinned
    if rn(n0,1) > r_pfz && rn(n0,4)~=7 && L1>1 && rn(n0,4)~=3
        
        vol=0.0; 
        rn0i1=vn(n0,:)*dt;
        
        % Loop throuh connected nodes and calculate swiped volume
        for i=1:numNbrs
            ii=conlist(n,i+1);        % connectionid for this connection
            linkid=connectivity(n0,2*ii);
            posinlink=connectivity(n0,2*ii+1);
            n1=links(linkid,3-posinlink);
            rt=rn(n1,1:3)-rn(n0,1:3); % calculate the length of the link and its tangent line direction. Always outward the node we are watching at
            L=norm(rt);
            if L>0.0
                % Calculate volume swiped my segment
                rn1i1=vn(n1,:)*dt;
                rti1=(rn(n1,1:3)+vn(n1,:)*dt)-(rn(n0,1:3)+vn(n0,:)*dt); 
                vI=0.5*det(cat(1,rt,cat(1,rn0i1,ppt_n)));
                vII=0.5*det(cat(1,rn1i1,cat(1,rti1,ppt_n)));
                vol=vol+0.5*norm(vI+vII);
            end
            

        end
               
        p_ppt=rand; % Draw random number
        Pi = vf*vol*mean;
        
        % Check if hitting a presipitate
        if p_ppt < Pi

            % Calculate precipitate radius
            rand_nr = rand;  
            length=exp(-sqrt(2*sigma2)*erfcinv(2*rand_nr)+mu)*10;
            
            if rn(n0,4)==3
                r_ppt=ppt(find(ppt(:,1)==n),2);
            else
                r_ppt = sqrt(length*corr_factor/pi); % Radius precipitate 
            end
            
            ppt_new=cat(1,ppt_new,[n0,r_ppt]);
            
            % Add stress nessescary to break presipitate
            s_ppt = 0.05*r_ppt; % Stress to break presipitate, need to be taken from table or function depending on presipitate size! UPDATE!  
            sn = s_ppt*burgv'; % Need to be a tensor with s_ppt in glide direction
            f_ppt = cross(sn',linedir);
            fn(n,:)=fn(n,:)-f_ppt;
            if norm(fn(n,:))>norm(f_ppt)
                fn(n,:)=[0 0 0];
            end
            rnnew(n0,4)=3; % Mark node as precipitate node
            
            % Calculate drag matrix
            for i=1:numNbrs
                ii=conlist(n,i+1);        % connectionid for this connection
                linkid=connectivity(n0,2*ii);
                posinlink=connectivity(n0,2*ii+1);
                n1=links(linkid,3-posinlink);
                rt=rn(n1,1:3)-rn(n0,1:3);% calculate the length of the link and its tangent line direction. Always outward the node we are watching at
                L=norm(rt);
                if L>0.0
                    %fsegn0=fseg(linkid,3*(posinlink-1)+[1:3]);
                    %fn(n,:)=fn(n,:)+fsegn0; % nodeid for the node that n0 is connected to
                    burgv=links(linkid,3:5); %(Enri modify)% burgers vector of the link                                                           
                    linedir=rt./L;

                    glideplane=cross(burgv,linedir);% Enrique Sep 2005
                    nmag=norm(glideplane);
                    if nmag<eps
                        nplane=links(linkid,6:8);
                        nmag=norm(nplane);
                        if nmag>eps
                            nplane=nplane./nmag;
                        end
                    else
                        nplane=glideplane/nmag;
                    end
                    % if abs(abs(nplane(1)*nplane(2))+abs(nplane(2)*nplane(3))+abs(nplane(1)*nplane(3))-1)>0.1 % this can be the first condition that checks for good glide planes
                    if abs(nplane(1)*nplane(2)*nplane(3))<0.1 % this is a second condition that checks for good glide planes
                        % this is a special glide plane that is not of 111 type signifying a junction dislocation
                        Btotal=Btotal+(2.0*L).*((Bclimb).*eye(3)+(Bline-Bclimb).*(linedir'*linedir));
                    else
                        cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
                        mdir=cross(nplane,linedir);
                        Bglide=2.081*(1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2));
                        Btotal=Btotal+(2.0*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Bclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );%( Bclimb ) .* ( nplane' * nplane ) +
                    end   
                end
            end
            
            % Calculate new velocities 
            if rcond(Btotal)<eps
                [evec,eval]=eig(Btotal);                    % find eigenvalues and eigen vectors of drag matrix
                evalmax=eval(1,1);
                eval=eval./evalmax;
                fvec=fn(n,:)'./evalmax;
                for i=2:3                                   % invert drag matrix and keep zero eigen values as zero
                    if eval(i,i)>eps
                        eval(i,i)=1/eval(i,i);
                    else
                        eval(i,i)=0.0d0;
                    end
                end
                vn(n,:)=(evec*eval*evec'*fvec)';  % calculate the velocity
            else
                vn(n,:)=(Btotal\fn(n,:)')';                 % Btotal was wellconditioned so just take the inverse
            end 
            
            % Check if driving force is high enough to break precipitate
            disp(sprintf('Ppt hit node %i: r=%d f_ppt=%d fn=%d ',n0,r_ppt,norm(f_ppt),norm(fn(n,:))));
            if norm(fn(n,:))>norm(f_ppt)
                disp(sprintf('Unpinned node %i!',n0));
                rnnew(n0,4)=0;
                continue;
            end
        end
    end
    
    lrn2=size(rn,2);
    if(rn(n0,lrn2)~=0)
        vn(n,:)=[0 0 0];
    end
   
end
ppt=ppt_new;