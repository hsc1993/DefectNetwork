function [vn,fn] = mobfcc0inga(fseg,rn,links,connectivity,nodelist,conlist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt)
%mobility law function (model: FCC0)

%Drag (Mobility) parameters
%Bscrew=12.84e-25;%IF THE UNIT LENGTH IS THE LATTICE PARAMETER (Enrique Sep 2005)
%Bedge=3.0214e-25;
%Bclimb=1e-12;%Almost no climb at all
Bscrew=1.9e-24; % Olmstead
Bedge=3.65e-24; % Olmstead
%Bscrew=2e-24;
%Bedge=6e-24;
Bclimb=1e-12;%Almost no climb at all
sf_barrier=1.2707e-10 + 0.204e-10;%1.2707e-10; %force per unit lenght to overcome a stacking fault plane + force due to the finite displacements. The nucleation barrier is missing though


%Bline=1.0e-2*min(Bscrew,Bedge);
Bline=1.0e-2*min(Bscrew,Bedge);
%numerical tolerance
eps=1e-16; %-18
TOL=1e-4;

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
    if(doinclusion)
        for k=1:size(inclusion_pos_rad,1)
            if(norm(rn(n0,1:3)-inclusion_pos_rad(k,1:3))<inclusion_pos_rad(k,4))
                Bscrew=1;
                Bedge=Bscrew;%1000*3.0214e-25;
                Bline=1.0e-2*min(Bscrew,Bedge);
            else
                Bscrew=1.9e-24; % Olmstead
                Bedge=3.65e-24; % Olmstead
                Bline=1.0e-0*min(Bscrew,Bedge);
            end
        end
    end
    numNbrs=conlist(n,1);           %numNbrs is the number of connections for node n0 in conlist
    fn(n,:)=zeros(1,3);             % initialize the total force and the total drag matrix
    Btotal=zeros(3,3);
    for i=1:numNbrs
        ii=conlist(n,i+1);                                                                      % connectionid for this connection
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);                                                               % calculate the length of the link and its tangent line direction
        L=norm(rt);
        if L>0.0
            fsegn0=fseg(linkid,3*(posinlink-1)+[1:3]);
            fn(n,:)=fn(n,:)+fsegn0; % nodeid for the node that n0 is connected to
            burgv=links(connectivity(n0,2*ii),3:5); % burgers vector of the link                                                           
            linedir=rt./L;
            nplane=links(linkid,6:8);
            nmag=norm(nplane);
            if nmag<eps
                % the normal plane is not defined try to define the normal plane
                Btotal=Btotal+(2.0*L).*((Bclimb).*eye(3)+(Bline-Bclimb).*(linedir'*linedir));
            else
                nplane=nplane./nmag;
                % if abs(abs(nplane(1)*nplane(2))+abs(nplane(2)*nplane(3))+abs(nplane(1)*nplane(3))-1)>0.1 % this can be the first condition that checks for good glide planes
                if (~((abs(dot(nplane,[1 0 0]/norm([1 0 0])))>0.95)|(abs(dot(nplane,[-1 0 0]/norm([-1 0 0])))>0.95)))%|...%The plane has to be a {111} plane % this is a second condition that checks for good glide planes
                % this is a special glide plane that is not of 111 type signifying a junction dislocation
                    Btotal=Btotal+(2.0*L).*((Bclimb).*eye(3)+(Bline-Bclimb).*(linedir'*linedir));
                else
                    cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
                    mdir=cross(nplane,linedir); %glidir
                    Bglide=1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2);
                    Btotal=Btotal+(2.0*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Bclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );
                end
            end
        end
    end
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
end