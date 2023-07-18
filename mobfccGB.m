function [vn,fn] = mobfccGB(fseg,rn,links,connectivity,nodelist,conlist)
%mobility law function (model: FCC0)
global Bscrew Bedge Beclimb Bline

%Drag (Mobility) parameters (should be specified by Input file)
%Bscrew=1e0;
%Bedge=1e0;
%Beclimb=1e8;
%Bline=1.0e-2*min(Bscrew,Bedge);

%numerical tolerance
eps=1e-12;

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

% For transformation on one side of GB
%coordinate system 1 (cubic coordinate system)
e1p = [1 0 0]; e2p = [0 1 0]; e3p = [0 0 1];
e1p=e1p/norm(e1p); e2p=e2p/norm(e2p); e3p=e3p/norm(e3p);            
%coordinate system 2
e1 =  [1 0 0]; e2  = [0 1 1]; e3 = [0 -1 1];
e1=e1/norm(e1); e2=e2/norm(e2); e3=e3/norm(e3);
%rotation matrix
% Q = [ dot(e1,e1p) dot(e2,e1p) dot(e3,e1p)
%       dot(e1,e2p) dot(e2,e2p) dot(e3,e2p)
%       dot(e1,e3p) dot(e2,e3p) dot(e3,e3p) ];
theta = 45;  
Q = [ cos(theta) sin(theta) 0
      -sin(theta) cos(theta) 0 
      0 0 1 ];  
% now cycle through all of the nodes for which the velocity must be calculated

for n=1:L1
    n0=nodelist(n);                 %n0 is the nodeid of the nth node in nodelist
    numNbrs=conlist(n,1);           %numNbrs is the number of connections for node n0 in conlist
    fn(n,:)=zeros(1,3);             % initialize the total force and the total drag matrix
    Btotal=zeros(3,3);
    % Transform burgers vector and glide direction on one side of GB
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
             if rn(n0,1) < 0 || rn(n1,1) < 0
                 burgv = Q*burgv';
                 nplane = Q*nplane';
                 burgv = burgv';
                 nplane = nplane';
             end
            nmag=norm(nplane);
            if nmag<eps
                % the normal plane is not defined try to define the normal plane
                Btotal=Btotal+(2.0*L).*((Beclimb).*eye(3)+(Bline-Beclimb).*(linedir'*linedir));
            else
                nplane=nplane./nmag;
                % if abs(abs(nplane(1)*nplane(2))+abs(nplane(2)*nplane(3))+abs(nplane(1)*nplane(3))-1)>0.1 % this can be the first condition that checks for good glide planes
                if abs(nplane(1)*nplane(2)*nplane(3))<0.01 % this is a second condition that checks for good glide planes
                    % this is a special glide plane that is not of 111 type signifying a junction dislocation
                    Btotal=Btotal+(2.0*L).*((Beclimb).*eye(3)+(Bline-Beclimb).*(linedir'*linedir));
                else
                    cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
                    mdir=cross(nplane,linedir);
                    Bglide=1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2);
                    Btotal=Btotal+(2.0*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Beclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );
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