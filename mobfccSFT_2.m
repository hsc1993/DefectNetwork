function [vn,fn] = mobfcc3(fseg,rn,links,connectivity,nodelist,conlist,mobility,dopartials,stackforcevec,rntol)
%mobility law function (model: FCC0)

norm111plane=norm([1 1 1]);
%Drag (Mobility) parameters
Bscrew=12.84e-25;%IF THE UNIT LENGTH IS THE LATTICE PARAMETER (Enrique Sep 2005)
Bedge=3.0214e-25;
%Bclimb=1e-12;%Almost no climb at all

Bline=1.0e-2*min(Bscrew,Bedge);
%numerical tolerance

SFT=8; %Introduction of SFT volume
SFTtol=0.01;

eps=1e-15;
%rntol=1e-2;

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
    Mtotal=zeros(3,3);
    velassigned=0;
    nplane=zeros(numNbrs,3);
    for i=1:numNbrs
        ii=conlist(n,i+1);        % connectionid for this connection
        linkid(i,:)=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid(i,:),3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);% calculate the length of the link and its tangent line direction. Always outward the node we are watching at
        L=norm(rt);
        if L>0.0
            fsegn0=fseg(linkid(i,:),3*(posinlink-1)+[1:3]);
            fn(n,:)=fn(n,:)+fsegn0; % nodeid for the node that n0 is connected to
            burgv=links(linkid(i,:),3:5); %(Enri modify)% burgers vector of the link                                                           
            linedir=rt./L;
            
            nplane(i,:)=links(linkid(i,:),6:8)/norm(links(linkid(i,:),6:8));
            
            if ~(((-rn(n0,1)-rn(n0,2)-rn(n0,3))<=(-(SFT+SFTtol)))&((rn(n0,1)-rn(n0,2)+rn(n0,3))<=(SFT-SFTtol))&((-rn(n0,1)+rn(n0,2)+rn(n0,3))<=(SFT-SFTtol))&((rn(n0,1)+rn(n0,2)-rn(n0,3))<=(SFT-SFTtol)))
                if(norm(stackforcevec(linkid(i,:),:))>eps)
                    stackingplane=stackforcevec(linkid(i,:),:)/norm(stackforcevec(linkid(i,:),:));
            
                    if(norm(cross(nplane(i,:),stackingplane))<rntol)
                        cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
                        mdir=cross(nplane(i,:),linedir);
                        Bglide=2.081*(1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2));
                        Mglide=1/Bglide;
                        %Btotal=Btotal+(2.0*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Bclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );%( Bclimb ) .* ( nplane' * nplane ) +
                        Mtotal=Mtotal+(1/(2.0*L)).*( ( Mglide ).* ( mdir' * mdir ) + ( 1/Bline ).*( linedir' * linedir ) );
                    else
                        Mtotal=Mtotal+(1/(2.0*L)).*((1/Bline).*(linedir'*linedir));
                    end
                elseif ((abs(dot(nplane(i,:),[1 1 1]/norm([1 1 1])))>(1-rntol))|(abs(dot(nplane(i,:),[-1 1 1]/norm([-1 1 1])))>(1-rntol))|(abs(dot(nplane(i,:),[1 -1 1]/norm([1 -1 1])))>(1-rntol))|(abs(dot(nplane(i,:),[1 1 -1]/norm([1 1 -1])))>(1-rntol)))
                    cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
                    mdir=cross(nplane(i,:),linedir);
                    Bglide=2.081*(1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2));
                    Mglide=1/Bglide;
                    %Btotal=Btotal+(2.0*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Bclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );%( Bclimb ) .* ( nplane' * nplane ) +
                    Mtotal=Mtotal+(1/(2.0*L)).*( ( Mglide ).* ( mdir' * mdir ) + ( 1/Bline ).*( linedir' * linedir ) );
                else
                    Mtotal=Mtotal+(1/(2.0*L)).*((1/Bline).*(linedir'*linedir));
                end
            else
               Mtotal=zeros(3,3);%Mtotal+(1/(2.0*L)).*((1/Bline).*(linedir'*linedir));
               break;
           end
        end
        %if velassigned==1
        %    break;
        %end
    end
%     if (rcond(Mtotal)<eps)
%         
%         [evec,eval]=eig(Mtotal);                    % find eigenvalues and eigen vectors of drag matrix
%         evalmax=eval(1,1);
%         eval=eval./evalmax;
%         fvec=fn(n,:)'./evalmax;
%         for i=2:3                                   % invert drag matrix and keep zero eigen values as zero
%             if eval(i,i)>eps
%                 eval(i,i)=1/eval(i,i);
%             else
%                 eval(i,i)=0.0d0;
%             end
%         end
%         vn(n,:)=(evec*eval*evec'*fvec)';  % calculate the velocity
%     else
    vn(n,:)=(Mtotal*fn(n,:)')';                 % Btotal was wellconditioned so just take the inverse
        
    for i=1:numNbrs
        if(norm(vn(n,:))~=0)
            if(abs(dot(vn(n,:)/norm(vn(n,:)),nplane(i,:)))>rntol)
                vn(n,:)=[0 0 0];
                break;
            end
        end
    end
    %end
    
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