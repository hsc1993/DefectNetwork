function [vn,fn] = mobfcc3(fseg,rn,links,connectivity,nodelist,conlist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad)
%mobility law function (model: FCC0)

norm111plane=norm([1 1 1]);
%Drag (Mobility) parameters
Bscrew=12.84e-25;%IF THE UNIT LENGTH IS THE LATTICE PARAMETER (Enrique Sep 2005)
Bedge=3.0214e-25;
%Bclimb=1e-12;%Almost no climb at all

%Bline=1.0e-2*min(Bscrew,Bedge);
Bline=1.0e1*min(Bscrew,Bedge);
%numerical tolerance
eps=1e-6;
%rntol=1e-2;

perfect_plane=[1 1 1; 1 1 -1;-1 1 1;1 -1 1]/norm([1 1 1]);
perfmatrix=[1 -1 0;-1 0 1;0 -1 1;0 1 1;1 0 1;1 1 0]./norm([1 1 0]);
numperf=size(perfmatrix,1); 
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
    fn(n,:)=zeros(1,3); % initialize the total force and the total drag matrix
    vn(n,:)=zeros(1,3);
    fn_newbase(n,:)=zeros(1,3);
    Mtotal=zeros(3,3);
    rt=zeros(numNbrs,3);
    nplane=zeros(numNbrs,3);
    for i=1:numNbrs
        ii=conlist(n,i+1);        % connectionid for this connection
        linkid(i,:)=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid(i,:),3-posinlink);
        rt(i,:)=(rn(n1,1:3)-rn(n0,1:3));% calculate the length of the link and its tangent line direction. Always outward the node we are watching at
        L(i)=norm(rn(n1,1:3)-rn(n0,1:3));
        if L(i)>0.0
            fsegn0=fseg(linkid(i,:),3*(posinlink-1)+[1:3]);
            fn(n,:)=fn(n,:)+fsegn0; % nodeid for the node that n0 is connected to
            burgv=links(linkid(i,:),3:5); %(Enri modify)% burgers vector of the link                                                           
            linedir=rt(i,:)/L(i);
            
            glideplane=cross(burgv,linedir);% Enrique Sep 2005
            nmag=norm(glideplane);
            if nmag<eps
                nplane(i,:)=AssignGlidePlaneScrewfcc(glideplane);
            else
                nplane(i,:)=glideplane/nmag;
            end
            if(doSFT)
                check=CheckSFT(rn,rn(n0,1:3),doSFT,SFT_plane);
            else
                check=0;
            end
            
            g_dir=cross(nplane(i,:),linedir)/norm(cross(nplane(i,:),linedir));
            n_dir=nplane(i,:);
            Base=[linedir',g_dir',n_dir'];
            fsegn0_newbase=(inv(Base)*fsegn0')';
            cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
           % mdir=cross(nplane(i,:),linedir);
            Bglide=2.081*(1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2));
            Mglide=1/Bglide;
            if(check==1)
                Mglide=0;
            end
            Mline=1/Bline;
                    %Btotal=Btotal+(2.0*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Bclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );%( Bclimb ) .* ( nplane' * nplane ) +
            Mtotal=(1/(2.0*L(i)))*[Mline 0 0;0 Mglide 0;0 0 0];
            vn0_newbase=(Mtotal*fsegn0_newbase')';
            vn0=(Base*vn0_newbase')';
            vn(n,:)=vn(n,:) + vn0;
            
        end
       
    end
    check_glideplane=0;
    for i=1:size(perfect_plane,1)
        if(abs(dot(vn(n,:)/norm(vn(n,:)),perfect_plane(i,:)))<eps)
            check_glideplane=1;
            break;
        end
    end
    if(check_glideplane==0)
        Powermax=0;
        vn_p=[0 0 0];
        for i=1:numNbrs
            compatible=0;
            if(L(i)>0.0)
                vn_proj(i,:)=(dot(vn(n,:),rt(i,:)/L(i)))*rt(i,:)/L(i);
            else
                vn_proj(i,:)=[0 0 0];
            end
            for j=1:numNbrs
                if(norm(vn_proj(i,:))>0.0)
                    if~((abs(dot(vn_proj(i,:)/norm(vn_proj),nplane(j,:)))<eps))
                        compatible=1;
                        break;
                    end
                end
            end
            if(compatible==0)
                Powertest=abs(dot(vn_proj(i,:),fn(n,:)));
                if(Powertest>Powermax)
                    vn_p=vn_proj(i,:);
                end
            end
        end
        vn(n,:)=vn_p;
    end
                
            
    lrn2=size(rn,2);
    if((rn(n0,lrn2)~=0))
        vn(n,:)=[0 0 0];
    end
end