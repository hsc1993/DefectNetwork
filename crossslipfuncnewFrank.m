%This function handle the cross-slip of the perfect and partials dislocations fro fcc metals

function [rn,links,connectivity,linksinconnect,fseg,stackforcevec,crossslip,printnode]=crossslipfuncnewFrank(rn,links,connectivity,linksinconnect,fseg,stackforcevec,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,crossslip,integrator,dt,dt0,rmax,totaltime,areamin,lmin,printnode,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding)
%Once we have constriction


maxconnections=8;
eps=1e-6;
partnorm=norm(1/6.*[1 1 2]);
perfnorm=norm(1/2.*[1 1 0]);
franknorm=norm(1/3.*[1 1 1]);
[lrn,lrn2]=size(rn);
[llink,llink2]=size(links);
reset=0;

for i=1:llink
    tangent=rn(links(i,2),1:3)-rn(links(i,1),1:3);
    length=norm(tangent);
    tangnorm=tangent./length;
    burgvec=links(i,3:5);
    nplane=links(i,6:8);
    if((abs(norm(burgvec)-franknorm)<eps)&(length>2*lmin))
        [rn,links,connectivity,linksinconnect,fseg,stackforcevec,printnode]=splitfrankpartialnew(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgvec,tangent,i,nplane,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,areamin,lmin,printnode,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);
    elseif((abs(norm(burgvec)-partnorm)<eps)&(length>3*lmin))
      %  [rn,links,connectivity,linksinconnect,fseg,stackforcevec,crossslip]=reactionpartialstairodlink(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgvec,tangent,i,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,crossslip,integrator,dt,dt0,rmax,totaltime);
        
    end
end

[lrn,lrn2]=size(rn);
i=1;
while (i<=lrn)
    reset=0;
    doesmove=rn(i,lrn2);
    if(doesmove==0)
    
        c=connectivity(i,1);
        
        linkid=zeros(1,c);
        tangent=zeros(c,3);
        burgv=zeros(c,3);
    
        for j=1:c
            linkid(j)=connectivity(i,2*j);
            tangent(j,:)=(rn(links(linkid(j),2),1:3)-rn(links(linkid(j),1),1:3));
            burgv(j,:)=links(linkid(j),3:5);
        end
    
        
   
    
        if((c==2)&(abs(norm(burgv(1,:))-partnorm)<eps)&(abs(norm(burgv(2,:))-partnorm)<eps)&(doesmove==0))
            
           %  [rn,links,connectivity,linksinconnect,fseg,stackforcevec,crossslip]=reactionpartialstairod(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgv,tangent,i,linkid,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,crossslip,integrator,dt,dt0,rmax,totaltime,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);
            
        elseif((c>2)&(doesmove==0))
            
            % [rn,links,connectivity,linksinconnect,fseg,stackforcevec,printnode]=crossmultinode(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgv,tangent,i,linkid,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,areamin,lmin,printnode,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);
            
        elseif((c==2)&(abs(norm(burgv(1,:))-perfnorm)<eps)&(abs(norm(burgv(2,:))-perfnorm)<eps)&(doesmove==0))
            
             [rn,links,connectivity,linksinconnect,fseg,stackforcevec,reset]=splitperfect(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgv,tangent,i,linkid,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);
             
                
        end
                
    end  
    
    if(reset==1)
        i=1;
        [lrn,lrn2]=size(rn);
    else
        i=i+1;
    end
    
end

        
        
    