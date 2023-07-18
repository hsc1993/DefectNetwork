function [crossslipstressplus]=CrossSlipTest(rn, links, stackforcevec, fseg,connectivity,linksinconnect, maxconnections,MU,NU,a,Ec,mobility,dopartials,appliedstress,integrator,dt,dt0,rmax,rntol,lmin,lmax,areamin,areamax,rann,gamma)


glideplane=[1 1 1]/norm([1 1 1]);
tangent=[1 -1 0]/norm([1 -1 0]);
eps=1e-15;
timesteps=30;
checknode=2;
forceproject=cross(glideplane,tangent)/norm(cross(glideplane,tangent));
%firsttry=[-2 2 4];
maxdist=0;%zeros(3,3);
mindist=0;%zeros(3,3);
equilibrium=0;
count=1;
rnold=rn;
linksold=links;
stackforcevecold=stackforcevec;
totaltime=0;
printnode=2;


%Cross-slip stress

crossslip=0;
crossslip2=0;
count=0;
crossmax=20000;
crossslipstressplus=crossmax;
othercrossplus=0;
appliedstressold=appliedstress;
appliedstress=appliedstress+crossmax*13.0676e-14.*(1).*[0.6124  0   -0.6124;0   -0.6124    0.6124;-0.6124    0.6124  0];
rnold=rn;
linksold=links;
stackforcevecold=stackforcevec;
fsegold=fseg;
connectivityold=connectivity;
linksinconnectold=linksinconnect;
while(crossslip==0)
    rn=rnold;
    links=linksold;
    stackforcevec=stackforcevecold;
    fseg=fsegold;
    connectivity=connectivityold;
    linksinconnect=linksinconnectold;
    [lrn,lrn2]=size(rn);
    %crossslipstressplus=crossslipstressplus+1;
    %appliedstress=appliedstress+13.0676e-14.*[1/12  -1/12  -2/12; -1/12  0  1/12; -2/12  1/12  3/12];
    [fseg]=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,[],mobility,dopartials,stackforcevec);
    [rn,links,connectivity,linksinconnect,fseg,stackforcevec,crossslip,printnode]=crossslipfuncnewFrank(rn,links,connectivity,linksinconnect,fseg,stackforcevec,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,crossslip,integrator,dt,dt0,rmax,totaltime,areamin,lmin,printnode);
    rn=rnold;
    links=linksold;
    stackforcevec=stackforcevecold;
    fseg=fsegold;
    connectivity=connectivityold;
    linksinconnect=linksinconnectold;
    [lrn,lrn2]=size(rn);
    
    appliedstress=appliedstress+13.0676e-14.*(1).*[0.6124  0   -0.6124;0   -0.6124    0.6124;-0.6124    0.6124  0];
    [fseg]=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),linksold,appliedstress,[],mobility,dopartials,stackforcevecold);
    [rn,links,connectivity,linksinconnect,fseg,stackforcevec,crossslip2,printnode]=crossslipfuncnewFrank(rn,linksold,connectivityold,linksinconnectold,fseg,stackforcevecold,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,crossslip2,integrator,dt,dt0,rmax,totaltime,areamin,lmin,printnode);
        
        %disp(crossslipstressplus);
    
    count=count+1;
    if(crossslip==crossslip2)
        
        if(crossslip==1)
            crossslip=0;
            crossslip2=0;
            crossmax=crossslipstressplus;
            crossslipstressplus=(crossmax+othercrossplus)/2;
            %crossslipstressplus=abs(crossslipstressplus);
            appliedstress=appliedstressold+crossslipstressplus*13.0676e-14.*(1).*[0.6124  0   -0.6124;0   -0.6124    0.6124;-0.6124    0.6124  0];
            
        
        elseif(crossslip==0)
            othercrossplus=crossslipstressplus;
            crossslipstressplus=(crossmax+othercrossplus)/2;
            appliedstress=appliedstressold+crossslipstressplus*13.0676e-14.*(1).*[0.6124  0   -0.6124;0   -0.6124    0.6124;-0.6124    0.6124  0];           
        end
    else
        crossslip=1;
    end
    disp(crossslipstressplus);
end
disp(sprintf('the Fleischer acute stress is %e',crossslipstressplus));


