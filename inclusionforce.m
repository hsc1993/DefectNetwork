%This function gives the force on the nodes of a segment due to the elastic inclusion depending
%on where is the segment


function [finclusion_1, finclusion_2]=inclusionforcevec(MU,NU,segments,linkid,inclusion_pos_rad)

%This function adds the force due to a elastic spherical inclusion

lseg=size(segments,1);
num_inclusions=size(inclusion_pos_rad,1);
YoungModulus_2=2*(1+NU)*MU;
B_2=-(YoungModulus_2*(1-NU))/((1+NU)*(1-(2*NU)));
a=(3/2)*((1/2)-NU);
b=(1-2*NU);
c=NU - (1/2);
e=1-NU;

%d=1-3*NU;
finclusion_1=zeros(lseg,3);
finclusion_2=zeros(lseg,3);
if length(linkid)==0
    
    for i=1:num_inclusions
        finclusion_1_i=zeros(lseg,3);
        finclusion_2_i=zeros(lseg,3);
        center_inclusion_i=inclusion_pos_rad(i,1:3);
        R_i=inclusion_pos_rad(i,4);
        MU_i=inclusion_pos_rad(i,5);
        NU_i=inclusion_pos_rad(i,6);
        rho_i=inclusion_pos_rad(i,7);
        YoungModulus_i=2*(1+NU_i)*MU_i;
        B_i=-(YoungModulus_i*(1-NU_i))/((1+NU_i)*(1-(2*NU_i)));
        alpha_i=1-2*NU_i;
        beta_i=3*((1/4) - (1/2)*NU_i);
        gamma_i=1+3*NU_i;
        delta_i=(3/2) + 3*NU_i;
        for j=1:lseg
            b1=segments(j,3:5);
            x1=segments(j,6:8)-inclusion_pos_rad(i,1:3);
            x2=segments(j,9:11)-inclusion_pos_rad(i,1:3);
            zbase=(x2-x1)/norm((x2-x1));
            ybase=cross(zbase,x2)/norm(cross(zbase,x2));
            xbase=cross(ybase,zbase)/norm(cross(ybase,zbase));
            Base=[xbase',ybase',zbase'];
            b1_newbase=(inv(Base)*b1')';
            x1_newbase=(inv(Base)*x1')';
            x2_newbase=(inv(Base)*x2')';
            d=x1_newbase(3);
            H_j=x1_newbase(1);
            L_j=norm(x2-x1);
            norm_x1_newbase=norm(x1_newbase);
            norm_x2_newbase=norm(x2_newbase);
            if((norm(x1_newbase)<=R_i)&(norm(x2_newbase)<=R_i)) %The segment is completely inside of the inclusion
                [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
            elseif((norm(x1_newbase)>R_i)&(norm(x2_newbase)>R_i))
                if((R_i^2-x1_newbase(1)^2-x1_newbase(2)^2)<=0)
                    [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                elseif((R_i^2-x1_newbase(1)^2-x1_newbase(2)^2)>0)
                    x_cut_1(1:2)=x1_newbase(1:2);
                    x_cut_1(3)=sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
                    x_cut_2(1:2)=x1_newbase(1:2);
                    x_cut_2(3)=-sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
                    if(((x_cut_1(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_1(3)>max(x1_newbase(3),x2_newbase(3))))&((x_cut_2(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_2(3)>max(x1_newbase(3),x2_newbase(3)))))
                        [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                    else
                        L_inside=norm(x_cut_2 - x_cut_1);
                        L_outside=L_j - L_inside;
                        [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                        [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                        finclusion_1_i(j,:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                        finclusion_2_i(j,:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
                    end
                end
            else
                x_cut(1:2)=x1_newbase(1:2);
                x_cut(3)=sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
                if((x_cut(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut(3)>max(x1_newbase(3),x2_newbase(3))))
                    x_cut(3)=-sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
                end
                if(norm(x1_newbase)<=R_i)
                    L_inside=norm(x_cut-x1_newbase);
                    L_outside=L_j - L_inside;
                elseif(norm(x2_newbase)<=R_i)
                    L_inside=norm(x_cut-x2_newbase);
                    L_outside=L_j - L_inside;
                end
                [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                finclusion_1_i(j,:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                finclusion_2_i(j,:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
            end
            
        end
        finclusion_1=finclusion_1 + finclusion_1_i;
        finclusion_2=finclusion_2 + finclusion_2_i;
    end
    finclusion_1=(Base*finclusion_1')';
    finclusion_2=(Base*finclusion_2')';

else
    for i=1:num_inclusions
        finclusion_1_i=zeros(lseg,3);
        finclusion_2_i=zeros(lseg,3);
        center_inclusion_i=inclusion_pos_rad(i,1:3);
        R_i=inclusion_pos_rad(i,4);
        MU_i=inclusion_pos_rad(i,5);
        NU_i=inclusion_pos_rad(i,6);
        rho_i=inclusion_pos_rad(i,7);
        YoungModulus_i=2*(1+NU_i)*MU_i;
        B_i=-(YoungModulus_i*(1-NU_i))/((1+NU_i)*(1-(2*NU_i)));
        alpha_i=1-2*NU_i;
        beta_i=3*((1/4) - (1/2)*NU_i);
        gamma_i=1+3*NU_i;
        delta_i=(3/2) + 3*NU_i;
        for j=1:length(linkid)
            b1=segments(linkid(j),3:5);
            x1=segments(linkid(j),6:8)-inclusion_pos_rad(i,1:3);
            x2=segments(linkid(j),9:11)-inclusion_pos_rad(i,1:3);
            zbase=(x2-x1)/norm((x2-x1));
            ybase=cross(zbase,x2)/norm(cross(zbase,x2));
            xbase=cross(ybase,zbase)/norm(cross(ybase,zbase));
            Base=[xbase',ybase',zbase'];
            b1_newbase=(inv(Base)*b1')';
            x1_newbase=(inv(Base)*x1')';
            x2_newbase=(inv(Base)*x2')';
            d=x1_newbase(3);
            H_j=x1_newbase(1);
            L_j=norm(x2-x1);
            norm_x1_newbase=norm(x1_newbase);
            norm_x2_newbase=norm(x2_newbase);
            if((norm(x1_newbase)<=R_i)&(norm(x2_newbase)<=R_i)) %The segment is completely inside of the inclusion
                [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
            elseif((norm(x1_newbase)>R_i)&(norm(x2_newbase)>R_i))
                if((R_i^2-x1_newbase(1)^2-x1_newbase(2)^2)<=0)
                    [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                elseif((R_i^2-x1_newbase(1)^2-x1_newbase(2)^2)>0)
                    x_cut_1(1:2)=x1_newbase(1:2);
                    x_cut_1(3)=sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
                    x_cut_2(1:2)=x1_newbase(1:2);
                    x_cut_2(3)=-sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
                    if(((x_cut_1(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_1(3)>max(x1_newbase(3),x2_newbase(3))))&((x_cut_2(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_2(3)>max(x1_newbase(3),x2_newbase(3)))))
                        [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                    else
                        L_inside=norm(x_cut_2 - x_cut_1);
                        L_outside=L_j - L_inside;
                        [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                        [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                        finclusion_1_i(linkid(j),:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                        finclusion_2_i(linkid(j),:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
                    end
                end
            else
                x_cut(1:2)=x1_newbase(1:2);
                x_cut(3)=sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
                if((x_cut(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut(3)>max(x1_newbase(3),x2_newbase(3))))
                    x_cut(3)=-sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
                end
                if(norm(x1_newbase)<=R_i)
                    L_inside=norm(x_cut-x1_newbase);
                    L_outside=L_j - L_inside;
                elseif(norm(x2_newbase)<=R_i)
                    L_inside=norm(x_cut-x2_newbase);
                    L_outside=L_j - L_inside;
                end
                [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                finclusion_1_i(linkid(j),:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                finclusion_2_i(linkid(j),:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
            end
            
        end
        finclusion_1=finclusion_1 + finclusion_1_i;
        finclusion_2=finclusion_2 + finclusion_2_i;
    end
    finclusion_1=(Base*finclusion_1')';
    finclusion_2=(Base*finclusion_2')';
    
end

function [finclusion_1, finclusion_2]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j)


    eps=1e-6;
    
    finclusion_1=zeros(1,3);
    finclusion_2=zeros(1,3);
    
    b1=b1_newbase(1);
    b2=b1_newbase(2);
    b3=b1_newbase(3);
    
    x=x1_newbase(1);
    z_1=x1_newbase(3);
    z_2=x2_newbase(3);
    
    norm_r1=sqrt(x^2+z_1^2);
    norm_r2=sqrt(x^2+z_2^2);
    
    const_general=rho_i/(1-NU_i);
    
    
        %%%%%%%%%%
        %
        %   This section calculates the total force on the segment due to the inclusion
        %
        %%%%%%%%%%



    integral_node_1=zeros(1,3);
    integral_node_2=zeros(1,3);
    finclusion_total=zeros(1,3);
    
    f_total_1(1)=-(b2/5)*((H_j^2/R_i)*alpha_i*(z_2-z_1) - beta_i*H_j^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))));
    f_total_2(1)=(b2/10)*((gamma_i/R_i)*((z_2^3-z_1^3)/3 + x^2*(z_2-z_1)) - (delta_i/2)*(x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))) + (z_2*norm_r2 - z_1*norm_r1)));
    f_total_3(1)=-(b3/5)*((H_j/R_i)*(alpha_i/2)*(z_2^2-z_1^2) - beta_i*H_j*(norm_r2-norm_r1));
    
    
    finclusion_total(1)=const_general*(f_total_1(1) + f_total_2(1) + f_total_3(1));
    
    finclusion_total(2)=const_general*(b1/10)*((gamma_i/R_i)*((z_2^3-z_1^3)/3 + x^2*(z_2-z_1)) - (delta_i/2)*(x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))) + (z_2*norm_r2 - z_1*norm_r1)));
    
    finclusion_total(3)=0;

    %%%%%%%%%%
    %
    %   This section calculates the force on node 2
    %
    %%%%%%%%%%

    integral_node_1=zeros(1,3);
    integral_node_2=zeros(1,3);
    finclusion_2=zeros(1,3);
    
    f_2_1(1)=-(b2/5)*(H_j^2/L_j)*((alpha_i/R_i)*((z_2^2-z_1^2)/2 - d*(z_2-z_1)) - beta_i*((norm_r2-norm_r1) - d*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2)))));
    f_2_2(1)=(b2/10)*((gamma_i/(R_i*L_j)*((z_2^4-z_1^4)/4 - d*(z_2^3-z_1^3)/3 + x^2*(z_2^2-z_1^2)/2 - d*x^2*(z_2-z_1))) - (delta_i/(6*L_j)*(norm_r2*(2*norm_r2^2 - 3*d*z_2) - norm_r1*(2*norm_r1^2 - 3*d*z_1) - 3*d*x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))))));
    f_2_3(1)=-(b3/5)*(H_j/L_j)*((alpha_i/R_i)*((z_2^3-z_1^3)/3 - d*(z_2^2-z_1^2)/2) - beta_i*(norm_r2*(z_2/2 - d) - norm_r1*(z_1/2 - d)) - (x^2/2)*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))));
    
    finclusion_2(1)=const_general*(f_2_1(1) + f_2_2(1) + f_2_3(1));
    
    finclusion_2(2)=const_general*(b1/10)*((gamma_i/(R_i*L_j)*((z_2^4-z_1^4)/4 - d*(z_2^3-z_1^3)/3 + x^2*(z_2^2-z_1^2)/2 - d*x^2*(z_2-z_1))) - (delta_i/(6*L_j)*(norm_r2*(2*norm_r2^2 - 3*d*z_2) - norm_r1*(2*norm_r1^2 - 3*d*z_1) - 3*d*x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))))));

    finclusion_2(3)=0;


    %%%%%%%%%%
    %
    %   This section calculates the force on node 1
    %
    %%%%%%%%%%

    finclusion_1=finclusion_total - finclusion_2;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [finclusion_1, finclusion_2]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j)

    eps=1e-6;
    
    finclusion_1=zeros(1,3);
    finclusion_2=zeros(1,3);
    
    b1=b1_newbase(1);
    b2=b1_newbase(2);
    b3=b1_newbase(3);
    
    x=x1_newbase(1);
    z_1=x1_newbase(3);
    z_2=x2_newbase(3);
    
    norm_r1=sqrt(x^2+z_1^2);
    norm_r2=sqrt(x^2+z_2^2);
    
    rho_o=(1/5)*rho_i*B_i/B_2;
    
    const_general=rho_o*R_i^3/(1-NU);

    %%%%%%%%%%
    %
    %   This section calculates the total force on the segment due to the inclusion
    %
    %%%%%%%%%%

    integral_node_1=zeros(1,3);
    integral_node_2=zeros(1,3);
    finclusion_total=zeros(1,3);
    
    f_total_1(1)=b2*H_j^2*((R_i*alpha_i/(2*x^3))*(x*(z_2/norm_r2^2 - z_1/norm_r1) + atan(z_2/x) - atan(z_1/x)) - (beta_i/x^2)*(z_2/norm_r2 - z_1/norm_r1));
    f_total_2(1)=(b2/2)*((R_i*gamma_i/x^2)*(z_2/norm_r2 - z_1/norm_r1) + (delta_i/x)*(atan(z_2/x) - atan(z_1/x)));
    f_total_3(1)=b3*H_j*(-(R_i*alpha_i/2)*((1/norm_r2^2) - (1/norm_r1^2)) + beta_i*((1/norm_r2) - (1/norm_r1)));
    
    finclusion_total(1)=const_general*(f_total_1(1) + f_total_2(1) + f_total_3(1));
    
    finclusion_total(2)=-const_general*(b1/2)*((R_i*gamma_i/x^2)*(z_2/norm_r2 - z_1/norm_r1) + (delta_i/x)*(atan(z_2/x) - atan(z_1/x)));
    
    finclusion_total(3)=0;

    %%%%%%%%%%
    %
    %   This section calculates the force on node 2
    %
    %%%%%%%%%%
    
    f_2_1(1)=b2*(H_j^2/L_j)*(-(R_i*alpha_i/(2*x^3))*((x^3 + d*x*z_2)/norm_r2^2 - (x^3 + d*x*z_1)/norm_r1) + d*(atan(z_2/x) - atan(z_1/x)) + (beta_i/x^2)*((x^2 + d*z_2)/norm_r2 - (x^2 + d*z_1)/norm_r1));
    f_2_2(1)=(b2/2)*(1/L_j)*(-(R_i*gamma_i/x^2)*((x^2 + d*z_2)/norm_r2 - (x^2 + d*z_1)/norm_r1) + delta_i*((log(norm_r2^2/norm_r1^2)/2) - (d/x)*(atan(z_2/x) - atan(z_1/x))));
    f_2_3(1)=b3*(H_j/L_j)*((R_i*alpha_i/2)*((d - z_2)/norm_r2^2 - (d - z_1)/norm_r1 + (atan(z_2/x) - atan(z_1/x))/x) - beta_i*((d - z_2)/norm_r2 - (d - z_1)/norm_r1 + log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2)))));


    finclusion_2(1)=const_general*(f_2_1(1) + f_2_2(1) + f_2_3(1));
    
    finclusion_2(2)=-const_general*(b1/2)*(1/L_j)*(-(R_i*gamma_i/x^2)*((x^2 + d*z_2)/norm_r2 - (x^2 + d*z_1)/norm_r1) + delta_i*((log(norm_r2^2/norm_r1^2)/2) - (d/x)*(atan(z_2/x) - atan(z_1/x))));
    
    finclusion_2(3)=0;

    %%%%%%%%%%
    %
    %   This section calculates the force on node 1
    %
    %%%%%%%%%%

    finclusion_1=finclusion_total - finclusion_2;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
