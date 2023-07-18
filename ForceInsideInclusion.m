%This function gives the force on the nodes of a segment due to the elastic inclusion depending
%on where is the segment


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
    
    const_general=(1/5)*rho_i/(1-NU_i);
    %const_general=(1/28)*rho_i/(1-NU_i);
    
    
        %%%%%%%%%%
        %
        %   This section calculates the total force on the segment due to the inclusion
        %
        %%%%%%%%%%



    integral_node_1=zeros(1,3);
    integral_node_2=zeros(1,3);
    finclusion_total=zeros(1,3);
    
    f_total_1(1)=-b2*H_j^2*((alpha_i/R_i^3)*((z_2^3 -z_1^3)/3 + x^2*(z_2 - z_1)) - beta_i*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))));
    f_total_2(1)=-b2*((gamma_i/R_i^3)*((z_2^5 -z_1^5)/5 + (2*x^2/3)*(z_2^3 -z_1^3) + x^4*(z_2 - z_1)) - (delta_i/2)*(x^2*(log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2)))) + (z_2*norm_r2 - z_1*norm_r1)));
    f_total_3(1)=-b3*H_j*((alpha_i/(4*R_i^3))*((z_2^4 -z_1^4) + 2*x^2*(z_2^2 -z_1^2)) - beta_i*(norm_r2 - norm_r1));
    
    
    %f_total_1(1)=-(b2/5)*((H_j^2/R_i)*alpha_i*(z_2-z_1) - beta_i*H_j^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))));
    %f_total_2(1)=(b2/10)*((gamma_i/R_i)*((z_2^3-z_1^3)/3 + x^2*(z_2-z_1)) - (delta_i/2)*(x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))) + (z_2*norm_r2 - z_1*norm_r1)));
    %f_total_3(1)=-(b3/5)*((H_j/R_i)*(alpha_i/2)*(z_2^2-z_1^2) - beta_i*H_j*(norm_r2-norm_r1));
    
    
    finclusion_total(1)=const_general*(f_total_1(1) + f_total_2(1) + f_total_3(1));
    
    %f_total_1(2)=(b1/10)*((gamma_i/R_i)*((z_2^3-z_1^3)/3 + x^2*(z_2-z_1)) - (delta_i/2)*(x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))) + (z_2*norm_r2 - z_1*norm_r1)));
    f_total_1(2)=b1*((gamma_i/R_i^3)*((z_2^5 -z_1^5)/5 + (2*x^2/3)*(z_2^3 -z_1^3) + x^4*(z_2 - z_1)) - (delta_i/2)*(x^2*(log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2)))) + (z_2*norm_r2 - z_1*norm_r1)));
    
    finclusion_total(2)=const_general*f_total_1(2);
    
    finclusion_total(3)=0;

    %%%%%%%%%%
    %
    %   This section calculates the force on node 2
    %
    %%%%%%%%%%

    integral_node_1=zeros(1,3);
    integral_node_2=zeros(1,3);
    finclusion_2=zeros(1,3);
    
    %f_2_1(1)=-(d/L_j)*f_total_1(1) - b2*(H_j^2/L_j)*((alpha_i/(4*R_i^3))*((z_2^4 - z_1^4) + 2*x^2*(z_2^2 - z_1^2)) - beta_i*(norm_r2 - norm_r1));
    %f_2_2(1)=-(d/L_j)*f_total_2(1) - (b2/L_j)*((gamma_i/(6*R_i^3))*((z_2^6 - z_1^6) + 3*x^2*(z_2^4 - z_1^4) + 3*x^4*(z_2^2 -z_1^2)) - (delta_i/3)*(norm_r2^3 - norm_r1^3));
    %f_2_3(1)=-(d/L_j)*f_total_3(1) - b3*(H_j/L_j)*((alpha_i/R_i^3)*((z_2^5 - z_1^5)/5 + (x^2/3)*(z_2^3 - z_1^3)) - (beta_i/2)*((z_2*norm_r2 - z_1*norm_r1) - x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2)))));
    
    
    %f_2_1(1)=-(b2/5)*(H_j^2/L_j)*((alpha_i/R_i)*((z_2^2-z_1^2)/2 - d*(z_2-z_1)) - beta_i*((norm_r2-norm_r1) - d*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2)))));
    %f_2_2(1)=(b2/10)*((gamma_i/(R_i*L_j)*((z_2^4-z_1^4)/4 - d*(z_2^3-z_1^3)/3 + x^2*(z_2^2-z_1^2)/2 - d*x^2*(z_2-z_1))) - (delta_i/(6*L_j)*(norm_r2*(2*norm_r2^2 - 3*d*z_2) - norm_r1*(2*norm_r1^2 - 3*d*z_1) - 3*d*x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))))));
    %f_2_3(1)=-(b3/5)*(H_j/L_j)*((alpha_i/R_i)*((z_2^3-z_1^3)/3 - d*(z_2^2-z_1^2)/2) - beta_i*(norm_r2*(z_2/2 - d) - norm_r1*(z_1/2 - d) - (x^2/2)*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2)))));
    
    %finclusion_2(1)=const_general*(f_2_1(1) + f_2_2(1) + f_2_3(1));
    
    %f_2_1(2)=(b1/10)*((gamma_i/(R_i*L_j)*((z_2^4-z_1^4)/4 - d*(z_2^3-z_1^3)/3 + x^2*(z_2^2-z_1^2)/2 - d*x^2*(z_2-z_1))) - (delta_i/(6*L_j)*(norm_r2*(2*norm_r2^2 - 3*d*z_2) - norm_r1*(2*norm_r1^2 - 3*d*z_1) - 3*d*x^2*log((z_2+sqrt(x^2+z_2^2))/(z_1+sqrt(x^2+z_1^2))))));
    %f_2_1(2)=-(d/L_j)*f_total_1(2) + (b1/L_j)*((gamma_i/(6*R_i^3))*((z_2^6 - z_1^6) + 3*x^2*(z_2^4 - z_1^4) + 3*x^4*(z_2^2 -z_1^2)) - (delta_i/3)*(norm_r2^3 - norm_r1^3));
    
    %finclusion_2(2)=const_general*f_2_1(2);

    %finclusion_2(3)=0;

    finclusion_2=finclusion_total/2;
    %%%%%%%%%%
    %
    %   This section calculates the force on node 1
    %
    %%%%%%%%%%

    %finclusion_1=finclusion_total - finclusion_2;
    finclusion_1=finclusion_total/2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


