function sf_force=StackingFaultForce(rn,rn0,doSFT,SFT_plane,sf_barrier,L)

SFTol=1e0;
dist_tol=0;
sf_force=zeros(1,3);

for i=0:(size(SFT_plane,1)/4)-1
    SFT_index=4*i+1;
    for j=0:3
        cut_point=CutLineSurface(rn0,SFT_plane(SFT_index+j,:));
        point_inside=PointInsideTriangle(cut_point,SFT_plane(SFT_index+j,:));
        dist=norm(rn0-cut_point);
        if((point_inside)&(dot((rn0-(SFT_plane(SFT_index+j,4:6))),SFT_plane(SFT_index+j,1:3))<=-dist_tol)&(dist<SFTol))
            sf_force=sf_force + sf_barrier*(SFT_plane(SFT_index+j,1:3)/norm(SFT_plane(SFT_index+j,1:3)))*(L/2);
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function cutpoint=CutLineSurface(rn0,SFT_plane)

n=SFT_plane(1:3)/norm(SFT_plane(1:3));
p=SFT_plane(4:6);

z=((n(1)^2+n(2)^2)*rn0(3) + n(3)^2*p(3) - n(1)*n(3)*(rn0(1)-p(1)) - n(2)*n(3)*(rn0(2)-p(2)))/(n(1)^2 + n(2)^2 + n(3)^2);
x=(n(1)/n(3))*(z - rn0(3)) + rn0(1);
y=(n(2)/n(3))*(z - rn0(3)) + rn0(2);

cutpoint=[x y z];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function point_inside=PointInsideTriangle(cut_point,SFT_plane)

tol=pi/100;

a=SFT_plane(4:6);
b=SFT_plane(7:9);
c=SFT_plane(10:12);

vec1=(a-cut_point)/norm(a-cut_point);
vec2=(b-cut_point)/norm(b-cut_point);
vec3=(c-cut_point)/norm(c-cut_point);

theta1=acos(dot(vec1,vec2));
theta2=acos(dot(vec2,vec3));
theta3=acos(dot(vec3,vec1));

if(((2*pi)-(theta1+theta2+theta3))<tol)
    point_inside=1;
else
    point_inside=0;
end
