function [rn,links] = generateHalfPDL(b,numSegs,shift)

rn_loop = [];
b_unitvec = b/norm(b)*numSegs;

n1 = 1/sqrt(2)*[ 0 -1 1 ];
n1 = 1/sqrt(2)*[ 1 -1 0 ];
t1 = cross(n1,b_unitvec);

% % make number of 
% t1 = t1*numSegs;
numSegs = 1;

links1 = [];
for i = 1:numSegs
    rn1(i,:) = [i*t1(1)+shift(1) i*t1(2)+shift(2) i*t1(3)+shift(3) 0];
    links1(i,:) = [i i+1 b n1];
end
links1 = links1(1:end-1,:);

n2_start = rn1(end,:);
n2 = 1/sqrt(2)*[1 -1 0];
n2 = 1/sqrt(2)*[0 -1 1];
t2 = cross(n2,b_unitvec);
links2 = [];
for i = 1:numSegs
    rn2(i,:) = n2_start+[i*t2(1) i*t2(2) i*t2(3) 0];
    links2(i,:) = [i i+1 b n2];
end
links2(:,1:2) = links2(:,1:2)+size(links1,1);

n3_start = rn2(end,:);
n3 = 1/sqrt(2)*[1 0 -1];
n3 = 1/sqrt(2)*[-1 0 1];
t3 = cross(n3,b_unitvec);
links3 = [];
for i = 1:numSegs
    rn3(i,:) = n3_start+[i*t3(1) i*t3(2) i*t3(3) 0];
    links3(i,:) = [i i+1 b n3];
end
links3(:,1:2) = links3(:,1:2)+size(links1,1)+size(links2,1);

n4_start = rn3(end,:);
n4 = 1/sqrt(2)*[0 1 -1];
t4 = cross(n4,b_unitvec);
links4 = [];
for i = 1:numSegs
    rn4(i,:) = n4_start+[i*t4(1) i*t4(2) i*t4(3) 0];
    links4(i,:) = [i i+1 b n4];
end
% 
% links4(:,1:2) = links4(:,1:2)+size(links1,1)+size(links2,1)+size(links3,1);
% 
% n5_start = rn4(end,:);
% n5 = 1/sqrt(2)*[-1 1 0];
% t5 = cross(n5,b_unitvec);
% links5 = [];
% for i = 1:numSegs
%     rn5(i,:) = n5_start+[i*t5(1) i*t5(2) i*t5(3) 0];
%     links5(i,:) = [i i+1 b n5];
% end
% links5(:,1:2) = links5(:,1:2)+size(links1,1)+size(links2,1)+size(links3,1)+size(links4,1);
% 
% n6_start = rn5(end,:);
% n6 = 1/sqrt(2)*[-1 0 1];
% t6 = cross(n6,b_unitvec);
% links6 = [];
% for i = 1:numSegs
%     rn6(i,:) = n6_start+[i*t6(1) i*t6(2) i*t6(3) 0];
%     links6(i,:) = [i i+1 b n6];
% end
% links6(:,1:2) = links6(:,1:2)+size(links1,1)+size(links2,1)+size(links3,1)+size(links4,1)+size(links5,1);
% 
% rn = [rn1;rn2;rn3;rn4;rn5;rn6];
% links = [links1;links2;links3;links4;links5;links6];
% links = [links;[size(rn,1) 1 b n1]];

rn = [rn1;rn2;rn3;];
links = [links1;links2;links3;];
links = [links;[size(rn,1) 1 b n1]];

end
