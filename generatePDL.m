function [rn,links] = generatePDL(nlist,b,numSegs,shift,id,seglength_multiply)

rn_loop = [];
% b_unitvec = b/norm(b);

n1 = nlist(1,:);
n2 = nlist(2,:);
n3 = nlist(3,:);
n4 = nlist(4,:);
n5 = nlist(5,:);
n6 = nlist(6,:);

% use seglength_multiplier to increase length of segments, while preventing
% discretization
if seglength_multiply == 1
    seglength_multiplier = numSegs;
    numSegs = 1;
else
    seglength_multiplier = 1;
end

t1 = cross(b,n1)*seglength_multiplier;
links1 = [];
for i = 1:numSegs
    rn1(i,:) = [i*t1(1)+shift(1) i*t1(2)+shift(2) i*t1(3)+shift(3) 0 id];
    links1(i,:) = [i i+1 b n1];
end
% links1 = links1(1:end-1,:);


n2_start = rn1(end,:);
t2 = cross(b,n2)*seglength_multiplier;
links2 = [];
for i = 1:numSegs
    rn2(i,:) = n2_start+[i*t2(1) i*t2(2) i*t2(3) 0 0];
    links2(i,:) = [i i+1 b n2];
end
links2(:,1:2) = links2(:,1:2)+size(links1,1);

n3_start = rn2(end,:);
t3 = cross(b,n3)*seglength_multiplier;
links3 = [];
for i = 1:numSegs
    rn3(i,:) = n3_start+[i*t3(1) i*t3(2) i*t3(3) 0 0];
    links3(i,:) = [i i+1 b n3];
end
links3(:,1:2) = links3(:,1:2)+size(links1,1)+size(links2,1);

n4_start = rn3(end,:);
t4 = cross(b,n4)*seglength_multiplier;
links4 = [];
for i = 1:numSegs
    rn4(i,:) = n4_start+[i*t4(1) i*t4(2) i*t4(3) 0 0];
    links4(i,:) = [i i+1 b n4];
end
links4(:,1:2) = links4(:,1:2)+size(links1,1)+size(links2,1)+size(links3,1);

n5_start = rn4(end,:);
t5 = cross(b,n5)*seglength_multiplier;
links5 = [];
for i = 1:numSegs
    rn5(i,:) = n5_start+[i*t5(1) i*t5(2) i*t5(3) 0 0];
    links5(i,:) = [i i+1 b n5];
end
links5(:,1:2) = links5(:,1:2)+size(links1,1)+size(links2,1)+size(links3,1)+size(links4,1);

n6_start = rn5(end,:);
t6 = cross(b,n6)*seglength_multiplier;
links6 = [];
for i = 1:numSegs
    rn6(i,:) = n6_start+[i*t6(1) i*t6(2) i*t6(3) 0 0];
    links6(i,:) = [i i+1 b n6];
end
links6 = links6(1:end-1,:);
links6(:,1:2) = links6(:,1:2)+size(links1,1)+size(links2,1)+size(links3,1)+size(links4,1)+size(links5,1);

rn = [rn1;rn2;rn3;rn4;rn5;rn6];
% rn = [rn1;rn2;rn3;rn4;rn5;];
links = [links1;links2;links3;links4;links5;links6];
% links = [links1;links2;links3;links4;links5;];
% links = [links;[size(rn,1) 1 b n1]];
links = [links;[size(rn,1) 1 b n6]];

end
