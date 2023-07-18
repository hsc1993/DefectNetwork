%segments=constructsegmentlist(rn,links);
[LINKMAX,m]=size(links);
n1s1=links(1,1);
n2s1=links(1,2);
burg1=links(1,3:5);
n1s2=links(2,1);
n2s2=links(2,2);
burg2=links(2,3:5);
[dist2,ddist2dt,L1,L2]=mindistcalc(rn(n1s1,1:3),rn(n2s1,1:3),rn(n1s2,1:3),rn(n2s2,1:3));
mindist=sqrt(dist2);
for i=2:LINKMAX
    for j=i+1:LINKMAX
        n1s1=links(i,1);
        n2s1=links(i,2);
        burg1=links(i,3:5);
        n1s2=links(j,1);
        n2s2=links(j,2);
        burg2=links(j,3:5);
        if(n1s1~=n1s2)&(n1s1~=n2s2)&(n2s1~=n1s2)&(n2s1~=n2s2)
            if(burg1~=burg2)
                [dist2,ddist2dt,L1,L2]=mindistcalc(rn(n1s1,1:3),rn(n2s1,1:3),rn(n1s2,1:3),rn(n2s2,1:3));
                dist(i,j)=sqrt(dist2);
                dist(j,i)=sqrt(dist2);
                if(dist(i,j)<mindist)
                    mindist=dist(i,j);
                end
            end

        end
    end
end
disp(sprintf('the minimum distance is %f',mindist));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55


