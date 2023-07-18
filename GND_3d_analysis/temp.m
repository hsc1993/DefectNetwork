for m=dt:dt:t

x1(end+1)=x1(end)+v1(end)*dt+(f(end)/m1)*(dt*dt/2);
 x2(end+1)=x2(end)+v2(end)*dt-
(f(end)/m2)*(dt*dt/2);
 l(end+1)=x2(end)-x1(end);
 f(end+1)=K*(l(end)-l0);
 v1(end+1)=(x1(end)-x1(end1))/dt+(f(end)*dt/(2*m1));
 v2(end+1)=(x2(end)-x2(end-1))/dt-
(f(end)*dt/(2*m2));
 T(end+1)=T(end)+dt;
 end
xr1=[x1(end)];
xr2=[x2(end)];
vr1=[v1(end)];
vr2=[v2(end)];
lr=[l(end)];
fr=[K*(l(end)-l0)];
Tr=[2];
for m=dt:dt:t
 xr1(end+1)=xr1(end)+vr1(end)*-dt+(fr(end)/m1)*(-
dt*-dt/2);
 xr2(end+1)=xr2(end)+vr2(end)*-dt-(fr(end)/m2)*(-
dt*-dt/2);
 lr(end+1)=xr2(end)-xr1(end);
 fr(end+1)=K*(lr(end)-l0);
 vr1(end+1)=(xr1(end)-xr1(end-1))/-dt+(fr(end)*-
dt/(2*m1));
 vr2(end+1)=(xr2(end)-xr2(end-1))/-dt-(fr(end)*-
dt/(2*m2));
 Tr(end+1)=Tr(end)-dt;
end
plot(Tr,xr1,Tr,xr2);
xlabel('t');
ylabel('x');
legend('m1','m2')