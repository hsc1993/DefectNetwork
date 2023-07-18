%%%%%%%%%
%  Try functions
%%%%%%%%%

x=linspace(0,100,1000);
C=10;
for i=1:1000
    %f(i)=x(i)^3/(x(i)^2+C^2)^(3/2);
    f(i)=1/(x(i)^2+C^2)^(3/2);
    %f(i)=-x(i)/(x(i)^2+C^2)^(1/2) + atan(x(i)/C);
    %f(i)=atan(x(i));
end
plot(x,f,'-r');