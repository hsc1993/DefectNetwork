% C11 = 1.18;
% C12 = 0.623;
% C44 = 0.325;

C11 = 169.9;
C12 = 122.6;
C44 = 76.2;

C = zeros(3,3,3,3);
% C(1,1,1,1)=C11; %C(11)
% C(2,2,2,2)=C11; %C(22)
% C(3,3,3,3)=C11; %C(33)
% C(1,1,2,2)=C12; %C(12)
% C(2,2,1,1)=C12; %C(21)
% C(3,3,2,2)=C12; %C(32)
% C(2,2,3,3)=C12; %C(23)
% C(1,1,3,3)=C12; %C(13)
% C(3,3,1,1)=C12; %C(31)
% C(1,2,1,2)=C44; %C(44)
% C(1,2,3,2)=C44; %C(47)
% C(3,2,1,2)=C44; %C(74)
% C(1,3,1,3)=C44; %C(55)
% C(2,1,3,1)=C44; %C(58)
% C(3,1,2,1)=C44; %C(85)
% C(2,3,2,3)=C44; %C(66)
% C(2,3,1,3)=C44; %C(69)
% C(1,3,2,3)=C44; %C(96)
% C(2,1,2,1)=C44; %C(77)
% C(3,1,3,1)=C44; %C(88)
% C(3,2,3,2)=C44; %C(99)
C(1,1,1,1) = C11; C(2,2,2,2) = C11; C(3,3,3,3) = C11;
C(1,1,2,2) = C12; C(1,1,3,3) = C12; C(2,2,1,1) = C12;
C(2,2,3,3) = C12; C(3,3,1,1) = C12; C(3,3,2,2) = C12;
C(1,2,1,2) = C44; C(2,1,2,1) = C44; C(1,3,1,3) = C44;
C(3,1,3,1) = C44; C(2,3,2,3) = C44; C(3,2,3,2) = C44;
C(1,2,2,1) = C44; C(2,1,1,2) = C44; C(1,3,3,1) = C44;
C(3,1,1,3) = C44; C(2,3,3,2) = C44; C(3,2,2,3) = C44;

T = [1/sqrt(6) 1/sqrt(6) -2/sqrt(6);
    1/sqrt(3) 1/sqrt(3) 1/sqrt(3);
    1/sqrt(2) -1/sqrt(2) 0];

e1p = [1 0 0]; e2p = [0 1 0]; e3p = [0 0 1];
%e1p=e1p/norm(e1p); e2p=e2p/norm(e2p); e3p=e3p/norm(e3p);            
%coordinate system 2
e1 =  [1 1 -2]; e2  = [1 1 1]; e3 = [-1 1 0];
e1=e1/norm(e1); e2=e2/norm(e2); e3=e3/norm(e3);

%rotation matrix
Q = [ dot(e1,e1p) dot(e1,e2p) dot(e1,e3p)
      dot(e2,e1p) dot(e2,e2p) dot(e2,e3p)
      dot(e3,e1p) dot(e3,e2p) dot(e3,e3p) ];

Cprime =zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for g=1:3
                    for h=1:3
                        for m=1:3
                            for n =1:3
                                Cprime(i,j,k,l)=Cprime(i,j,k,l)+T(i,g)*T(j,h)*T(k,m)*T(l,n)*C(g,h,m,n);
                            end
                        end
                    end
                end
            end
        end
    end
end
disp(sprintf('Cp(1,1,1,1) = %g',Cprime(1,1,1,1)));
disp(sprintf('Cp(2,2,2,2) = %g',Cprime(2,2,2,2)));
disp(sprintf('Cp(3,3,3,3) = %g',Cprime(3,3,3,3)));

disp(sprintf('Cp(1,1,2,2) = %g',Cprime(1,1,2,2)));
disp(sprintf('Cp(1,1,3,3) = %g',Cprime(1,1,3,3)));
disp(sprintf('Cp(2,2,3,3) = %g',Cprime(2,2,3,3)));

disp(sprintf('Cp(2,3,2,3) = %g',Cprime(2,3,2,3)));
disp(sprintf('Cp(3,1,3,1) = %g',Cprime(3,1,3,1)));
disp(sprintf('Cp(1,2,1,2) = %g',Cprime(1,3,1,3)));
% CC=zeros(9,9);
% CC(1,1)=C11; %C(11)
% CC(2,2)=C11; %C(22)
% CC(3,3)=C11; %C(33)
% CC(1,2)=C12; %C(12)
% CC(2,1)=C12; %C(21)
% CC(3,2)=C12; %C(32)
% CC(2,3)=C12; %C(23)
% CC(1,3)=C12; %C(13)
% CC(3,1)=C12; %C(31)
% CC(4,4)=C44; %C(44)
% CC(4,7)=C44; %C(47)
% CC(7,4)=C44; %C(74)
% CC(5,5)=C44; %C(55)
% CC(5,8)=C44; %C(58)
% CC(8,5)=C44; %C(85)
% CC(6,6)=C44; %C(66)
% CC(6,9)=C44; %C(69)
% CC(9,6)=C44; %C(96)
% CC(7,7)=C44; %C(77)
% CC(8,8)=C44; %C(88)
% CC(9,9)=C44; %C(99)
% A = zeros(3,3,3,3);
% for a=1:3
%     for b=1:3
%         for c=1:3
%             for d=1:3
%                 A(a,b,c,d) = CC(3*(a-1)+b,3*(c-1)+d);
%             end
%         end
%     end
% end
% Aprime =zeros(3,3,3,3);
% for i=1:3
%     for j=1:3
%         for k=1:3
%             for l=1:3
%                 for g=1:3
%                     for h=1:3
%                         for m=1:3
%                             for n =1:3
%                                 Aprime(i,j,k,l)=Aprime(i,j,k,l)+T(i,g)*T(j,h)*A(g,h,m,n)*T(k,m)*T(l,n);
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% CCprime =zeros(9,9);
% for a=1:3
%     for b=1:3
%         for c=1:3
%             for d=1:3
%                 CCprime(3*(a-1)+b,3*(c-1)+d) = Aprime(a,b,c,d);
%             end
%         end
%     end
% end
% CCprime