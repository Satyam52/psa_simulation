function [Ai,Ax,Bx,z] = Get_matrix(M,Pe)

z = [0  0.0469 0.2308 0.5000 0.7692 0.9531 1];


A = zeros(7,7);
D = zeros(7,7);
C= zeros(7,7);

for j = 1:7
    A(j,:) = [1,z(j),z(j)^2,z(j)^3,z(j)^4,z(j)^5,z(j)^6];
    C(j,:) = [0,1,2*z(j),3*z(j)^2,4*z(j)^3,5*z(j)^4,6*z(j)^5];
    D(j,:) = [0,0,2,6*z(j),12*z(j)^2,20*z(j)^3,30*z(j)^4];
    
end

Ainv = inv(A);
Ax = C*Ainv;
Bx = D*Ainv;


% Kent Knaebel , PSA Textbook
A3=(Ax(1,1)-Pe)/Ax(M+2,1);
A4=1/Ax(M+2,1);
A1=1/(Ax(1,M+2)-(A3*Ax(M+2,M+2)));
A2=Ax(M+2,M+2)*A1;
A5=A2*A4;
Ai=[A1,A2,A3,A4,A5];

end