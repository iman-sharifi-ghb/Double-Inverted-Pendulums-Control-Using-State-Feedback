function [A,B,C,D]=State_Space()
    global M m1 m2 L1 L2 g d1 d2 d3 w1 w2 w3 
    syms u1 u2
    syms q1 q2 q3
    syms dq1 dq2 dq3

    M=[M+m1+m2               L1*(m1+m2)*cos(q2)    m2*L2*cos(q3);
       L1*(m1+m2)*cos(q2)    L1^2*(m1+m2)         L1*L2*m2*cos(q2-q3);
       L2*m2*cos(q3)         L1*L2*m2*cos(q2-q3)  L2^2*m2];
    M = subs(M,[q1, q2, q3],[0 0 0]);

    B1 =[L1*(m1+m2)*(dq2^2)*sin(q2)+m2*L2*(dq2^2);
        -L1*L2*m2*(dq2^2)*sin(q2-q3)+g*(m1+m2)*L1*sin(q2);
        L1*L2*m2*(dq2^2)*sin(q2-q3)+g*L2*m2*sin(q3)];

    B2=[d1*dq1;d2*dq2;d3*dq3];
    B3=[u1;u2;0];
    B4=[w1;w2;w3];

    B = B1-B2+B3+B4;
    ddq = simplify(M\B);

    f1 = dq1;
    f2 = ddq(1);
    f3 = dq2;
    f4 = ddq(2);
    f5 = dq3;
    f6 = ddq(3);

    q = [q1;dq1;q2;dq2;q3;dq3];
    f = [f1;f2;f3;f4;f5;f6];
    u = [u1;u2];

    A = jacobian(f,q);% A = simplify(A);% consuming time
    B = jacobian(f,u);% B = simplify(B);

    A = subs(A,[q1,q2,q3,dq1,dq2,dq3,u1,u2],[0,0,0,0,0,0,0,0]);
    B = subs(B,[q1,q2,q3,dq1,dq2,dq3,u1,u2],[0,0,0,0,0,0,0,0]);

    A = vpa(A,6);B = vpa(B,6);
    A = double(A);B = double(B);
    C = [1 0 0 0 0 0;0 0 1 0 0 0];D = zeros(2,2);
end

