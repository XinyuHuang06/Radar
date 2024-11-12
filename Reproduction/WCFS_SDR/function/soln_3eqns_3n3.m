
%% this function is to solve the following two equations with 3 vars.
%%                      real(a1*conj(s1)*s2)+real(a2*conj(s2)*s3)+real(a3*conj(s3)*s1)+a4*|z|^2 =0
%%                      real(b1*conj(s1)*s2)+real(b2*conj(s2)*s3)+real(b3*conj(s3)*s1)+b4*|z|^2 =0
%% cm1*|s1|^2+c0*|s2|^2+real(c1*conj(s1)*s2)+real(c2*conj(s2)*s3)+real(c3*conj(s3)*s1)+c4*|z|^2 =0,
%% where cm1>0, c0<0, a4,b4,c4 are real,

%% In the above, a1, a2, a3, b1, b2, b3, c1, c2, c3 are complex numbers;
%% they, together with cm1, c0, a4, b4, c4, are all inputs.


%% input: a -- 14-dim vector (9 complex numbers, and 5 real number)
%%        a1=a(1);a2=a(2);a3=a(3);
%%        a4=a(4);  % real number
%%        b1=a(5);b2=a(6);b3=a(7);
%%        b4=a(8);  % real number
%%        cm1=a(9);   % real number
%%        c0=a(10);   % real number
%%        c1=a(11);c2=a(12);c3=a(13);
%%        c4=a(14);   % real number



%% output: s1,s2,s3: a complex-valued solution of the equation system
%%                  (as long as cm1*c0<0, there must be a NONZERO
%%                  solution)

%% To run this function, you need: soln_1eqn_3n.m
%%                                soln_2eqns_3n.m
%%                                     solneqn3.m
%%                                   solneqnsys.m
%%                               soln_3eqns_3n2.m

function s=soln_3eqns_3n3(a)

epsi=10^(-8);
epsii=10^(-8);

a1=a(1);a2=a(2);a3=a(3);a4=real(a(4));
b1=a(5);b2=a(6);b3=a(7);b4=real(a(8));
cm1=real(a(9));c0=real(a(10));c1=a(11);c2=a(12);c3=a(13);c4=real(a(14));

s=zeros(3,1);
s1=s(1);s2=s(2);s3=s(3);

%% check in input
if (length(a)~=14)||(cm1*c0>=0)||(abs(imag(a(4)))>epsii)||(abs(imag(a(8)))>epsii)||(abs(imag(a(14)))>epsii)||(abs(imag(a(9)))>epsii)||(abs(imag(a(10)))>epsii)
    disp('coefficents assumption not satisfied - Stop 1 from soln_3eqns_3n3')
    return
end
%% end of check in input

if cm1<0
    cm1=-cm1;c0=-c0;c1=-c1;c2=-c2;c3=-c3;c4=-c4;
end
%% coefficient sign adjust.

Delta=real(a1)*(-imag(b1))-real(b1)*(-imag(a1));
if abs(Delta)<epsii
    s3=0;
else
    s3=1;
end

if abs(Delta)<epsii
    if (abs(imag(a1))<epsii)&&(abs(imag(b1))<epsii)
        x=imag(c1)/sqrt(-c0*cm1);y=0;z=1;
        r1p=sqrt((x+sqrt(x^2+4*(y^2+z^2)))/2);
        r2p=sqrt((-x+sqrt(x^2+4*(y^2+z^2)))/2);
        tt2=pi/2;tt1=0;
        r1=r1p/sqrt(cm1);r2=r2p/sqrt(-c0);
        s1=r1*exp(i*tt1);
        s2=r2*exp(i*tt2);
    elseif abs(imag(a1))>=epsii
        x=-real(c1)/sqrt(-c0*cm1)+imag(c1)*real(a1)/(sqrt(-c0*cm1)*imag(a1));
        y=1;z=real(a1)/imag(a1);
        r1p=sqrt((x+sqrt(x^2+4*(y^2+z^2)))/2);
        r2p=sqrt((-x+sqrt(x^2+4*(y^2+z^2)))/2);
        tt1=0;
        if z<=0
            tt2=-acos(y/sqrt(y^2+z^2));
        else
            tt2=acos(y/sqrt(y^2+z^2));
        end
        r1=r1p/sqrt(cm1);r2=r2p/sqrt(-c0);
        s1=r1*exp(i*tt1);
        s2=r2*exp(i*tt2);
    else %% this branch imag(b1)~=0, imag(a1)==0
        x=-real(c1)/sqrt(-c0*cm1)+imag(c1)*real(b1)/(sqrt(-c0*cm1)*imag(b1)); 
        y=1;z=real(b1)/imag(b1);
        r1p=sqrt((x+sqrt(x^2+4*(y^2+z^2)))/2);
        r2p=sqrt((-x+sqrt(x^2+4*(y^2+z^2)))/2);
        tt1=0;
        if z<=0
            tt2=-acos(y/sqrt(y^2+z^2));
        else
            tt2=acos(y/sqrt(y^2+z^2));
        end
        r1=r1p/sqrt(cm1);r2=r2p/sqrt(-c0);
        s1=r1*exp(i*tt1);
        s2=r2*exp(i*tt2);
    end
else
    a1p=a1/sqrt(-c0*cm1);b1p=b1/sqrt(-c0*cm1);c1p=c1/sqrt(-c0*cm1);
    a2p=a2/sqrt(-c0);b2p=b2/sqrt(-c0);c2p=c2/sqrt(-c0);
    a3p=a3/sqrt(cm1);b3p=b3/sqrt(cm1);c3p=c3/sqrt(cm1);
    a4p=a4;b4p=b4;c4p=c4; %a4p,b4p,c4p real numbers
    Aa=[-imag(a1p),real(a1p),0;-imag(b1p),real(b1p),0;-imag(c1p),real(c1p),1]; 
    app=inv(Aa)*[a2p,a3p,a4p;b2p,b3p,b4p;c2p,c3p,c4p]; %app(1,3),app(2,3),app(3,3) real numbers
    na=zeros(9,1);
    na(1)=app(1,1);na(2)=app(1,2);na(3)=real(app(1,3));
    na(4)=app(2,1);na(5)=app(2,2);na(6)=real(app(2,3));
    na(7)=app(3,1);na(8)=app(3,2);na(9)=real(app(3,3));
    ss=soln_3eqns_3n2(na);
    s1=ss(1)/sqrt(cm1);
    s2=ss(2)/sqrt(-c0);
end

s(1)=s1;s(2)=s2;s(3)=s3;

eqn1=real(a1*conj(s1)*s2)+real(a2*conj(s2)*s3)+real(a3*conj(s3)*s1)+a4*(abs(s3))^2;
eqn2=real(b1*conj(s1)*s2)+real(b2*conj(s2)*s3)+real(b3*conj(s3)*s1)+b4*(abs(s3))^2;
eqn3=cm1*abs(s1)^2+c0*abs(s2)^2+real(c1*conj(s1)*s2)+real(c2*conj(s2)*s3)+real(c3*conj(s3)*s1)+c4*(abs(s3))^2;
  
%% Test (s3=1)
%% ia=randn(14,1);ib=randn(14,1);a=complex(ia,ib);a(4)=real(a(4));a(8)=real(a(8));a(9)=abs(a(9));a(10)=-abs(a(10));a(14)=real(a(14));
%% soln_3eqns_3n3(a)
% an extreme data: a=[0.0140-0.0063i 1.2129+0.0644i 1.4109-0.1683i 0.0865 2.5453-0.0280i -0.0858-1.5420i 1.2868-0.1584i 0.2477 0.4728 -2.0605 0.2871-1.1715i -0.2692-1.2130i -0.1822+2.0161i 0.3272] 

%% Test (s3=0)
%% ia=randn(14,1);ib=randn(14,1);a=complex(ia,ib);a(4)=real(a(4));a(8)=real(a(8));a(9)=abs(a(9));a(10)=-abs(a(10));a(14)=real(a(14));
%% d=randn;d2=randn;a(1)=d+d*i;a(5)=1/2*d2+i/2*d2; soln_3eqns_3n3(a)
