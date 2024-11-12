%% this function is to solve the following two equations with 3 vars.
%%                      real(a1*conj(s1)*s2)+real(a2*conj(s2)*s3)+real(a3*conj(s3)*s1)         =0
%%                      real(b1*conj(s1)*s2)+real(b2*conj(s2)*s3)+real(b3*conj(s3)*s1)         =0
%% cm1*|s1|^2+c0*|s2|^2+real(c1*conj(s1)*s2)+real(c2*conj(s2)*s3)+real(c3*conj(s3)*s1)+c4*|z|^2=0,
%% where cm1>0, c0<0, c4 is real,

%% In the above, a1, a2, a3, b1, b2, b3, c1, c2, c3 are complex numbers; they together with cm1, c0,
%% are all inputs.


%% input: a -- 12-dim vector (9 complex numbers, and 3 real number)
%%             a1=a(1);a2=a(2);a3=a(3);
%%             b1=a(4);b2=a(5);b3=a(6);
%%             cm1=a(7);  % real number
%%             c0=a(8);   % real number
%%             c1=a(9);c2=a(10);c3=a(11);
%%             c4=a(12);  % real number


%% output: s1,s2,s3: a complex-valued solution of the equation system
%%                  (as long as cm1*c0<0, there must be a NONZERO
%%                  solution)

%% To run this function, you need: solneqn3.m
%%                                 solneqnsys.m


function s=solneqnsys3(a)

epsi=10^(-8);
a1=a(1);a2=a(2);a3=a(3);
b1=a(4);b2=a(5);b3=a(6);
cm1=real(a(7));c0=real(a(8));c1=a(9);c2=a(10);c3=a(11);c4=real(a(12));

s=zeros(3,1);
s1=s(1);s2=s(2);s3=s(3);

if (length(a)~=12)||(cm1*c0>=0)||(abs(imag(a(12)))>epsi)||(abs(imag(a(7)))>epsi)||(abs(imag(a(8)))>epsi)
    disp('coefficents assumption not satisfied - Stop 1 from solneqnsys3')
    return
end

if cm1<0
    cm1=-cm1;c0=-c0;c1=-c1;c2=-c2;c3=-c3;c4=-c4;
end

Delta=real(a1)*(-imag(b1))-real(b1)*(-imag(a1));
if abs(Delta)<epsi
    s3=0;
else
    s3=1;
end

if abs(Delta)<epsi
    if (abs(imag(a1))<epsi)&&(abs(imag(b1))<epsi)
        x=imag(c1)/sqrt(-c0*cm1);y=0;z=1;
        r1p=sqrt((x+sqrt(x^2+4*(y^2+z^2)))/2);
        r2p=sqrt((-x+sqrt(x^2+4*(y^2+z^2)))/2);
        tt2=pi/2;tt1=0;
        r1=r1p/sqrt(cm1);r2=r2p/sqrt(-c0);
        s1=r1*exp(i*tt1);
        s2=r2*exp(i*tt2);
    elseif abs(imag(a1))>=epsi
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
    Aa=[-imag(a1p),real(a1p),0;-imag(b1p),real(b1p),0;-imag(c1p),real(c1p),1];
    app=inv(Aa)*[a2p,a3p,0;b2p,b3p,0;c2p,c3p,c4];
    na=zeros(7,1);
    na(1)=app(1,1);na(2)=app(1,2);na(3)=app(2,1);na(4)=app(2,2);
    na(5)=app(3,1);na(6)=app(3,2);na(7)=real(app(3,3));
    ss=solneqn3(na);
    s1=ss(1)/sqrt(cm1);
    s2=ss(2)/sqrt(-c0);
end

s(1)=s1;s(2)=s2;s(3)=s3;



% eqn1=real(a1*conj(s1)*s2)+real(a2*conj(s2)*s3)+real(a3*conj(s3)*s1)
% eqn2=real(b1*conj(s1)*s2)+real(b2*conj(s2)*s3)+real(b3*conj(s3)*s1)
% eqn3=cm1*abs(s1)^2+c0*abs(s2)^2+real(c1*conj(s1)*s2)+real(c2*conj(s2)*s3)+real(c3*conj(s3)*s1)+c4*(abs(s3))^2

%Test: (s3=1)
%ar=randn(12,1);ai=randn(12,1);a=complex(ar,ai);a(7)=abs(a(7));a(8)=-abs(a(8));a(12)=real(a(12));
%s=solneqnsys3(a);
% or Test (s3=0)
%ar=randn(12,1);ai=randn(12,1);a=complex(ar,ai);a(7)=abs(a(7));a(8)=-abs(a(8));a(12)=real(a(12));
%d=randn;d2=randn;a(1)=d+d*i;a(4)=1/2*d2+i/2*d2;s=solneqnsys3(a)