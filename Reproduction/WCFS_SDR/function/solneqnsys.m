%% this function is to solve the following two equations with 3 vars.
%%                  a1*s1*s2+a2*s2*s3+a3*s1*s3        =0
%% bm1*s1^2+b0*s2^2+b1*s1*s2+b2*s2*s3+b3*s1*s3+b4*s3^2=0,
%% where bm1>0 and b0<0;

%% input: a -- 9-dim vector (real number)

%% output: s1,s2,s3: a solution of the equation system
%%                   (as long as bm1*b0<0, there must be a NONZERO solution)


function s=solneqnsys(a)

epsi=10^(-8);
a1=a(1);a2=a(2);a3=a(3);
bm1=a(4);b0=a(5);b1=a(6);b2=a(7);b3=a(8);b4=a(9);

s=zeros(3,1);
s1=s(1);s2=s(2);s3=s(3);

if (length(a)~=9)||(bm1*b0>=0)
    disp('coefficents assumption not satisfied - from solneqnsys')
    return
end

if bm1<0
    bm1=-bm1;b0=-b0;b1=-b1;b2=-b2;b3=-b3;b4=-b4;
end

if abs(a1)<=epsi
    s3=0;s2=1;s1=(-b1+sqrt(b1^2-4*bm1*b0))/(2*bm1);
    return
else 
    s3=1;
end

if (abs(a2)<=epsi)&&(abs(a3)<=epsi)
    if b4<=0
        s2=0;s1=(-b3+sqrt(b3^2-4*bm1*b4))/(2*bm1);
    else
        s1=0;s2=(-b2+sqrt(b2^2-4*b0*b4))/(2*b0);
    end
elseif (abs(a2)>epsi)&&(abs(a3)<=epsi)
    if b4<=(b3^2)/(4*bm1)
        s2=0;s1=(-b3+sqrt(b3^2-4*bm1*b4))/(2*bm1);
    else
        s1=-a2/a1;
        c1=b2+b1*s1;
        c2=bm1*(s1)^2+b3*s1+b4;
        s2=(-c1+sqrt(c1^2-4*b0*c2))/(2*b0);
    end
elseif (abs(a2)<=epsi)&&(abs(a3)>epsi)
    if b4>=(b2^2)/(4*b0)
        s1=0;s2=(-b2+sqrt(b2^2-4*b0*b4))/(2*b0);
    else
        s2=-a3/a1;
        c1=b3+b1*s2;
        c2=b0*(s2)^2+b2*s2+b4;
        s1=(-c1+sqrt(c1^2-4*bm1*c2))/(2*bm1);
    end
elseif (abs(a2)>epsi)&&(abs(a3)>epsi)
    lam1=(a1^2)*b0;
    lam2=2*a1*a3*b0-a1*a2*b1+(a1^2)*b2;
    lam3=(a2^2)*bm1+(a3^2)*b0-a2*a3*b1+2*a1*a3*b2-a1*a2*b3+(a1^2)*b4;
    lam4=(a3^2)*b2-a2*a3*b3+2*a1*a3*b4;
    lam5=b4*a3^2;
    rs=roots([lam1,lam2,lam3,lam4,lam5]);
    for ii=1:length(rs)
        if (abs(imag(rs(ii)))<=epsi)&&(abs(a3+a1*real(rs(ii)))>epsi)
            s2=real(rs(ii));
            s1=-a2*s2/(a3+a1*s2);
            break;
        end
    end
end
        
s(1)=s1;s(2)=s2;s(3)=s3;
% Test:
% a=randn(9,1);a(4)=abs(a(4));a(5)=-abs(a(5));
% s=solneqnsys(a);
%eqn1=a(1)*s1*s2+a(2)*s2*s3+a(3)*s1*s3
%eqn2=a(4)*s1^2+a(5)*s2^2+a(6)*s1*s2+a(7)*s2*s3+a(8)*s1*s3+a(9)*s3^2


