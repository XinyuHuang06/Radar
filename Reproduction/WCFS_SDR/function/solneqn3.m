%% this function is to solve the following two equations with 3 vars.
%%                      imag(conj(s1)*s2)+real(a2*conj(s2))+real(a3*s1)    =0
%%                      real(conj(s1)*s2)+real(b2*conj(s2))+real(b3*s1)    =0
%%    |s1|^2 - |s2|^2+                    real(c2*conj(s2))+real(c3*s1)+c4 =0.

%% In the above, a2, a3, b2, b3, c2, c3 are complex numbers, c4 is real
%% number; they are all inputs.

%% input: a -- 7-dim vector (6 complex numbers, and 1 real number)
%%       a2=a(1);a3=a(2);b2=a(3);b3=a(4);c2=a(5);c3=a(6);c4=a(7);

%% output: s1,s2: a complex-valued solution of the equation system
%%  
%% To run this function, you need: solneqnsys.m


function s=solneqn3(a)

epsi=10^(-8);
a2=a(1);a3=a(2);
b2=a(3);b3=a(4);
c2=a(5);c3=a(6);c4=real(a(7));

u2=abs(a2);xi2=angle(a2);
u3=abs(a3);xi3=angle(a3);
v2=abs(b2);et2=angle(b2);
v3=abs(b3);et3=angle(b3);
w2=abs(c2);zt2=angle(c2);
w3=abs(c3);zt3=angle(c3);

s(1)=0;s(2)=0;
s1=s(1);s2=s(2);
r1=0;tt1=0;
r2=0;tt2=0;

if (length(a)~=7)||(abs(imag(a(7)))>epsi)
    disp('coefficents assumption not satisfied -1 from solneqn3')
    return
end

p3=v2*v3*cos(et2+et3);
p0=-u2*u3*sin(xi2+xi3);
p2=-u2*v3*cos(xi2+et3)-v2*u3*cos(et2+xi3)-v2*v3*sin(et2+et3);
p1=u2*u3*cos(xi2+xi3)+u2*v3*sin(xi2+et3)+v2*u3*sin(et2+xi3);

if abs(p3)>=epsi
    rt=roots([p3,p2,p1,p0]);
    for ii=1:length(rt)
        if abs(imag(rt(ii)))<=epsi
            t=rt(ii);
            break;
        end
    end
    tt20=atan((t*v2*cos(et2)-u2*cos(xi2))/(u2*sin(xi2)-t*v2*sin(et2)));
    tt10=atan((t*v3*cos(et3)-u3*cos(xi3))/(t*v3*sin(et3)-u3*sin(xi3)));
    %%by design, cos(tt20-tt10)~=0
    if abs(cos(tt20-tt10))<epsi
        disp('report -2 from solneqn3') %  cos(theta_2^0-theta_1^0)==0 
        return
    else
        aa=zeros(9,1);
        aa(1)=cos(tt20-tt10);aa(2)=v2*cos(et2-tt20);aa(3)=v3*cos(et3+tt10);aa(4)=1;
        aa(5)=-1;aa(6)=0;aa(7)=w2*cos(zt2-tt20);aa(8)=w3*cos(zt3+tt10);aa(9)=c4;
        ss=solneqnsys(aa);
        %% here ss(3) should be 1, due to cos(tt20-tt10)~=0;
        r1=ss(1);r2=ss(2);
        if r1<0
            tt10=tt10+pi;r1=-r1;
        end
        if r2<0
            tt20=tt20+pi;r2=-r2;
        end
        s1=r1*exp(i*tt10);
        s2=r2*exp(i*tt20);
    end
elseif abs(p0)>=epsi
    rs=0;
    tt20=atan(cos(et2)/(-sin(et2)));
    tt10=atan(cos(et3)/sin(et3));
    %%by design, cos(tt20-tt10)=0, and sin(tt20-tt10)~=0
    if abs(sin(tt20-tt10))<epsi
        disp('report -3 from solneqn3')  % sin(theta_2^0-theta_1^0)==0
        return
    else
        aa=zeros(9,1);
        aa(1)=sin(tt20-tt10);aa(2)=u2*cos(xi2-tt20);aa(3)=u3*cos(xi3+tt10);aa(4)=1;
        aa(5)=-1;aa(6)=0;aa(7)=w2*cos(zt2-tt20);aa(8)=w3*cos(zt3+tt10);aa(9)=c4;
        ss=solneqnsys(aa);
        %% here ss(3) should be 1, due to sin(tt20-tt10)~=0;
        r1=ss(1);r2=ss(2);
        if r1<0
            tt10=tt10+pi;r1=-r1;
        end
        if r2<0
            tt20=tt20+pi;r2=-r2;
        end
        s1=r1*exp(i*tt10);
        s2=r2*exp(i*tt20);
    end
else  %% abs(p3)==0 && abs(p0)==0
    t=0;
    tt20=atan((t*v2*cos(et2)-u2*cos(xi2))/(u2*sin(xi2)-t*v2*sin(et2)));
    tt10=atan((t*v3*cos(et3)-u3*cos(xi3))/(t*v3*sin(et3)-u3*sin(xi3)));
    %%by design, cos(tt20-tt10)~=0
    if abs(cos(tt20-tt10))<epsi
        disp('report -4 from solneqn3') %  cos(theta_2^0-theta_1^0)==0 
        return
    else
        aa=zeros(9,1);
        aa(1)=cos(tt20-tt10);aa(2)=v2*cos(et2-tt20);aa(3)=v3*cos(et3+tt10);aa(4)=1;
        aa(5)=-1;aa(6)=0;aa(7)=w2*cos(zt2-tt20);aa(8)=w3*cos(zt3+tt10);aa(9)=c4;
        ss=solneqnsys(aa);
        %% here ss(3) should be 1, due to cos(tt20-tt10)~=0;
        r1=ss(1);r2=ss(2);
        if r1<0
            tt10=tt10+pi;r1=-r1;
        end
        if r2<0
            tt20=tt20+pi;r2=-r2;
        end
        s1=r1*exp(i*tt10);
        s2=r2*exp(i*tt20);
    end
end

s(1)=s1;s(2)=s2;

% Test:
% ra=randn(7,1);ia=randn(7,1);a=complex(ra,ia);a(7)=real(a(7))
%eqn1=imag(conj(s1)*s2)+real(a2*conj(s2))+real(a3*s1)
%eqn2=real(conj(s1)*s2)+real(b2*conj(s2))+real(b3*s1)
%eqn3=r1^2-r2^2+real(c2*conj(s2))+real(c3*s1)+c4


