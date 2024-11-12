%% This function is to solve the following system of equations: 
%%       imag(conj(s(1))*s(2))+real(a(1)*conj(s(2)))+real(a(2)*s(1))+a(3)=0
%%       real(conj(s(1))*s(2))+real(a(4)*conj(s(2)))+real(a(5)*s(1))+a(6)=0
%% (abs(s(1)))^2-(abs(s(2)))^2+real(a(7)*conj(s(2)))+real(a(8)*s(1))+a(9)=0.

%% Inputs(coefficients of the equation): a(1),a(2): complex or real number;
%%                                            a(3): real;
%%                                       a(4),a(5): complex or real number;
%%                                            a(6): real;
%%                                       a(7),a(8): complex or real number;
%%                                            a(9): real;


%% Output(solution): s(1),s(2).

%% To run this function, you need: soln_1eqn_3n.m
%%                                soln_2eqns_3n.m
%%                                     solneqn3.m
%%                                   solneqnsys.m



function s=soln_3eqns_3n2(a)

epsi=10^(-8);
epsii=10^(-8);

a2=a(1);a3=a(2);a4=real(a(3));
b2=a(4);b3=a(5);b4=real(a(6));
c2=a(7);c3=a(8);c4=real(a(9));

s=zeros(2,1);
s1=s(1);s2=s(2);

%% check in input
if (length(a)~=9)||(abs(imag(a(3)))>=epsii)||(abs(imag(a(6)))>=epsii)||(abs(imag(a(9)))>=epsii)
    disp('coefficents assumption not satisfied - Stop 1 from soln_3eqns_3n2')
    return
end
%% end of check in input

s0=zeros(2,1);
s0=soln_2eqns_3n([a2,a3,a4,b2,b3,b4]);
s01=s0(1);s02=s0(2);

a2pr=a2+i*s01;a3pr=a3+i*conj(s02);
b2pr=b2+s01;b3pr=b3+conj(s02);
c2pr=c2-2*s02;c3pr=c3+2*conj(s01);c4pr=(abs(s01))^2-(abs(s02))^2+real(c2*conj(s02))+real(c3*s01)+c4;

spr=zeros(2,1);
spr=solneqn3([a2pr,a3pr,b2pr,b3pr,c2pr,c3pr,c4pr]);
spr1=spr(1);spr2=spr(2);

s1=spr1+s01;s2=spr2+s02;

s(1)=s1;s(2)=s2;


%Test:
%ia=randn(9,1);ib=randn(9,1);a=complex(ia,ib);a(3)=real(a(3));a(6)=real(a(6));a(9)=real(a(9));
%s=soln_3eqns_3n2(a);

%imag(conj(spr1)*spr2)+real(a2pr*conj(spr2))+real(a3pr*spr1)
%real(conj(spr1)*spr2)+real(b2pr*conj(spr2))+real(b3pr*spr1)+b4pr

%eqn1=imag(conj(s(1))*s(2))+real(a(1)*conj(s(2)))+real(a(2)*s(1))+a(3) %% is it equal to zero?     
%eqn2=real(conj(s(1))*s(2))+real(a(4)*conj(s(2)))+real(a(5)*s(1))+a(6) %% is it equal to zero?     
%eqn3=(abs(s(1)))^2-(abs(s(2)))^2+real(a(7)*conj(s(2)))+real(a(8)*s(1))+a(9) %% is it equal to zero? 


