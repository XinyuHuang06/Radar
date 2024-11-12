%% This function is to solve the following system of equations: 
%%  imag(conj(s(1))*s(2))+real(a(1)*conj(s(2)))+real(a(2)*s(1))+a(3)=0.
%%  real(conj(s(1))*s(2))+real(a(4)*conj(s(2)))+real(a(5)*s(1))+a(6)=0.

%% Inputs(coefficients of the equation): a(1),a(2): complex or real number;
%%                                            a(3): real;
%%                                       a(4),a(5): complex or real number;
%%                                            a(6): real;
%% Output(solution): s(1),s(2).

%% To run this function, you need: soln_1eqn_3n.m





function s=soln_2eqns_3n(a)

epsi=10^(-8);
epsii=10^(-8);
a2=a(1);a3=a(2);a4=real(a(3));
b2=a(4);b3=a(5);b4=real(a(6));

s=zeros(2,1);
s1=s(1);s2=s(2);

%% check in input
if (length(a)~=6)||(abs(imag(a(3)))>=epsii)||(abs(imag(a(6)))>=epsii)
    disp('coefficents assumption not satisfied - Stop 1 from soln_twoeqns_3n')
    return
end
%% end of check in input

s0=zeros(2,1);
s0=soln_1eqn_3n([a2,a3,a4]);
s01=s0(1);s02=s0(2);

a2pr=a2+i*s01; a3pr=a3+i*conj(s02);
b2pr=b2+s01; b3pr=b3+conj(s02); b4pr=real(conj(s01)*s02)+real(b3*s01)+real(b2*conj(s02))+b4;

u2=abs(a2pr);
xi2=angle(a2pr);
u3=abs(a3pr);
xi3=angle(a3pr);
v2=abs(b2pr);eta2=angle(b2pr);v3=abs(b3pr);eta3=angle(b3pr);

r1=0;theta1=0;r2=0;theta2=0;
spr1=0;spr2=0;

if abs(b4pr)<=epsii % b4pr==0
    r1=0;theta1=0;r2=0;theta2=0;
end

if b4pr<-epsii % b4pr<0
        if u2<=epsii % u2=0
            theta2=pi/2-xi3;
            theta1=theta2;
            r2=(-(v2*cos(theta2-eta2)+v3*cos(theta1+eta3))+sqrt((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))^2-4*b4pr))/2;
            r1=r2;
        end
        if (u2>epsii)&&(u3<=epsii) % u2>0 and u3=0
            theta2=pi/2+xi2;
            theta1=theta2;
            r2=(-(v2*cos(theta2-eta2)+v3*cos(theta1+eta3))+sqrt((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))^2-4*b4pr))/2;
            r1=r2;
        end
        if (u2>epsii)&&(u3>epsii) % u2>0 and u3>0
                sxi=xi2+xi3;
                if (abs(sxi-2*pi)<=epsii)||(abs(sxi-pi)<=epsii)||(abs(sxi)<=epsii)||(abs(sxi+pi)<=epsii)||(abs(sxi+2*pi)<=epsii)
                    theta2=pi/2-xi3;
                    theta1=theta2;
                    r2=(-(v2*cos(theta2-eta2)+v3*cos(theta1+eta3))+sqrt((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))^2-4*b4pr))/2;
                    r1=r2;
                else
                    theta2=pi/2+(xi2-xi3)/2;
                    theta1=theta2;
                    lambda=u2/u3;
                    r2=(-(v2*cos(theta2-eta2)+lambda*v3*cos(theta1+eta3))+sqrt((v2*cos(theta2-eta2)+lambda*v3*cos(theta1+eta3))^2-4*lambda*b4pr))/(2*lambda);
                    r1=lambda*r2;
                end
            
        end
end

if b4pr>epsii % b4pr>0
        if u2<=epsii % u2=0
            theta1=pi/2-xi3;
            theta2=theta1-pi;
            r2=((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))+sqrt((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))^2+4*b4pr))/2;
            r1=r2;
        end
        if (u2>epsii)&&(u3<=epsii) % u2>0 and u3=0
            u3<=epsii % u3=0
            theta2=pi/2+xi2;
            theta1=theta2+pi;
            r2=((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))+sqrt((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))^2+4*b4pr))/2;
            r1=r2;
        end
        if (u2>epsii)&&(u3>epsii) % u2>0 and u3>0
                sxi=xi2+xi3;
                if (abs(sxi-2*pi)<=epsii)||(abs(sxi-pi)<=epsii)||(abs(sxi)<=epsii)||(abs(sxi+pi)<=epsii)||(abs(sxi+2*pi)<=epsii)
                    theta1=pi/2-xi3;
                    theta2=theta1-pi;
                    r2=((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))+sqrt((v2*cos(theta2-eta2)+v3*cos(theta1+eta3))^2+4*b4pr))/2;
                    r1=r2;
                else
                    theta1=(xi2-xi3)/2;
                    theta2=theta1-pi;
                    lambda=u2/u3;
                    r2=((v2*cos(theta2-eta2)+lambda*v3*cos(theta1+eta3))+sqrt((v2*cos(theta2-eta2)+lambda*v3*cos(theta1+eta3))^2+4*lambda*b4pr))/(2*lambda);
                    r1=lambda*r2;
                end
        end
end

 
spr1=r1*exp(i*theta1);spr2=r2*exp(i*theta2);
 
s1=spr1+s01;s2=spr2+s02;

s(1)=s1;s(2)=s2;

%Test:
%ia=randn(6,1);ib=randn(6,1);a=complex(ia,ib);a(3)=real(a(3));a(6)=real(a(6));
%s=soln_2eqns_3n(a);

%imag(conj(spr1)*spr2)+real(a2pr*conj(spr2))+real(a3pr*spr1)
%real(conj(spr1)*spr2)+real(b2pr*conj(spr2))+real(b3pr*spr1)+b4pr

%eqn1=imag(conj(s(1))*s(2))+real(a(1)*conj(s(2)))+real(a(2)*s(1))+a(3) %% is it equal to zero?     
%eqn2=real(conj(s(1))*s(2))+real(a(4)*conj(s(2)))+real(a(5)*s(1))+a(6) %% is it equal to zero?     
                
        
                
                