%% This function is to solve the following equation: 
%%  imag(conj(s(1))*s(2))+real(a(1)*conj(s(2)))+real(a(2)*s(1))+a(3)=0.

%% Inputs(coefficients of the equation): a(1),a(2): complex or real number;
%%                                            a(3): real;
%% Output(solution): s(1),s(2).



function s=soln_1eqn_3n(a)

epsi=10^(-8);
a2=a(1);a3=a(2);a4=real(a(3));

s=zeros(2,1);
s1=s(1);s2=s(2);

%% check in input
if (length(a)~=3)||(abs(imag(a(3)))>=epsi)
    disp('coefficents assumption not satisfied - Stop 1 from soln_1eqn_3n')
    return
end
%% end of check in input

u2=abs(a2);xsi2=angle(a2);
u3=abs(a3);xsi3=angle(a3);

if abs(a4)<epsi % or abs(a4)==0
    s1=0;s2=0;
else
    if (u2<epsi)&&(u3<epsi) % or (u2==0)&&(u3==0)
        s1=1;s2=-a4*i;
    end
    if (u2>=epsi) % or u2~=0
        s1=0;s2=-a4/u2*exp(i*xsi2);
    end
    if (u3>=epsi)&&(u2<epsi) % or u3~=0 && u2==0
        s1=-a4/u3*exp(i*(-xsi3));s2=0;
    end
end
s(1)=s1;s(2)=s2;

% test: 
% ia=randn(3,1);ib=randn(3,1);a=complex(ia,ib);a(3)=real(a(3));
% s=soln_singleeqn_3n(a);
%eqn=imag(conj(s(1))*s(2))+real(a(1)*conj(s(2)))+real(a(2)*s(1))+a(3) %% is it equal to zero?     