%% Matlab Code for rank-one decomposition - by Yongwei Huang and Shuzhong Zhang
%% April 2007

%% Plus the additional functions:    fcndcmps12.m
%%                                   fcndcmps13.m
%%                                     solneqn3.m
%%                                   solneqnsys.m
%%                                  solneqnsys3.m
%%                                 soln_1eqn_3n.m
%%                                soln_2eqns_3n.m
%%                               soln_3eqns_3n2.m
%%                               soln_3eqns_3n3.m

%% input: A (real symmetric, Hermitian)
%%        B (real symmetric, Hermitian)
%%        C (real symmetric, Hermitian)
%%        X (psd, complex Hermitian symmetric, and could be real symmetric)

%% output: Xs: the desired vector satisfying Xs'*A*Xs=trace(A*X), Xs'*B*Xs=trace(B*X), Xs'*C*Xs=trace(C*X)
%%          r: rank of X;
%%       rdel= trace(AX)
%%       rdelb=trace(BX)
%%       redlc=trace(CX)


%% Given problem size: n
%% Methods to randomly generate A:
    %% construct A (Hermitian):

    %% method 1:
    %A1=randn(n); 
    %A2=randn(n); 
    %A=(A1+i*A2+(A1+i*A2)')/2;
    
    %% method 2:
    %A1=randn(n);
    %UA1=triu(A1); %get A1's upper triangle matrix including diagonal elements
    %LA1=tril(A1);
    %A=UA1+UA1'-diag(diag(UA1))+i*(LA1-LA1');
    
    %% construct A (real symmetric)
    %A=real(A);
    %or
    %A=(A1+A1')/2;
    
    %% generate integer-entry A:
    %% A=round(10*A);    or A=fix(10*A);
    
%% construct B
    %% method 1:
    %B1=randn(n); 
    %B2=randn(n); 
    %B=(B1+i*B2+(B1+i*B2)')/2;
    
    %% method 2:
    %B1=randn(n);
    %UB1=triu(B1); % get B1's upper triangle matrix including diagonal elements
    %LB1=tril(B1);
    %B=UB1+UB1'-diag(diag(UB1))+i*(LB1-LB1');
    
    %% construct B (real symmetric)
    %B=real(B);
    %or
    %B=(B1+B1')/2;
    
    %% generate integer-entry B:
    %% B=round(10*B);   or   B=fix(10*B);
    
%% construct C
    %% method 1:
    %C1=randn(n); 
    %C2=randn(n); 
    %C=(C1+i*C2+(C1+i*C2)')/2;
    
    %% method 2:
    %C1=randn(n);
    %UC1=triu(C1); % get C1's upper triangle matrix including diagonal elements
    %LC1=tril(C1);
    %C=UC1+UC1'-diag(diag(UC1))+i*(LC1-LC1');
    
    %% construct C (real symmetric)
    %C=real(C);
    %or
    %C=(C1+C1')/2;
    
    %% generate integer-entry C:
    %% C=round(10*C);   or   C=fix(10*C);

%% construct X: method 1
%X1=randn(n);X2=randn(n);
%X=(X1+i*X2)*(X1+i*X2)';

% generate integer-entry X:
% X=(round(X1+i*X2))*(round(X1+i*X2))' % fix, floor, ceil 

% construct X (psd, real symmetric)
% X=real(X);

%% construct X: method 2
%rr=n-1; % rank of X, this parameter is only for the use of generating X, controlling the rank of X.
%X=zeros(n);
%for ii=1:rr
%    rdv=randn(n,1);
%    idv=randn(n,1);
%    X=X+(rdv+i*idv)*(rdv+i*idv)';
%end

% generate integer-entry X:
%rr=n-1; 
%X=zeros(n);
%for ii=1:rr
%    rdv=round(2*randn(n,1));
%    idv=round(2*randn(n,1));  % fix, floor, ceil
%    X=X+(rdv+i*idv)*(rdv+i*idv)';
%end

% construct X (psd, real, rank controlled)
%rr=n-1;
%X=zeros(n);
%for ii=1:rr
%    rdv=randn(n,1);
%   X=X+(rdv)*(rdv)';
%end


%% function starts...

function [Xs,r,rdel,rdelb,rdelc]=fcndcmps13_3n(A,B,C,X)


%epsi=10^(-8);
epsii=10^(-8);

%% begin of checking inputs


[rowA,colA]=size(A);
[rowB,colB]=size(B);
[rowC,colC]=size(C);
[rowX,colX]=size(X);
if (rowA==colA)&&(rowB==colB)&&(rowC==colC)&&(rowX==colX)&&(rowX==rowA)&&(rowA==rowB)&&(rowB==rowC)
    n=rowA;
else
    error('A, B, C and X must be square matrices and the same size.');
    %disp('A, B, C and X must be square matrices and the same size.');
    %return;
end
%n=length(diag(X));
if n<=2
    error('Size must be >=3.');
    %disp('Size must be >=3.');
    %return;
end

% check Hermitian or symmetric of A,B,C and X 
for ii=1:n-1
    for jj=ii+1:n
        if (abs(real(A(ii,jj)-A(jj,ii)'))>epsii) || (abs(imag(A(ii,jj)-A(jj,ii)'))>epsii) ||(abs(real(B(ii,jj)-B(jj,ii)'))>epsii) || (abs(imag(B(ii,jj)-B(jj,ii)'))>epsii) || (abs(real(C(ii,jj)-C(jj,ii)'))>epsii) || (abs(imag(C(ii,jj)-C(jj,ii)'))>epsii) || (abs(real(X(ii,jj)-X(jj,ii)'))>epsii)||(abs(imag(X(ii,jj)-X(jj,ii)'))>epsii) 
            error('A, B, C and X must be Hermitian.');
            %disp('A, B, C and X must be Hermitian.');
            %return;
        end
    end
end

% check positive semidefinitness of X
[EE,Di]=eig(X);
for ii=1:n
    if Di(ii,ii)<=-epsii
        error('X must be positive semidefinite.');
        %disp('X must be positive semidefinite.');
        %return;
    end
end

% rank of X
%r=rank(X);

r=0;
for ii=1:n
    if abs(Di(ii,ii))>=epsii
        r=r+1;
    end
end


%% end of checking input




%define Xs to save the desired vector
Xs=zeros(n,1);

%% delta=trace(A'*X)=trace(X'*A)=trace(A*X')=trace(A*X)
rdel=real(trace(A*X'));

rdelb=real(trace(B*X'));

rdelc=real(trace(C*X'));



if r==1    
%    [C,D]=eig(X);
    Xs=sqrt(Di(n,n))*EE(:,n);
    return;
end

if r>=3
    Xss=fcndcmps13(A,B,C,X);
    Xs=sqrt(r)*Xss(:,1);
    return;
end

if r==2
    Xss=fcndcmps12(A,B,X);
    Xs1=Xss(:,1);Xs2=Xss(:,2);
    trc=real(Xs1'*C*Xs1);  
    if abs(trc-rdelc/r)<=epsii   % r=2
        Xs=sqrt(r)*Xs1;
        return;
    end
    % randomly generate Xs3 linearly independent of Xs1, Xs2
    Xs3=zeros(n,1);
    rr=rank(Xs1*Xs1'+Xs2*Xs2'); %% should be 2
    trynb=1000;
    for ii=1:trynb
        rvt=randn(n,1)+i*randn(n,1);
        rrn=rank(Xs1*Xs1'/(norm(Xs1))^2+Xs2*Xs2'/(norm(Xs2))^2+rvt*rvt'/(norm(rvt))^2);
        if rrn>rr % rrn=3
            Xs3=rvt*(norm(Xs1)+norm(Xs2)+norm(rvt))/3;
            break;
        end
    end
    if (ii==trynb)&&(norm(Xs3)<=epsii)
        error('no Xs3 found up to 1000 tries!');
    end
    % Xs3 generated!
    a=zeros(14,1);
    a(1)=2*Xs1'*A*Xs2;a(2)=2*Xs2'*A*Xs3;a(3)=2*Xs3'*A*Xs1;a(4)=real(Xs3'*A*Xs3-rdel/r); % r=2
    a(5)=2*Xs1'*B*Xs2;a(6)=2*Xs2'*B*Xs3;a(7)=2*Xs3'*B*Xs1;a(8)=real(Xs3'*B*Xs3-rdelb/r);
    a(9)=real(Xs1'*C*Xs1-rdelc/r); % a(9)*a(10)<0
    a(10)=real(Xs2'*C*Xs2-rdelc/r);
    a(11)=2*Xs1'*C*Xs2;a(12)=2*Xs2'*C*Xs3;a(13)=2*Xs3'*C*Xs1;a(14)=real(Xs3'*C*Xs3-rdelc/r);
    s=soln_3eqns_3n3(a);
    newy=(s(1)*Xs1+s(2)*Xs2+s(3)*Xs3)/(norm(s));  %% s must be a non-zero soln.
    Xs=newy*sqrt(r); % r=2
end


r;
rdel;
Xs'*A*Xs;
rdelb;
Xs'*B*Xs;
rdelc;
Xs'*C*Xs;
