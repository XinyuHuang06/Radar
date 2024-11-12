%% Matlab Code for rank-one decomposition - by Yongwei Huang and Shuzhong Zhang
%% April 2007

%% Plus the additional functions: fcndcmps13_3n.m
%%                                   fcndcmps12.m
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
%%        D (real symmetric, Hermitian)
%%        X (psd, complex Hermitian symmetric, and could be real symmetric)

%% output: Xs: the desired vector satisfying Xs'*A*Xs=trace(A*X), Xs'*B*Xs=trace(B*X), Xs'*C*Xs=trace(C*X)
%%          r: rank of X;
%%       rdel= trace(AX)
%%       rdelb=trace(BX)
%%       redlc=trace(CX)
%%       redld=trace(DX)


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
    
 %% construct D
    %% method 1:
    %D1=randn(n); 
    %D2=randn(n); 
    %D=(D1+i*D2+(D1+i*D2)')/2;
    
    %% method 2:
    %D1=randn(n);
    %UD1=triu(D1); % get D1's upper triangle matrix including diagonal elements
    %LD1=tril(D1);
    %D=UD1+UD1'-diag(diag(UD1))+i*(LD1-LD1');
    
    %% construct D (real symmetric)
    %D=real(D);
    %or
    %D=(D1+D1')/2;
    
    %% generate integer-entry D:
    %% D=round(10*D);   or   D=fix(10*D);
    
    %% generate a positive definite D:
    % rr=n; 
    % D=zeros(n);
    % for ii=1:rr
    %     rdv=randn(n,1);
    %     idv=randn(n,1);
    %     D=D+(rdv+i*idv)*(rdv+i*idv)';
    % end

    % integer-entry D positive definite:
    % rr=n; 
    % D=zeros(n);
    % for ii=1:rr
    %     rdv=round(2*randn(n,1));
    %     idv=round(2*randn(n,1));  % fix, floor, ceil
    %     D=D+(rdv+i*idv)*(rdv+i*idv)';
    % end

    % construct D (pd, real, rank controlled)
    % rr=n;
    % D=zeros(n);
    % for ii=1:rr
    %     rdv=randn(n,1);
    %     D=D+(rdv)*(rdv)';
    % end

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

function [Xs,r,rdel,rdelb,rdelc,rdeld]=fcndcmps14(A,B,C,D,X)

%epsi=10^(-8);
epsii=10^(-8);

%% begin of checking inputs


[rowA,colA]=size(A);
[rowB,colB]=size(B);
[rowC,colC]=size(C);
[rowD,colD]=size(D);
[rowX,colX]=size(X);
if (rowA==colA)&&(rowB==colB)&&(rowC==colC)&&(rowD==colD)&&(rowX==colX)&&(rowX==rowA)&&(rowA==rowB)&&(rowB==rowC)&&(rowC==rowD)
    n=rowA;
else
    error('A, B, C, D and X must be square matrices and the same size.');
    %disp('A, B, C, D and X must be square matrices and the same size.');
    %return;
end
%n=length(diag(X));
if n<=2
    error('Size must be >=3.');
    %disp('Size must be >=3.');
    %return;
end

% check Hermitian or symmetric of A,B,C,D and X 
for ii=1:n-1
    for jj=ii+1:n
        if (abs(real(A(ii,jj)-A(jj,ii)'))>epsii) || (abs(imag(A(ii,jj)-A(jj,ii)'))>epsii) ||(abs(real(B(ii,jj)-B(jj,ii)'))>epsii) || (abs(imag(B(ii,jj)-B(jj,ii)'))>epsii) || (abs(real(C(ii,jj)-C(jj,ii)'))>epsii) || (abs(imag(C(ii,jj)-C(jj,ii)'))>epsii) ||(abs(real(D(ii,jj)-D(jj,ii)'))>epsii) || (abs(imag(D(ii,jj)-D(jj,ii)'))>epsii) || (abs(real(X(ii,jj)-X(jj,ii)'))>epsii)||(abs(imag(X(ii,jj)-X(jj,ii)'))>epsii) 
            error('A, B, C, D and X must be Hermitian.');
            %disp('A, B, C, D and X must be Hermitian.');
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


%% end of checking input


% rank of X
%r=rank(X);
r=0;
for ii=1:n
    if abs(Di(ii,ii))>=epsii
        r=r+1;
    end
end

%define Xs to save the desired vector
Xs=zeros(n,1);

%% delta=trace(A'*X)=trace(X'*A)=trace(A*X')=trace(A*X)
rdel=real(trace(A*X'));
rdelb=real(trace(B*X'));
rdelc=real(trace(C*X'));
rdeld=real(trace(D*X'));

%disp('Ensure that the input data satisfies that there are real numbers a1,a2, a3, a4 such that a1*A+a2*B+a3*C+a4*D is positive definite.');
%disp('Otherwise, the output is unexpected.');


rde=zeros(4,1);
rde(1)=rdel;rde(2)=rdelb;rde(3)=rdelc;rde(4)=rdeld;

flag1=0;
for ii=1:4
    if abs(rde(ii))>=epsii
        flag1=1;
        ii0=ii;
        break;
    end
end

if flag1==0
    error('The input data problem: at least one of tr(AX),tr(BX),tr(CX),tr(DX) is non-zero.');
end

if ii0==1
    Xsy=fcndcmps13_3n(B-rdelb/rdel*A,C-rdelc/rdel*A,D-rdeld/rdel*A,X);
    t=real(Xsy'*A*Xsy)/rdel;
    if t<=0
        error('unexpected exit -0.');
    end
    Xs=Xsy/sqrt(t);
end

if ii0==2
    Xsy=fcndcmps13_3n(A-rdel/rdelb*B,C-rdelc/rdelb*B,D-rdeld/rdelb*B,X);
    t=real(Xsy'*B*Xsy)/rdelb;
    if t<=0
        error('unexpected exit -0.');
    end
    Xs=Xsy/sqrt(t);
end

if ii0==3
    Xsy=fcndcmps13_3n(A-rdel/rdelc*C,B-rdelb/rdelc*C,D-rdeld/rdelc*C,X);
    t=real(Xsy'*C*Xsy)/rdelc;
    if t<=0
        error('unexpected exit -0.');
    end
    Xs=Xsy/sqrt(t);
end

if ii0==4
    Xsy=fcndcmps13_3n(A-rdel/rdeld*D,B-rdelb/rdeld*D,C-rdelc/rdeld*D,X);
    t=real(Xsy'*D*Xsy)/rdeld;
    if t<=0
        error('unexpected exit -0.');
    end
    Xs=Xsy/sqrt(t);
end

Xs;
rdel;
Xs'*A*Xs;
rdelb;
Xs'*B*Xs;
rdelc;
Xs'*C*Xs;
rdeld;
Xs'*D*Xs;






