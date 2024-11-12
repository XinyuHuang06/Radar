%% Matlab Code for rank-one decomposition - by Yongwei Huang and Shuzhong Zhang
%% April 2007

%% Plus the additional functions: solneqnsys.m
%%                                solneqn3.m
%%                                solneqnsys3.m
%%                                fcndcmps12.m

%% input: A (real symmetric, Hermitian)
%%        B (real symmetric, Hermitian)
%%        C (real symmetric, Hermitian)
%%        X (psd, complex Hermitian symmetric, and could be real symmetric)

%% output: Xs: each column of Xs is the decomposed vector;
%%          r: rank of X;
%%       rdel= trace(AX)/r
%%       rdelb=trace(BX)/r
%%       redlc=trace(CX)/r
%%       clas1=Xs(:,r-1)'*C*Xs(:,r-1)
%%       clas2=Xs(:,r)'*C*Xs(:,r)

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
    
    %% construct integer-entry A
    % A=round(10*A) % or A=fix(10*A);
    
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
    
    %% construct integer-entry B
    % B=round(10*B) % or B=fix(10*B);
    
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
    
    %% construct integer-entry C
    % C=round(10*C) % or C=fix(10*C);

%% construct X: method 1
%X1=randn(n);X2=randn(n);
%X=(X1+i*X2)*(X1+i*X2)';
% construct X (psd, real symmetric)
%X=real(X);

%% construct X: method 2
%rr=n-1; % rank of X, this parameter is only for the use of generating X, controlling the rank of X.
%X=zeros(n);
%for ii=1:rr
%    rdv=randn(n,1);
%    idv=randn(n,1);
%    X=X+(rdv+i*idv)*(rdv+i*idv)';
%end

%% construct integer-entry X:
   %rr=n-1; 
   %X=zeros(n);
   %for ii=1:rr
   %    rdv=floor(randn(n,1));  % ceil, round, fix
   %    idv=floor(randn(n,1));
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

function [Xs,r,rdel,rdelb,rdelc,clas1,clas2]=fcndcmps13(A,B,C,X)


epsi=10^(-8);

%% begin of checking inputs

% check size of inputs

%if isreal(X)
%    disp('X should be complex-valued matrix.');
%    return;
%end

%% when inputs A,B,C,X are real, it can works as well (outputting rank-2
%% decomposition of X);

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

% check Hermitian or symmetric of A, B, C and X 
for ii=1:n-1
    for jj=ii+1:n
        if (abs(real(A(ii,jj)-A(jj,ii)'))>10^(-10)) || (abs(imag(A(ii,jj)-A(jj,ii)'))>10^(-10)) ||(abs(real(B(ii,jj)-B(jj,ii)'))>10^(-10)) || (abs(imag(B(ii,jj)-B(jj,ii)'))>10^(-10)) || (abs(real(C(ii,jj)-C(jj,ii)'))>10^(-10)) || (abs(imag(C(ii,jj)-C(jj,ii)'))>10^(-10)) || (abs(real(X(ii,jj)-X(jj,ii)'))>10^(-10)) || (abs(imag(X(ii,jj)-X(jj,ii)'))>10^(-10)) 
            error('A, B, C and X must be Hermitian.');
            %disp('A, B, C and X must be Hermitian.');
            %return;
        end
    end
end

% check positive semidefinitness of X
eigX=eig(X);
for ii=1:n
    if eigX(ii)<=-10^(-10)
        error('X must be positive semidefinite.');
        %disp('X must be positive semidefinite.');
        %return;
    end
end

% rank of X
%r=rank(X);

r=0;
for ii=1:n
    if abs(eigX(ii))>=10^(-10)
        r=r+1;
    end
end

if r<=2
    error('The rank of X should be >= 3');
    %disp('The rank of X should be >= 3');
    %return;
end

%% end of checking input




%define Xs to save decomposed vectors
Xs=zeros(n);
%inner products 
%vA=vect(A);
%vB=vect(B);
%vC=vect(C);
%vX=vect(X);
%delta=real((vA)'*vX);
%% delta=trace(A'*X)=trace(X'*A)=trace(A*X')=trace(A*X)
delta=real(trace(A*X'));
rdel=delta/r;
%delta=real((vB)'*vX);
delta=real(trace(B*X'));
rdelb=delta/r;
%delta=real((vC)'*vX);
delta=real(trace(C*X'));
rdelc=delta/r;
clas1=0;
clas2=0;

%% if X is Hermitian, then please use function fcndcmps12.m

%if norm(imag(vect(X)))>=epsi % or norm(imag(X),'fro')>=epsi, if no vect.m
%    Xs=fcndcmps12(A,B,X);
%    return
%end

%% if B=0, please use function fcndcmps1.m

%if norm(vect(B))<=epsi  % or norm(B,'fro')<=epsi, if no vect.m
%    Xs=fcndcmps1(A,X);
%    return
%end


%% if rank of X <=2, done return solution infomation

%if (r<=2)
%    Xs=fcndcmps12(A,B,X);
%    disp('rank of X <=2, only return solution satisfying A,B, but not for C')
%    return
%end


rt=r; % rank test, rt from r to 2; rt>=3
Xs1=zeros(n);
Xu=X;
for ii=1:n
    Xs1=fcndcmps12(A,B,Xu); %% rank(Xs1)=rank(Xu)=r-rank(Xs)
    %% find v1 and v2 s.t. v1'*C*v1>rdelc, v2'*C*v2<rdelc, and
    %% v1'*A*v1=v2'*A*v2=rdel,v1'*B*v1=v2'*B*v2=rdelb
    v1=zeros(n,1);
    v2=zeros(n,1);
    for jj=1:n
        vv=Xs1(:,jj);
        if norm(vv)>=epsi
           deltemp=real(vv'*C*vv); %or real(vC'*vect(vv*vv'))   or real(trace(C*(vv*vv')))
           if (deltemp>rdelc)&&(abs(deltemp-rdelc)>=epsi)&&(norm(v1)==0)
               v1=vv;
               dltp1=deltemp-rdelc;
               jj1=jj;
           elseif (deltemp<rdelc)&&(abs(deltemp-rdelc)>=epsi)&&(norm(v2)==0)
               v2=vv;
               dltp2=deltemp-rdelc;
               jj2=jj;
           end
        else
            continue
        end
        if (norm(v1)>0)&&(norm(v2)>0)
            break; 
        end
    end
    
    if  (norm(v1)==0)&&(norm(v2)==0) %%could be norm(v1)==0 || norm(v2)==0
        break;    %% break 1, rarely happen; after the break, rt>=3, norm(v1)==0, norm(v2)==0, i.e., Xs already equalized with A,B,C.
    end
    %%end find v1 and v2
    
    %% find v3 from Xs1 s.t. v3~=0, v3~v1, v3~=v2; rank(Xs1)=rt >=3 at the
    %% moment
    v3=zeros(n,1);
    for jj=1:n
        if (jj~=jj1)&&(jj~=jj2)&&(norm(Xs1(:,jj))>epsi)
            v3=Xs1(:,jj);
            break;
        end
    end
    %% v3 found; if v3=0, the only case is rank(Xs1)=2, which will not
    %% happen here. 
    
    %% find a new x
    a=zeros(12,1);
    a(1)=v1'*A*v2;a(2)=v2'*A*v3;a(3)=v3'*A*v1;
    a(4)=v1'*B*v2;a(5)=v2'*B*v3;a(6)=v3'*B*v1;
    a(7)=dltp1;a(8)=dltp2;a(9)=2*v1'*C*v2;a(10)=2*v2'*C*v3;a(11)=2*v3'*C*v1;a(12)=real(v3'*C*v3-rdelc);
    s=solneqnsys3(a);
    newx=(s(1)*v1+s(2)*v2+s(3)*v3)/(norm(s));  %% s must be a non-zero soln.
    
    Xs(:,(r-rt+1))=newx;
    Xu=Xu-newx*newx';
    rt=rt-1;
    if rt==2
        break; %% break 2, must happen; after the break, rt==2, but norm(v1)>0, and norm(v2)>0.
    end
end

%continue break 1: rt>=3, rarely happen
if rt>=3
    disp('f13 Done-1'); %% no v1, v2 found, this case happens only when Xs1=0, or Xs1 equalized with C
        % at the moment the rank of Xs1 is rt>=3; rank(Xs1)=rank(Xu)=rt,
        % rank(Xs)=r-rt
    DeXs1=zeros(n,rt);
    jj=1;
    for ii=1:n
        if norm(Xs1(:,ii))>epsi
            DeXs1(:,jj)=Xs1(:,ii); % rank of Xs1 = rt
            jj=jj+1;
        end
     end
     Xs(:,(r-rt+1):r)=DeXs1;
end

%continue break 2: rt==2, must happen
if rt==2
    DeXs1=zeros(n,rt);
    Xs1=fcndcmps12(A,B,Xu);
    jj=1;
    for ii=1:n
        if norm(Xs1(:,ii))>epsi
            DeXs1(:,jj)=Xs1(:,ii); % rank of Xs1 = rt
            jj=jj+1;
        end
     end
     Xs(:,(r-rt+1):r)=DeXs1;
end

r;
rdel;
Xs(:,r-2)'*A*Xs(:,r-2);
Xs(:,r-1)'*A*Xs(:,r-1);
Xs(:,r)'*A*Xs(:,r);
rdelb;
Xs(:,r-2)'*B*Xs(:,r-2);
Xs(:,r-1)'*B*Xs(:,r-1);
Xs(:,r)'*B*Xs(:,r);
rdelc;
Xs(:,r-2)'*C*Xs(:,r-2);
clas1=real(Xs(:,r-1)'*C*Xs(:,r-1));
clas2=real(Xs(:,r)'*C*Xs(:,r));
