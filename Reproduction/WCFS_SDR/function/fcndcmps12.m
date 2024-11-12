%% Matlab Code for rank-one decomposition - by Yongwei Huang and Shuzhong Zhang
%% April 2007


%% input: A (Hermitian, or real symmetric), 
%%        B (Hermitian, or real symmetric), 
%%        X (psd, Hermitian, complex, or real matrix)

%% output: Xs1: each column of Xs1 is the decomposed vector;
%%          r: rank of X;
%%       rdel= trace(AX)/r
%%       rdelb=trace(BX)/r

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
    
%% Two methods to randomly generate X:
   %% construct X: method 1
   %X1=randn(n);X2=randn(n);
   %X=(X1+i*X2)*(X1+i*X2)';
   
   %% construct X (psd, real symmetric)
   %X=real(X);
  
   %% construct X: method 2
   %rr=n-1; %% rank of X, this parameter is only for the use of generating X, controlling the rank of X.
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
   
   %% construct X (psd, real symmetric)
   %rr=n-1;
   %X=zeros(n);
   %for ii=1:rr
   %    rdv=randn(n,1);
   %    X=X+(rdv)*(rdv)';
   %end


%% function starts

function [Xs,r,rdel,rdelb]=fcndcmps12(A,B,X)


epsi=10^(-8);

%% begin of checking inputs

% check size of inputs
[rowA,colA]=size(A);
[rowB,colB]=size(B);
[rowX,colX]=size(X);
if (rowA==colA)&&(rowB==colB)&&(rowX==colX)&&(rowX==rowA)&&(rowA==rowB)
    n=rowA;
else
    error('A, B and X must be square matrices and the same size.');
    %disp('A, B and X must be square matrices and the same size.');
    %return;
end
%n=length(diag(X));
if n<=1
    error('Size must be >=2.');
    %disp('Size must be >=2.');
    %return;
end

% check Hermitian or symmetric of A, B and X 
for ii=1:n-1
    for jj=ii+1:n
        if (abs(real(A(ii,jj)-A(jj,ii)'))>10^(-10)) || (abs(imag(A(ii,jj)-A(jj,ii)'))>10^(-10)) ||(abs(real(B(ii,jj)-B(jj,ii)'))>10^(-10)) || (abs(imag(B(ii,jj)-B(jj,ii)'))>10^(-10)) || (abs(real(X(ii,jj)-X(jj,ii)'))>10^(-10)) || (abs(imag(X(ii,jj)-X(jj,ii)'))>10^(-10)) 
            error('A, B and X must be Hermitian or symmetric.');
            %disp('A, B and X must be Hermitian or symmetric.');
            %return;
        end
    end
end

% check positive semidefinitness of X
%% X must be positive semidefinite

eigX=eig(X);
for ii=1:n
    if eigX(ii)<=-10^(-10)
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
    if abs(eigX(ii))>=10^(-8)   % accuracy adjustable according to applications.
        r=r+1;
    end
end



%define Xs to save decomposed vectors for 1st decomposition, Xs1 for 2nd
%decomposition
Xs=zeros(n);
Xs1=zeros(n);
rdel=0;
rdelb=0;


%% X must be complex, otherwise terminate the function
%if norm(imag(X),'fro')<10^(-10)  %% if norm(vect(imag(X)))<epsi
%    disp('X must be complex Hermitian!')
%    return
%end


%inner product of A and X
%vA=vect(A);
%vX=vect(X);
%delta=real((vA)'*vX);  
%% delta=trace(A'*X)=trace(X'*A)=trace(A*X')=trace(A*X)
delta=real(trace(A*X'));
rdel=delta/r;
%inner product of B and X
%vB=vect(B);
%deltb=real((vB)'*vX);
deltb=real(trace(B*X'));
rdelb=deltb/r;


%% if rank of X <=1, done easily and return quickly

if r<=1
    if r==1
       [C D]=eig(X);
       lmda=diag(D);
%       for jj=1:n
%           if abs(lmda(jj))>epsi
%              jj0=jj;
%              break;
%           end
%       end
% an alternative method to construct Xs when r==1
       [mmmmx,jj0]=max(lmda);
       
       Xs(:,r)=(sqrt(lmda(jj0)))*(C(:,jj0));
       Xs;
       disp('rank-1 of X')
       return
    else
        disp('X=0')
        return
    end
end


%% cholesky decom of X=R'*R=L*L', R upper triangle matrix, L=R', L=R'=(v1,v2),L'=R=(v1';v2'), X=v1*v1'+v2*v2'
%% Ltri=(chol(X))';
%% eigenvalue decom of X:=B*D*B',B=(v1,v2),X=l1*v1*v1'+l2*v2*v2'

%% First step decomposition start

Xu=X;
rt=r; %rank test
for ii=1:(n-1)
    %% find v1 and v2 s.t. v1'*A*v1<rdel, v2'*A*v2>rdel
    [C D]=eig(Xu);
    lmda=diag(D);
    v1=zeros(n,1);
    v2=zeros(n,1);
    for jj=1:n
        if abs(lmda(jj))>=epsi
           vv=sqrt(lmda(jj))*C(:,jj);
           deltemp=real(vv'*A*vv); %or real(vA'*vect(vv*vv')) or real(trace(A*(vv*vv')))
           if (deltemp<rdel)&&(abs(deltemp-rdel)>=epsi)&&(norm(v1)==0)
               v1=vv;
               dltp1=deltemp-rdel;
           elseif (deltemp>rdel)&&(abs(deltemp-rdel)>=epsi)&&(norm(v2)==0)
               v2=vv;
               dltp2=deltemp-rdel;
           end
        else
            continue
        end
        if (norm(v1)>0)&&(norm(v2)>0)
            break; 
        end
    end
    
    if  (norm(v1)==0)||(norm(v2)==0) %%could be norm(v1)==0 && norm(v2)==0, but rarely happen. all equalized, 
        disp('f12 Done-1'); %% no v1, v2 found
        %Xs1=Xu; %to be improved, Xs1 should be not covered totally by Xu, at
        %at the moment the rank of Xu is rt.
        DeXu=zeros(n,rt);
        [C D]=eig(Xu);
        lmda=diag(D);
        
        jj=1;
        for ii=1:n
            if abs(lmda(ii))>=epsi
                DeXu(:,jj)=(sqrt(lmda(ii)))*(C(:,ii)); % rank of Xu = rt
                jj=jj+1;
            end
        end
% an alternative method to fill DeXu:
%        for jj=1:rt                 % rt>=2
%            [mmmmx,jj0]=max(lmda);  % lmda >=0
%            DeXu(:,jj)=(sqrt(mmmmx))*(C(:,jj0));
%            lmda(jj0)=0;
%        end
                
 
        Xs1(:,(r-rt+1):r)=DeXu;
        break;    %% break 1, rt>=2, after the break, norm(v1)==0, norm(v2)==0, 
                  %% i.e., [C D]=eig(Xu), Xu already equalized. -- rarely happen
    end
    %%end find v1 and v2
    
    %%update x1,x2 s.t. x1*x1'+x2*x2'=v1*v1'+v2*v2'
    dltaeqn=-1;
    dltaeqn=4*(real(v2'*A*v1))^2-4*dltp2*dltp1; %% >0
    gama1=(-2*real(v2'*A*v1)-sqrt(dltaeqn))/(2*dltp2);
    gama2=(-2*real(v2'*A*v1)+sqrt(dltaeqn))/(2*dltp2);
    if gama1>0
        gama=gama1;
    else
        gama=gama2;
    end
    newx=(v1+gama*v2)/sqrt(1+gama^2);
    Xs1(:,(r-rt+1))=newx;
    Xu=Xu-newx*newx';
    rt=rt-1;   %    rt: rank(Xu); 
    if rt==1
        break; %% break 2, after the break, rt==1, norm(v1)>0, and norm(v2)>0. -- must happen
    end
    clear C D lmda vv
end

if (norm(v1)>0)&&(norm(v2)>0)&&(rt>=1)  %% continue from break 2
    [C D]=eig(Xu);
    lmda=diag(D);
    
    for jj=1:n
        if abs(lmda(jj))>epsi
            jj0=jj;
            break;
        end
    end
% an alternative method to the last column of Xs
%    [mmmmx,jj0]=max(lmda);
%    Xs1(:,r)=(sqrt(mmmmx))*(C(:,jj0));
    
    Xs1(:,r)=(sqrt(lmda(jj0)))*(C(:,jj0));
end
Xs1;
%% first step by Sturm and Zhang paper, done here.


%% now 2nd step, involving A and B (both Hertitian)
%% Note that A and B could be real symmetric, but X must be complex Hermitian


Xu=Xs1;
rt=r; %rt: rank test, initialized to be r, the rank of X, or Xs
for ii=1:(n-1)
    %% find v1 and v2 s.t. v1'*B*v1>rdelb, v2'*B*v2<rdelb, and
    %% v1'*A*v1=v2'*A*v2=rdel
    v1=zeros(n,1);
    v2=zeros(n,1);
    for jj=1:n
        vv=Xu(:,jj);
        if norm(vv)>=epsi
           deltemp=real(vv'*B*vv); %or real(vB'*vect(vv*vv')) or real(trace(B*(vv*vv')))
           if (deltemp>rdelb)&&(abs(deltemp-rdelb)>=epsi)&&(norm(v1)==0)
               v1=vv;
               dltp1=deltemp-rdelb;
               jj1=jj;
           elseif (deltemp<rdelb)&&(abs(deltemp-rdelb)>=epsi)&&(norm(v2)==0)
               v2=vv;
               dltp2=deltemp-rdelb;
               jj2=jj;
           end
        else
            continue
        end
        if (norm(v1)>0)&&(norm(v2)>0)
            break; 
        end
    end
    
    if  (norm(v1)==0)||(norm(v2)==0) %%could be norm(v1)==0 && norm(v2)==0, rarely happen, all equalized
        disp('f12 Done-2'); %% no v1, v2 found, this case happens only when Xu=0, or Xu equalized with B
        % at the moment the rank of Xu is rt. -- rarely happen
        DeXu=zeros(n,rt);
        %[C D]=eig(Xu*Xu');
        %lmda=diag(D);
        jj=1;
        for ii=1:n
            if norm(Xu(:,ii))>=epsi
                DeXu(:,jj)=Xu(:,ii); % rank of Xu = rt
                jj=jj+1;
            end
        end
        Xs(:,(r-rt+1):r)=DeXu;
        break;    %% break 3, after the break, norm(v1)==0, norm(v2)==0, i.e., Xs1 already equalized with A and B.
    end
    %%end find v1 and v2
    
    %%update x1,x2 s.t. x1*x1'+x2*x2'=v1*v1'+v2*v2'
    dltaeqn=-1;
    w1=v1'*A*v2; w2=v1'*B*v2;
    modu1=abs(w1); modu2=abs(w2);
    alpha1=angle(w1); alpha2=angle(w2);
    dltaeqn=4*(modu2*sin(alpha2-alpha1))^2-4*dltp2*dltp1; %% >0
    gama1=(-2*modu2*sin(alpha2-alpha1)-sqrt(dltaeqn))/(2*dltp1);
    gama2=(-2*modu2*sin(alpha2-alpha1)+sqrt(dltaeqn))/(2*dltp1);
    if gama1>0
        gama=gama1;
    else
        gama=gama2;
    end
    w=gama*exp(i*(alpha1+pi/2));
    newx1=(w*v1+v2)/sqrt(1+gama^2);
    newx2=(-v1+conj(w)*v2)/sqrt(1+gama^2);
%    rdel
%    newx1'*A*newx1
%    newx2'*A*newx2
%    rdelb
%    newx1'*B*newx1
    Xs(:,(r-rt+1))=newx1;
    Xu(:,jj1)=zeros(n,1);
    Xu(:,jj2)=newx2;
%    rank(Xu); 
    rt=rt-1;
    if rt==1
        break; %% break 4, after the break, rt==1, Xu has rank 1, norm(v1)>0, and norm(v2)>0. -- must happen
    end
    clear vv
end

if (norm(v1)>0)&&(norm(v2)>0)&&(rt>=1)  %% continue break 4
    for jj=1:n
        if norm(Xu(:,jj))>epsi
            jj0=jj;
            break;
        end
    end
    Xs(:,r)=Xu(:,jj0);
end
Xs;

%% when the input X is real symmetric, meaningful to test the following:
%testX=zeros(n);
%for ii=1:r
%     testX=testX+real(Xs(:,ii))*imag(Xs(:,ii))'-imag(Xs(:,ii))*real(Xs(:,ii))';
%end
%testX
%%this testX should be zero
