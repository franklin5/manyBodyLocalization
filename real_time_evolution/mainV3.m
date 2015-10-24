clear all,
clc

load Hilbert0.mat;

% Basic Parameter Definition

%judge=0;% Clean background!---This background serves as a heat bath! 
judge=1;% Disordered background!---'energy' bath due to interaction

%boundary=1; % @Open Boundary
 boundary=0; % @Periodic Boundary

L=4;
N=2;
N2=1;

W=10;
U=1;


tmax0=10;
Nt0=50;
dt0=log(tmax0)/Nt0;

Ndis=1;

position=1;

% ------ Derived parameters -------------
dim=factorial(L)/factorial(N)/factorial(L-N); % too large L or N will give non-number: NaN or Inf.
dim2=factorial(L)/factorial(N2)/factorial(L-N2);

initial;

for jdis=1:Ndis,
    
    % -- START: variable declare for each realization of disorder ----------
    density1f0=zeros(L,Nt0);
    density2f0=zeros(L,Nt0);
    vector0=zeros(dim*dim2,1); % initial state vector
    % -- END: variable declare for each realization of disorder ----------
    
    tic; % CLOCK I: time for each loop!!
    
    %randV=-W/2+W*rand([1,L]); % generate a random on-site potential within [-W/2,W/2];
    randV=[ 4.99741748906672        -3.370901246089488      -2.173821947071701      4.472010820172727];
    Construction;
    
    % --------- rescaled Hamiltonian -------------
    row2=1:dim*dim2;
    col2=1:dim*dim2;
    ele2=1+0*row2;
    unit=sparse(row2,col2,ele2,dim*dim2,dim*dim2);
    Hpolaron2=(Hpolaron-b*unit)/a;% Take care!! how to operate the non-zero element is crucial!!
    %---------------------------------------------------
    % -------- End of Cosntruction of Hamiltonian ------
    %=============================================================


    %==============================================================
    % -------- time-dependent evolution --------------
    %-------------------------------------------------
    delta=2*a/dim/L;
    t_c=1/delta;
    fprintf('Characteristic Coherance Time: t_c= %5.2f \n', t_c);
    
    %initial state---zero temperature
    label=1:N;
    for jdim=1:dim,
        slater=zeros(N,N);
        for j1=1:N,%row
            for j2=1:N,%col
                slater(j1,j2)=phi(basis1(j1,jdim),label(j2));
            end
        end
        vector0((p2-1)*dim+jdim)=det(slater);%Note that there is  no permutation of N particle in Fock state
        %representation, so there is no need to times 1/sqrt(N!);
    end
    
%     %initial state---infinite temperature   
%     for jdim=1:dim,% product state between spin up and down: spin down is at L/2;
%         p=(p2-1)*dim+jdim;
%         vector0(p)=1/sqrt(dim);% normalization!!!
%     end
    
    %--------------------------------------------------------------
    % --------- Kernal Polynomial Method ---------------
    %------------------------------------------
    
    
    % ------- first stage: [0,tmax0] ------------ 
    tic;
    cutoff0=round(3*a*tmax0);
    %cutoff0=21;Nt0=2;dt0=1.700598690831078;
    listY=0:dt0:(Nt0-1)*dt0;
    listX=1:1:cutoff0;listX=listX';
    timelist=kron(listX*0+1,listY);
    cutofflist=kron(listX,listY*0+1);
    coeff0=(-1i).^(cutofflist-1).*besselj(cutofflist-1,a*timelist);%coeff(cutoff,Nt)
    
    %-----------------------------------------------------
    % --------- START: time evolution --------------
    WFt0=zeros(dim*dim2,Nt0);
    X1=vector0;%Initial Vector
    WFt0=WFt0+X1*coeff0(1,:);
    X2=Hpolaron2*X1;
    WFt0=WFt0+2*X2*coeff0(2,:);
    
    for j=3:cutoff0,
        X3=2*Hpolaron2*X2-X1;
        WFt0=WFt0+2*X3*coeff0(j,:);
        X1=X2;X2=X3;
    end
    %  ---------- END: time evolution ---------------
    %-------------------------------------------------------
    
    %-------------------------------------------------------
    %--- Measurement I start -------
    for jt=1:Nt0,
        vectort=WFt0(:,jt);
        % ## departure ##
        ALLdepart(jt,jdis)=sqrt(sum(rr.^2.*abs(vectort).^2));
        
        % ## entropy ##
        RDM=zeros(dim2,dim2);
        for row2=1:dim2,
            for col2=row2:dim2,
                for jjj=1:dim,
                    RDM(row2,col2)=RDM(row2,col2)+vectort((row2-1)*dim+jjj)*vectort((col2-1)*dim+jjj)';
                end
                if col2~=row2,
                    RDM(col2,row2)=RDM(row2,col2)';
                end
            end
        end
        [v,d]=eig(RDM);
        ALLentropy(jt,jdis)=-trace(d.*log(d));
        %{
        % ## density distribution of impurity fermion
        for j1=1:dim,
            for j2=1:dim2,
                for jn2=1:N2,
                    density2f0(basis2(jn2,j2),jt)=density2f0(basis2(jn2,j2),jt)+abs(vectort((j2-1)*dim+j1))^2;
                end
            end
        end
        % ## density distribution of majority fermions
        for j2=1:dim2,
            for j1=1:dim,
                for jn=1:N,
                    density1f0(basis1(jn,j1),jt)=density1f0(basis1(jn,j1),jt)+abs(vectort((j2-1)*dim+j1))^2;
                end
            end
        end
%}    
    end
    %---- Measurement I end -----------
    %--------------------------------------------
    clock3=toc;
    fprintf('## Time of First Stage Evolution: %10.2f \n',clock3);
    
    
end

ALLdepart

%{
space=1:L;
time0=0:dt0:tmax0-dt0;
time1=tmax0:dt:tmax0+tmax-dt;
time=zeros(Nt0+Nt,1);
for j=1:Nt0+Nt,
    if j<=Nt0,
        time(j)=time0(j);
    else
        time(j)=time1(j-Nt0);
    end
end

Mdensity1f=Mdensity1f/Ndis;
Mdensity2f=Mdensity2f/Ndis;
Mdensity1f0=Mdensity1f0/Ndis;
Mdensity2f0=Mdensity2f0/Ndis;
Mdepart=sum(ALLdepart,2)/Ndis;
Mentropy=sum(ALLentropy,2)/Ndis;
vecX=1:Ndis;
aaa=kron(Mdepart,vecX);
bbb=kron(Mentropy,vecX);
Ddepart=sqrt(sum((ALLdepart-aaa).^2)/Ndis);
Dentropy=sqrt(sum((ALLentropy-bbb).^2)/Ndis);

%particle numbers of the second half of the lattice
HalfCount=zeros(1,L)+1;
HalfCount=(1+sign(space-L/2))/2.*HalfCount;
HalfNum0=HalfCount*Mdensity2f0;
HalfCount=zeros(1,L)+1;
HalfCount=(1+sign(space-L/2))/2.*HalfCount;
HalfNum1=HalfCount*Mdensity2f;
HalfNum=zeros(1,Nt0+Nt);
for j=1:Nt0+Nt,
    if j<=Nt0,
        HalfNum(j)=HalfNum0(j);
    else
        HalfNum(j)=HalfNum1(j-Nt0);
    end
end

save SystemParameters L N N2 W U judge boundary Ndis position Nt tmax Nt0 tmax0 ...
    dim dim2 randV space time0 time1 time;
save Measurment Mdensity1f Mdensity2f Mdensity1f0 Mdensity2f0 ...
    ALLdepart Mdepart Ddepart ALLentropy Mentropy Dentropy HalfNum EimpKin EimpPot EmajKin EmajPot Eint;
%save('/scratch/ly15/WaveFunctiongT01.mat','WFt');
%}









