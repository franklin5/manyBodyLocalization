clear all,
clc

load Hilbert0.mat;

% Basic Parameter Definition

%judge=0;% Clean background!---This background serves as a heat bath! 
judge=1;% Disordered background!---'energy' bath due to interaction

%boundary=1; % @Open Boundary
 boundary=0; % @Periodic Boundary

L=2;
N=1;
N2=1;

W=10;
U=1;

tmax=2990;%Calculation time is proportional to tmax!!!
Nt=200;
dt=tmax/Nt;

tmax0=10;
Nt0=200;
dt0=tmax0/Nt0;

Ndis=200;

position=1;

% ------ Derived parameters -------------
dim=factorial(L)/factorial(N)/factorial(L-N); % too large L or N will give non-number: NaN or Inf.
dim2=factorial(L)/factorial(N2)/factorial(L-N2);

% ------ Matrix Definition --------------
Mdensity1f0=zeros(L,Nt0);
Mdensity2f0=zeros(L,Nt0);
Mdensity1f=zeros(L,Nt);
Mdensity2f=zeros(L,Nt);

ALLdepart=zeros(Nt0+Nt,Ndis);
ALLentropy=zeros(Nt0+Nt,Ndis);

site2=zeros(1,N2);
site2(1)=position;
p2=index(site2,N2,L);%index of impurity particle, also the position on the lattice

rr=1-p2:1:L-p2;
mm=zeros(dim,1)+1;
rr=kron(rr',mm);%departure distance!

H0=zeros(L,L);
En=zeros(1,L);

%======================================================================
% ------------ Single Particle State (used to construct initial state) --------------------
%----------------------------------------
for j=1:L,
    if j==1,
        H0(j,j+1)=-1;
    elseif j==L,
        H0(j,j-1)=-1;
    else
        H0(j,j+1)=-1;
        H0(j,j-1)=-1;
    end
end
[phi,d]=eig(H0);
for j=1:L,
    En(j)=d(j,j);
end

fprintf('===================================================================\n');
fprintf(' --- Time Dependent Evolution of Disordered Impurity fermion --- \n');
fprintf('===================================================================\n');
fprintf('Parameters of System: \n');
fprintf('---------------------------------------------------------- \n');
fprintf('SystemSize L=%3d \n Number of Particle N=%3d \n Disorder W=%5.2f \n Interaction U=%5.2f \n', L, N, W, U);
fprintf('Total Loops of Disorder: %3d \n', Ndis);
fprintf('Dimension of Hilbert Space: %d \n', dim*dim2)
fprintf('Total time tmax=%5d',tmax0+tmax)
if boundary==0,
    fprintf('#Periodic Boundary Conditions\n');
else
    fprintf('#Open Boundary Conditions\n');
end
fprintf('---------------------------------------------------------- \n');
for jdis=1:Ndis,
    
    % -- START: variable declare for each realization of disorder ----------
    density1f0=zeros(L,Nt0);
    density2f0=zeros(L,Nt0);
    density1f=zeros(L,Nt);
    density2f=zeros(L,Nt);
    vector0=zeros(dim*dim2,1); % initial state vector
    % -- END: variable declare for each realization of disorder ----------
    
    tic; % CLOCK I: time for each loop!!
    
    randV=-W/2+W*rand([1,L]); % generate a random on-site potential within [-W/2,W/2];
    
    %============================================================
    % ----- Construction of many particle Hamiltonian-----
    %------------------------------------------------------------
    row=[];
    col=[];
    ele=[];
    count=0;
    
    tic; % CLOCK 2: time of constructing Hamiltonian!!
    for jdim2=1:dim2,
        for jdim1=1:dim,
            
            % ----- random potential --------
            count=count+1;
            row(count)=(jdim2-1)*dim+jdim1;
            col(count)=(jdim2-1)*dim+jdim1;
            sumAA=0;
            if judge==1,
                for jpar=1:N,
                    sumAA=sumAA+randV(round(basis1(jpar,jdim1)));
                end
            end
            for jpar2=1:N2,
                sumAA=sumAA+randV(basis2(jpar2,jdim2));
            end
            ele(count)=sumAA;
            
            % ------ interaction term ------------
            count=count+1;
            row(count)=(jdim2-1)*dim+jdim1;
            col(count)=(jdim2-1)*dim+jdim1;
            sumAA=0;
            for jpar=1:N,
                for jpar2=1:N2,
                    if basis1(jpar,jdim1)==basis2(jpar2,jdim2),
                        sumAA=sumAA+U;
                    end
                end
            end
            ele(count)=sumAA;
            
            % --------- hopping of spin down -------------
            for jpar2=1:N2,
                
                % @@@@@@@@@@@@@@@@ Stert: Boundary term @@@@@@@@@@@@@@@@@@@@@@@@@
                if boundary==0,
                    % -------- boundary at site L -----------------
                    if basis2(jpar2,jdim2)==L&&basis2(1,jdim2)~=1,
                        if N2==1,
                            site=1;
                            q=index(site,N2,L);
                        else
                            nnn=2:N2;
                            site=basis2(:,jdim2);
                            site(nnn)=site(nnn-1);
                            site(1)=1;
                            q=index(site,N2,L);
                        end
                        count=count+1;
                        row(count)=(q-1)*dim+jdim1;
                        col(count)=(jdim2-1)*dim+jdim1;
                        ele(count)=-1;
                    end
                    
                    %-------- boundary at site 1 --------------------------
                    if basis2(jpar2,jdim2)==1&&basis2(N2,jdim2)~=L,
                        if N2==1,
                            site=L;
                            q=index(site,N2,L);
                        else
                            nnn=2:N2;
                            site=basis2(:,jdim2);
                            site(nnn-1)=site(nnn);
                            site(N2)=L;
                            q=index(site,N2,L);
                        end
                        count=count+1;
                        row(count)=(q-1)*dim+jdim1;
                        col(count)=(jdim2-1)*dim+jdim1;
                        ele(count)=-1;
                    end
                end
                % @@@@@@@@@@@@@@@@ End: Boundary term @@@@@@@@@@@@@@@@@@@@@@@@@
                
                if jpar2<N2&&basis2(jpar2,jdim2)+1<basis2(jpar2+1,jdim2)||jpar2==N2&&basis2(jpar2,jdim2)<L;
                    site=basis2(:,jdim2);
                    site(jpar2)=site(jpar2)+1;
                    q=index(site,N2,L);
                    count=count+1;
                    row(count)=(q-1)*dim+jdim1;
                    col(count)=(jdim2-1)*dim+jdim1;
                    ele(count)=-1;
                end
                if jpar2==1&&basis2(jpar2,jdim2)>1||jpar2>1&&basis2(jpar2,jdim2)-1>basis2(jpar2-1,jdim2),
                    site=basis2(:,jdim2);
                    site(jpar2)=site(jpar2)-1;
                    q=index(site,N2,L);
                    count=count+1;
                    row(count)=(q-1)*dim+jdim1;
                    col(count)=(jdim2-1)*dim+jdim1;
                    ele(count)=-1;
                end
            end
            
            
            % ----------- hopping of spin up ---------------
            for jpar=1:N,
                
                % @@@@@@@@@@@@@@@@ Start: Boundary term @@@@@@@@@@@@@@@@@@@@@@@@@
                if boundary==0,
                    % -------- boundary terms (at site L) -----------------
                    if basis1(jpar,jdim1)==L&&basis1(1,jdim1)~=1,
                        if N==1,
                            site=1;
                            q=index(site,N2,L);
                        else
                            nnn=2:N;
                            site=basis1(:,jdim1);
                            site(nnn)=site(nnn-1);
                            site(1)=1;
                            q=index(site,N,L);
                        end
                        count=count+1;
                        row(count)=(jdim2-1)*dim+q;
                        col(count)=(jdim2-1)*dim+jdim1;
                        ele(count)=-1;
                    end
                    
                    %-------- boundary terms (at site 1) --------------------------
                    if basis1(jpar,jdim1)==1&&basis1(N,jdim1)~=L,
                        if N==1,
                            site=L;
                            q=index(site,N,L);
                        else
                            nnn=2:N;
                            site=basis1(:,jdim1);
                            site(nnn-1)=site(nnn);
                            site(N)=L;
                            q=index(site,N,L);
                        end
                        count=count+1;
                        row(count)=(jdim2-1)*dim+q;
                        col(count)=(jdim2-1)*dim+jdim1;
                        ele(count)=-1;
                    end
                end
                % @@@@@@@@@@@@@@@@ End: Boundary term @@@@@@@@@@@@@@@@@@@@@@@@@
                
                if jpar<N&&basis1(jpar,jdim1)+1<basis1(jpar+1,jdim1)||jpar==N&&basis1(jpar,jdim1)<L;
                    site=basis1(:,jdim1);
                    site(jpar)=site(jpar)+1;
                    q=index(site,N,L);
                    count=count+1;
                    row(count)=(jdim2-1)*dim+q;
                    col(count)=(jdim2-1)*dim+jdim1;
                    ele(count)=-1;
                end
                if jpar==1&&basis1(jpar,jdim1)>1||jpar>1&&basis1(jpar,jdim1)-1>basis1(jpar-1,jdim1),
                    site=basis1(:,jdim1);
                    site(jpar)=site(jpar)-1;
                    q=index(site,N,L);
                    count=count+1;
                    row(count)=(jdim2-1)*dim+q;
                    col(count)=(jdim2-1)*dim+jdim1;
                    ele(count)=-1;
                end
            end
            
        end
    end
    clock2=toc; % END OF CLOCK II
    fprintf('Loop %d: Time of Hamiltonian Cosntruction: %10.2f \n',jdis,clock2);
    
    % ------- maximum and minimum of energy spectrum ----------
    Hpolaron=sparse(row,col,ele,dim*dim2,dim*dim2);
    [v1,d1]=eigs(Hpolaron,3,'sa');
    [v2,d2]=eigs(-Hpolaron,3,'sa');
    Emax=-d2(1,1);
    Emin=d1(1,1);
    epsilon=0.01;
    a=(Emax-Emin)/(2-epsilon);
    b=(Emax+Emin)/2;
    fprintf('---------------------------------------------------------- \n');
    fprintf('Half Wides of Many-Body Spectrum: a=(Emax-Emin)/2=%5.2f \n', a);
    fprintf('Center of Many-Body Spectrum: b=(Emax+Emin)/2=%5.2f \n', b);
    
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
    end
    %---- Measurement I end -----------
    %--------------------------------------------
    clock3=toc;
    fprintf('## Time of First Stage Evolution: %10.2f \n',clock3);
    
    
    % -------- second stage: [tmax0,tmax] ----------
    cutoff=round(1.5*a*tmax);
    fprintf('Cutoff of Chebyshev Expansion: cutoff= %d \n', cutoff);
    tic; %CLOCK 3: time of Coeff
    listY=0:dt:(Nt-1)*dt;
    listX=1:1:cutoff;listX=listX';
    timelist=kron(listX*0+1,listY);
    cutofflist=kron(listX,listY*0+1);
    coeff=(-1i).^(cutofflist-1).*besselj(cutofflist-1,a*timelist);%coeff(cutoff,Nt)
    clock4=toc; %END CLOCK 3
    fprintf('## Time of Second Stage---Coeffecients: %10.2f \n',clock4);
    
    %-----------------------------------------------------
    % --------- START: time evolution --------------
    tic;
    X1=vectort;
    WFt=zeros(dim*dim2,Nt);
    WFt=WFt+X1*coeff(1,:);
    X2=Hpolaron2*X1;
    WFt=WFt+2*X2*coeff(2,:);
    
    for j=3:cutoff,
        X3=2*Hpolaron2*X2-X1;
        WFt=WFt+2*X3*coeff(j,:);
        X1=X2;X2=X3;
        if mod(j,4000)==0
           toc
           disp('j=')
           j
           tic
           save('/scratch/ly15/WaveFunctiongT01.mat','WFt','X1','X2','j');
        end
    end
    save('/scratch/ly15/WaveFunctiongT01.mat','WFt');
    clock5=toc;
    fprintf('## Time of Second Stage Time Evolution: %10.2f \n',clock5);
    %  ---------- END: time evolution ---------------
    %-------------------------------------------------------
    
    %-------------------------------------------------------
    %  ---------- Measurement of density, entropy, departure  ---------------
    tic;
    for jt=1:Nt,
        vectort=WFt(:,jt);
        % departure
        ALLdepart(jt+Nt0,jdis)=sqrt(sum(rr.^2.*abs(vectort).^2));
        % entropy
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
        ALLentropy(jt+Nt0,jdis)=-trace(d.*log(d));
        % impurity fermion density
        for j1=1:dim,
            for j2=1:dim2,
                for jn2=1:N2,
                    density2f(basis2(jn2,j2),jt)=density2f(basis2(jn2,j2),jt)+abs(vectort((j2-1)*dim+j1))^2;
                end
            end
        end
        % majority fermions density
        for j2=1:dim2,
            for j1=1:dim,
                for jn=1:N,
                    density1f(basis1(jn,j1),jt)=density1f(basis1(jn,j1),jt)+abs(vectort((j2-1)*dim+j1))^2;
                end
            end
        end
        
    end
    clock6=toc;
    fprintf('## Time of Second Stage----Measurement: %10.2f \n',clock6);
    
    Mdensity2f0=Mdensity2f0+density2f0;
    Mdensity1f0=Mdensity1f0+density1f0;
    
    Mdensity2f=Mdensity2f+density2f;
    Mdensity1f=Mdensity1f+density1f;
    
    %-------------------------------------------------------------------------------
    % -------------- Start measurement of Energy ------------
    fprintf('=======================================================================\n');
    fprintf('----------- Start Measurement of Energy -----------------------\n');
    
    EimpKin=zeros(1,Nt0+Nt);
    EimpPot=zeros(1,Nt0+Nt);
    EmajKin=zeros(1,Nt0+Nt);
    EmajPot=zeros(1,Nt0+Nt);
    Eint=zeros(1,Nt0+Nt);
    
    % construction of impurity Hamiltonian
    HimpKin=HimpurityKinetic(dim,dim2,L,N2,boundary,basis2);
    HimpPot=HimpurityPotential(randV,dim,dim2,N2,basis2);
    
    HmajKin=HmajorityKinetic(dim,dim2,L,N,boundary,basis1);
    HmajPot=HmajorityPotential(randV,dim,dim2,N,basis1,judge);
    
    Hint=Hinteraction(U,dim,dim2,N,N2,basis1,basis2);
    for jt0=1:Nt0,
        vectort=WFt0(:,jt0);
        EimpKin(jt0)=vectort'*HimpKin*vectort;
        EimpPot(jt0)=vectort'*HimpPot*vectort;
        EmajKin(jt0)=vectort'*HmajKin*vectort;
        EmajPot(jt0)=vectort'*HmajPot*vectort;
        Eint(jt0)=vectort'*Hint*vectort;
    end
    for jt=1:Nt,
        vectort=WFt(:,jt);
        EimpKin(jt+Nt0)=vectort'*HimpKin*vectort;
        EimpPot(jt+Nt0)=vectort'*HimpPot*vectort;
        EmajKin(jt+Nt0)=vectort'*HmajKin*vectort;
        EmajPot(jt+Nt0)=vectort'*HmajPot*vectort;
        Eint(jt+Nt0)=vectort'*Hint*vectort;
    end
    fprintf('----------- Measurement of Energy Finished! -----------\n')
    %========================================================
    %   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %====================================================================
    
    clock1=toc; % END OF CLOCK I
    fprintf('$$ Loop %d Finished!', jdis);
    fprintf('---------------------------------------------------------- \n');

end

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










