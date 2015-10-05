clear all,
clc

load Hilbert0.mat;

%judge=0;% Clean background!---This background serves as a heat bath! 
judge=1;% Disordered background!---'energy' bath due to interaction

boundary=1; % @Open Boundary
% boundary=0; % @Periodic Boundary

L=6;
N=L/2;
N2=1;

W=10; %stength of random potential.
U=1;

E=0.0000000000000000001000002;

% configuration of spin up
vector=zeros(1,N);
for j=1:N,
    vector(j)=2*j;
end
jtar=index(vector,N,L);

Ndis=5000;

% ------ Derived parameters -------------
dim=factorial(L)/factorial(N)/factorial(L-N); % too large L or N will give non-number: NaN or Inf.
dim2=factorial(L)/factorial(N2)/factorial(L-N2);

Veff=zeros(Ndis,1);% Storing the final effective hopping coefficients V_{eff} after decimation
for jdis=1:Ndis,

	   tic;
    
    randV=W*(rand([1,L])-1/2); % generate a random on-site potential within [-W,W];

%     %=====================================================
%     % single particle eigenstate
%     H0=zeros(L,L);
%     En=zeros(1,L);
%     for j=1:L,
%         if j==1,
%             H0(j,j)=randV(j);
%             H0(j,j+1)=-1;
%         elseif j==L,
%             H0(j,j)=randV(j);
%             H0(j,j-1)=-1;
%         else
%             H0(j,j)=randV(j);
%             H0(j,j+1)=-1;
%             H0(j,j-1)=-1;
%         end
%     end
%     [phi,d]=eig(H0);
%     for j=1:L,
%         En(j)=d(j,j);
%     end
    
    %============================================================
    % ----- Construction of many particle Hamiltonian-----
    %------------------------------------------------------------
    row=[];
    col=[];
    ele=[];
    count=0;
    
%    tic; % CLOCK 2: time of constructing Hamiltonian!!
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
%    clock2=toc; % END OF CLOCK II
%     fprintf('Loop %d: Time of Hamiltonian Cosntruction: %10.2f \n',jdis,clock2);
    
    % solve the maximum and minimum of the energy spectrum!!
    Hpolaron=sparse(row,col,ele,dim*dim2,dim*dim2);
    
    
    %============================================================
    % ----- decimation procedure -----
    %------------------------------------------------------------
    %==========================================================
    
    D0=dim*dim2;
    Df=dim2;
    % reorganization of matrix!!
    X0=1:D0;
    Xf=zeros(D0,1);
    for jd=1:D0,
        if jd<=Df,
            Xf(jd)=(jd-1)*dim+jtar;
        else
            point=ceil((jd-Df)/(dim-1));
            point2=(point-1)*dim+jtar+mod(jd-Df-1,dim-1)+1;
            if point2<=D0,
                Xf(jd)=point2;
            else
                Xf(jd)=mod(point2,D0);
            end
        end
    end
    
    H0=zeros(D0,D0);
    for jd1=1:D0,
        for jd2=jd1:D0,
            if jd1==jd2,
                H0(jd1,jd2)=Hpolaron(Xf(jd1),Xf(jd2));
            else
                H0(jd1,jd2)=Hpolaron(Xf(jd1),Xf(jd2));
                H0(jd2,jd1)=H0(jd1,jd2)';
            end
        end
    end
    
    Hreg=H0;
       
    % decimation
    dimB=D0;
    for jloop=1:round(D0-Df),
        
%         jloop
        
        if jloop>1,
            H0=H1;
        end
        dimB=dimB-1;
        H1=zeros(dimB,dimB);
        for jd1=1:dimB,
            for jd2=1:dimB,
                H1(jd1,jd2)=H0(jd1,jd2)+H0(jd1,dimB+1)*H0(dimB+1,jd2)/(E-H0(dimB+1,dimB+1));
            end
        end
    end
    
    Xf=zeros(Df,1);
    for jd=1:Df,
        if jd==1,
            Xf(jd)=jd;
        elseif jd==2,
            Xf(jd)=Df;
        else
            Xf(jd)=jd-1;
        end
    end
    H02=zeros(Df,Df);
    for jd1=1:Df,
        for jd2=jd1:Df,
            if jd2==jd1,
                H02(jd1,jd2)=H1(jd1,jd2);
            else
                H02(jd1,jd2)=H1(Xf(jd1),Xf(jd2))';
                H02(jd2,jd1)=H02(jd1,jd2)';
            end
        end
    end
    
    dimB=Df;
    for jloop=1:round(Df-2),        
        if jloop>1,
            H02=H1;
        end
        dimB=dimB-1;
        H1=zeros(dimB,dimB);
        for jd1=1:dimB,
            for jd2=1:dimB,
                H1(jd1,jd2)=H02(jd1,jd2)+H02(jd1,dimB+1)*H02(dimB+1,jd2)/(E-H02(dimB+1,dimB+1));
            end
        end
    end

    % the final H1 is the renormalized Hamiltonian in the reduced spin down
    % subspace!!!
    
    Veff(jdis)=H1(1,2);

    clock=toc;
    if mod(jdis,100)==0,
        fprintf('%d /%d is finished! Time cost= %10.2f \n', jdis,Ndis,clock);
        save datagDT040 jdis Veff;
    end
    
end

z=log(abs(Veff));
lnVmean=sum(z)/Ndis;
m1=-20;m2=5;dm=(m2-m1)/20;
x=m1:dm:m2;
y=histc(z,x)/Ndis;

% figure;stem(x,y);

x2=m1:(m2-m1)/1000:m2;
z2=round((x2-m1)/dm)+1;
y2=zeros(length(x2),1);
for j=1:length(x2),
    y2(j)=y(z2(j));
end

figure;
plot(x2,y2,'-');

save datagDT04 x y x2 y2 Veff lnVmean U N N2 L Ndis;


% z=Veff;
% m1=-8;m2=4;dm=(m2-m1)/20;
% x=m1:dm:m2;
% y=histc(z,x)/Ndis;
% figure;stem(x,y);

