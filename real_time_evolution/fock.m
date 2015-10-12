clear all,
clc

% many-body Fock states construction
% time: for L=20,N=10, we need 12min to find out all these basis.
% Further: It should be further simplified by build a index matrix to avoid repeated
% calculation!!

L=3;
N=1;
N2=1;

dim=factorial(L)/factorial(N)/factorial(L-N); % too large L or N will give non-number: NaN or Inf.
basis1=zeros(N,dim);
tic;
site=zeros(N+1,1);
for p=1:dim,    
    loop=1;
    if N==1,
        loop=loop+1;
        site(loop)=p;
    else        
        sumAA=0;
        while loop<N,
            len=L-site(loop);
            num=N-loop;            
            for jcar=len-1:-1:num,
                sum2=sumAA;
                sumAA=sumAA+factorial(jcar)/factorial(num)/factorial(jcar-num);
                if sumAA>=p,
                    loop=loop+1;
                    sumAA=sum2;
                    site(loop)=site(loop-1)+len-jcar;
                    break;
                end
            end            
        end
        site(loop+1)=site(loop)+p-sum2;
    end    
    for j=1:N,
        basis1(j,p)=site(j+1);
    end
end
toc

% basis for spin down fermion
dim2=factorial(L)/factorial(N2)/factorial(L-N2);
basis2=zeros(N2,dim2);

site=zeros(N2+1,1);
for p=1:dim2,
    
    loop=1;
    if N2==1,
        loop=loop+1;
        site(loop)=p;
    else
        
        sumAA=0;
        while loop<N2,
            len=L-site(loop);
            num=N2-loop;
            
            for jcar=len-1:-1:num,
                sum2=sumAA;
                sumAA=sumAA+factorial(jcar)/factorial(num)/factorial(jcar-num);
                if sumAA>=p,
                    loop=loop+1;
                    sumAA=sum2;
                    site(loop)=site(loop-1)+len-jcar;
                    break;
                end
            end
            
        end
        site(loop+1)=site(loop)+p-sum2;
    end
    
    for j=1:N2,
        basis2(j,p)=site(j+1);
    end
end

save Hilbert0 L N N2 dim basis1 basis2;

% %-----------------------------
% % find out the index
% 
% site=basis1(:,1);
% if N==1,
%     p=site;
% else
%     loop=1;
%     p=0;
%     while loop<N,
%         if loop==1,
%             for i=1:1:site(loop)-1,
%                 p=p+factorial(L-i)/factorial(N-loop)/factorial(L-i-N+loop);
%             end
%             loop=loop+1;
%         else
%             for i=site(loop-1)+1:1:site(loop)-1,
%                 p=p+factorial(L-i)/factorial(N-loop)/factorial(L-i-N+loop);
%             end
%             loop=loop+1;
%         end
%     end
%     p=p+site(N)-site(N-1); % that's fine!
% end

























        



