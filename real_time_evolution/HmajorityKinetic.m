function z=HmajorityKinetic(dim,dim2,L,N,boundary,basis1)

row=[];
col=[];
ele=[];
count=0;
tic; % majority kinetic energy
for jdim2=1:dim2,
    for jdim1=1:dim,
        for jpar=1:N,
            
            % @@@@@@@@@@@@@@@@ Stert: Boundary term @@@@@@@@@@@@@@@@@@@@@@@@@
            if boundary==0,
                % -------- boundary terms (at site L) -----------------
                if basis1(jpar,jdim1)==L&&basis1(1,jdim1)~=1,
                    if N==1,%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        site=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        q=index(site,N,L);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
z=sparse(row,col,ele,dim*dim2,dim*dim2);
clock2=toc; % END OF CLOCK II
fprintf('Hamiltonian for kinetic part of majority ->time = %10.2f \n',clock2);

