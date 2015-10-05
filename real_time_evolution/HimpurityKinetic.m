function z=HimpurityKinetic(dim,dim2,L,N2,boundary,basis2)

row=[];
col=[];
ele=[];
count=0;
tic;
for jdim2=1:dim2,
    for jdim1=1:dim,
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
    end
end
z=sparse(row,col,ele,dim*dim2,dim*dim2);
clock2=toc; % END OF CLOCK II
fprintf('Hamiltonian for kinetic part of impurity -> time= %10.2f \n',clock2);
    