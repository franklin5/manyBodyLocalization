function z=Hinteraction(U,dim,dim2,N,N2,basis1,basis2)

row=[];
col=[];
ele=[];
count=0;
tic;
for jdim2=1:dim2,
    for jdim1=1:dim,
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
        
    end
end
z=sparse(row,col,ele,dim*dim2,dim*dim2);
clock2=toc; % END OF CLOCK II
fprintf('Hamiltonian for interaction -> time= %10.2f \n',clock2);