function z=HmajorityPotential(randV,dim,dim2,N,basis1,judge)

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
        ele(count)=sumAA;
    end
end
z=sparse(row,col,ele,dim*dim2,dim*dim2);
clock2=toc; % END OF CLOCK II
fprintf('Hamiltonian for potential part of majority -> time= %10.2f \n',clock2);