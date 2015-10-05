function z=HimpurityPotential(randV,dim,dim2,N2,basis2)

row=[];
col=[];
ele=[];
count=0;
tic;% potential energy of impurity
for jdim2=1:dim2,
    for jdim1=1:dim,
        
        % ----- random potential --------
        count=count+1;
        row(count)=(jdim2-1)*dim+jdim1;
        col(count)=(jdim2-1)*dim+jdim1;
        sumAA=0;
        for jpar2=1:N2,
            sumAA=sumAA+randV(basis2(jpar2,jdim2));
        end
        ele(count)=sumAA;
        
    end
end
z=sparse(row,col,ele,dim*dim2,dim*dim2);
clock2=toc; % END OF CLOCK II
fprintf('Hamiltonian for potential part of impurity -> time= %10.2f \n',clock2);