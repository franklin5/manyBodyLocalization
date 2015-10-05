function z=index(site,N,L)

if N==1,
    z=site;
else
    loop=1;
    p=0;
    while loop<N,
        if loop==1,
            for i=1:1:site(loop)-1,
                p=p+factorial(L-i)/factorial(N-loop)/factorial(L-i-N+loop);
            end
            loop=loop+1;
        else
            for i=site(loop-1)+1:1:site(loop)-1,
                p=p+factorial(L-i)/factorial(N-loop)/factorial(L-i-N+loop);
            end
            loop=loop+1;
        end
    end
    p=p+site(N)-site(N-1); % that's fine!
    
    z=p;
end