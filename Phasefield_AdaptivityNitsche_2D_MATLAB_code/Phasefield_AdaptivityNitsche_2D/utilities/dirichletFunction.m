function u_D = dirichletFunction(X,k,test)

 x = X(:,1); y = X(:,2);
 
u_D = zeros(length(x),2); %[ux,uy]

if (test==1)
    dirichletA = find(abs(y-0.5)<1.e-6);
    u_D(dirichletA,1) = k;    
elseif (test==5)
    u_D(find(y>80),2) = k;    
end


