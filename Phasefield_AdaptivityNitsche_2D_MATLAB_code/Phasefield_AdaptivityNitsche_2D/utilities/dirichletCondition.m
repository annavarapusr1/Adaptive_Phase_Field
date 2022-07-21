function uCCD = dirichletCondition(X,Xref,nodesCCDStd,nodesCCDRef,k,test)

XdirichletStd = X(nodesCCDStd,:);
XdirichletRef = Xref(nodesCCDRef,:);

uCCDStd = dirichletFunction(XdirichletStd,k,test);

if(XdirichletRef)
    uCCDRef = dirichletFunction(XdirichletRef,k,test);
    uCCD = [uCCDStd(:,1); uCCDRef(:,1); uCCDStd(:,2); uCCDRef(:,2)];
else
    uCCD = [uCCDStd(:,1); uCCDStd(:,2)];
end

