%this code runs the GA optimization alorithm
for i=1:1
ObjFun=@objectmain;
nvars=4;
LB=[2,0.05,0.05,3];
UB=[20,1.0,1.0,50];
%nonlcon = @constraints;
IntCon=4;
options = optimoptions('ga','FunctionTolerance', 0.001);

[x,fval,exitflag,output,population,scores] = ga(ObjFun,nvars,[],[],[],[],LB,UB,[],IntCon,options);

save(['run_',num2str(i)]);
end