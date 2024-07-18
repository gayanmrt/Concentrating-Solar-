%this code calls optical nd thermal models amd calculate objective fuction
function eff_sys_av=objectmain(In)



DNI_array=[700;600;300]; %W/m2, DNI at 0, 30 , 60 transverse angles
Trnv_array=[0;30;60]; %tranverse angles west from zenith, degrees

eff_sys=zeros(3,20);

for i=1:3
    
eff_sys(i,:)=Thermalmodel_v8(In,DNI_array(i),Trnv_array(i));

end

if eff_sys(1,7)<0
eff_sys(1,7)=0;
end

if eff_sys(2,7)<0
eff_sys(2,7)=0;
end

if eff_sys(3,7)<0
eff_sys(3,7)=0;
end

eff_sys_av=-((eff_sys(1,7)*1+eff_sys(2,7)*4+eff_sys(3,7)*4)/9); % toolbox minimises, gx=-fx
end


%}
