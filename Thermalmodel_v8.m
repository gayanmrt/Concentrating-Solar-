
function Quarray_out = Thermalmodel_v8(In,DNI,theta_T)
%This code incldes thermal and optical model of LFR system with
%CPC,compound parabolic concentrator, conductive resistance is not considered, one cover temperature value  
%to save data in each run first create v 7.3 mat file %save saveddata.mat Quarray -v7.3;

h=In(1);  
wm=In(2);
sm=In(3); 
nmax=In(4);
T_r=In(5);
nmax=2*nmax;


%Tr=250 + 273.15; %300C receiver temperature

T_a=20 + 273.15; %20C ambient temperature for australia
ste_co=5.67*10^(-8);

ep_p=0.095;
ep_c=0.9;% glass cover tube same emissivity value for both sides
g=9.81; % graviT_ational acceleration m/s2
D_co=0.125; %glass cover diameter m
D_ci=D_co-0.006; %glass cover tube thickness is taken as 3mm
%k_gc=1.2; %glass cover therml conductivity,W/mK
Dr=0.07; %fixed absorber diameter for this design     
PL=1;%plant lenght for simulation
%AsR=3; %ratio between receiver width/gap, w/d
Ar=pi()*Dr*PL; %receiver area 
%sm1=1.1; % gap between first two mirrors, wm+2*sm for check with odd number of mirrors
sm1=sm;

Quarray=zeros(1,17);
            
%totAp=PL*total_apature(h,wm,sm,sm1,nmax);

%qmfw= (sm1/2 + (nmax/2-1)*(wm+sm) + wm/2)/2; %mid point of one side of the mirror field
%fn=sqrt(qmfw*qmfw + h*h);%calculate a single focal length
fn=10.6; %use one focal lenght

mirrorgen_v2(h,nmax,Dr,sm,wm,theta_T,fn); %generate mirror field

system('"C:\SolTrace\2012.7.9\SolTrace.exe" -s "C:\Users\sirimanm\Desktop\DR\4th year\codes\20211001_LCOE_Yr_33%_Aus\20201112_CPC.lk"'); % run soltrace

%read from file
filename = 'Soltrace_Data.txt'; 
fileID = fopen(filename);
textin=textscan(fileID,'%f');
Soltrace_Data=cell2mat(textin);
fclose(fileID);

Opt_eff=Soltrace_Data(1,1); %optical efficiency
Receiver_Hits=Soltrace_Data(2,1); %rays absorbed by receiver
Sun_area=Soltrace_Data(3,1); %Area of ray generation
Ray_count=Soltrace_Data(4,1); %number of rays generated

Power_per_ray=DNI*Sun_area/Ray_count; %Watts per ray taking DNI in W/m2
Irabs_Ar=Power_per_ray*Receiver_Hits;

%wg=Dr*(1+(2/AsR/tan(37/180*pi())));
%Ag=wg*PL;

%h_cp and h_rp are convective heat transfer coefficient and
%radiative heat transfer coefficient from receiver plate/tube respectively.
%h_co and h_ro are convective heat transfer coefficient and radiative heat
%transfer coefficient from glass cover to air respectively. A_p and A_cg
%are surface areas of plate/tube facing towards cavity and glass cover
%respectively.

%k_im=0.04; %W/mK 
%t_im=0.065; %insulation thickness in m
%UL2=k_im/t_im; %receiver top heat loss
%Where, k_im and t_im are thermal conductivity and thickness of insulation material respectively

%   call airprop funtion for air properties, unit is in K
%   y=[Cp (kJ/kgK); Cv (kJ/kgK); (cp/cv); dynamic viscosity(10-5 kg/m s);
%   Prandtl Number; Kinematic Viscosity (10-5 m2/s); Density (kg/m3)];
%
%{
Tr_0=T_a+200;
Tc_0=T_a;
x0 = [Tr_0,Tc_0];
%while abs(T_diff)>0.1 %tepmerature iteration T_po(j+1)
%options=optimset('Display','iter');  
x=fsolve(@root2d,x0);
      
    function F=root2d(x)
           
%p refers to absorber plate while c refers to glass cover
%T_p=x(1); T_c=x(2)

%%%%%%%%%%%%%%%%
%for outer surface of the glass cover tube; convection
con=airprop((x(2)+T_a)/2); % getting air propertis for glass cover outer surface and ambient air
k_a=0.01*con(5);
Pr_co=con(6);
vi_k=10^(-5)*con(7);%(10-5 m2/s)
beta= 1/((x(2)+T_a)/2); %coeficient of volume expansion 

Ra_co=g*beta*abs((x(2)-T_a))*(D_co^3)*Pr_co/(vi_k^2);% Raleigh number, absolute value of Grashof number  

if (Ra_co>10^12 )
fprintf('Ra_co( %d ) is in out of range, 10^12.\n',Ra_co);
else
Nu_co=(0.6 + 0.387*((Ra_co)^(1/6))/((1+(0.559/Pr_co)^(9/16))^(8/27)))^2; %10^3<Re<2.6*10^5
end

h_co=Nu_co*k_a/D_co ; % heat loss from the glass cover 


%%%%%%%%%%%%%%%%
%Q_loss=Q_co_conv + Q_co_rad;

R_ca= 1/((pi()*D_co*PL*h_co) + (pi()*D_co*PL*ste_co*ep_c*((x(2))^2+T_a^2)*(x(2)+T_a))); %thermal resistance between outer glass surface and ambient

%%%%%%%%%%%%%%%%
%for glass tube conduction
%Q_c_cond=2*pi()*k_gc*PL*(T_ci(j)-T_co(j))/log(D_co/D_ci);

%R_cico=log(D_co/D_ci)/(2*pi()*k_gc*PL); %thermal resistance

%calculate glass cover inner temperature
%T_ci(j)= Q_loss*R_cico + T_co(j);



%%%%%%%%%%%%%%%%
%for absorber tube and glass cover inner surface; radiation
R_pc=((1/ep_p) + (1-ep_c)/ep_c*Dr/D_ci)/(pi()*Dr*PL*ste_co*((x(1))^2+x(2)^2)*(x(1)+x(2))); %thermal resistance


R_T=R_pc + R_ca; % total thermal resistance

UL=1/(pi()*Dr*PL)/R_T;

Trmax= T_a + Irabs_Ar/(UL*Ar);


F(1)=x(2)*(R_ca+R_pc)- R_ca*x(1) - R_pc*T_a ; %T_ci(j) - Q_po_rad*log(D_co/D_ci)/(2*pi()*k_gc*PL);

F(2)=x(1)-(Trmax*T_a)^(0.5);
%{ 
Quarray_it=zeros(1,3);
Quarray_it(1,1)=x(1);
Quarray_it(1,2)=x(2);
Quarray_it(1,3)=UL;

%save iteration.mat Quarray_it -v7.3;

%
matObject = matfile('iteration.mat','Writable',true);
[nrowsB,ncolsB] = size(matObject,'Quarray_it');
matObject.Quarray_it((nrowsB+1),:) = Quarray_it(1,:);
%}


end

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%calculate thermal parameters again based on iterated temperatures

%%%%%%%%%%%%%%%%
%for outer surface of the glass cover tube; convection
con=airprop((x(2)+T_a)/2); % getting air propertis for glass cover outer surface and ambient air
k_a=0.01*con(5);
Pr_co=con(6);
vi_k=10^(-5)*con(7);%(10-5 m2/s)
beta= 1/((x(2)+T_a)/2); %coeficient of volume expansion 

Ra_co=g*beta*abs((x(2)-T_a))*(D_co^3)*Pr_co/(vi_k^2);% Raleigh number, absolute value of Grashof number  

if (Ra_co>10^12 )
fprintf('Ra_co( %d ) is in out of range, 10^12.\n',Ra_co);
else
Nu_co=(0.6 + 0.387*((Ra_co)^(1/6))/((1+(0.559/Pr_co)^(9/16))^(8/27)))^2; %10^3<Re<2.6*10^5
end

h_co=Nu_co*k_a/D_co ; % heat loss from the glass cover 


%%%%%%%%%%%%%%%%
%Q_loss=Q_co_conv + Q_co_rad;

R_ca= 1/((pi()*D_co*PL*h_co) + (pi()*D_co*PL*ste_co*ep_c*((x(2))^2+T_a^2)*(x(2)+T_a))); %thermal resistance between outer glass surface and ambient

%%%%%%%%%%%%%%%%
%for glass tube conduction
%Q_c_cond=2*pi()*k_gc*PL*(T_ci(j)-T_co(j))/log(D_co/D_ci);

%R_cico=log(D_co/D_ci)/(2*pi()*k_gc*PL); %thermal resistance

%calculate glass cover inner temperature
%T_ci(j)= Q_loss*R_cico + T_co(j);


%%%%%%%%%%%%%%%%
%for absorber tube and glass cover inner surface; radiation
R_pc=((1/ep_p) + (1-ep_c)/ep_c*Dr/D_ci)/(pi()*Dr*PL*ste_co*((x(1))^2+x(2)^2)*(x(1)+x(2))); %thermal resistance


R_T=R_pc + R_ca; % total thermal resistance
UL=1/(pi()*Dr*PL)/R_T;

Trmax= T_a + Irabs_Ar/(UL*Ar);
T_r=(Trmax*T_a)^(0.5); 
%
%}

T_r_k=T_r;
[UL,T_c]=TempReceiver(T_r_k,T_a,g,D_co,PL,ste_co,ep_c,ep_p,Dr,D_ci);%getting UL and cover temperature


Qu=Irabs_Ar-UL*Ar*(T_r-T_a);
power_on_mir=Irabs_Ar/Opt_eff;

eff_col=Qu/power_on_mir; % collector efficiency

eff_car=1-T_a/T_r; % carnot efficiency
eff_sys_i=eff_col*eff_car; %system efficiency



%
i=1;%save array 
Quarray(i,1)=h;
Quarray(i,2)=Dr;
Quarray(i,3)=wm;
Quarray(i,4)=sm;
Quarray(i,5)=nmax;
Quarray(i,6)=Opt_eff;
Quarray(i,7)=eff_sys_i; %system efficiency
Quarray(i,8)=eff_car; %carnot efficiency
Quarray(i,9)=eff_col; %collector efficiency
Quarray(i,10)=Qu;
Quarray(i,11)=UL; %overall heat transfer coeficient
Quarray(i,12)=T_r; % receiver(optimum receiver) temperaure
Quarray(i,13)=T_c; %cover temperaure
Quarray(i,14)=power_on_mir; %apature area
Quarray(i,15)=DNI; %DNI
Quarray(i,16)=fn; %focal length
Quarray(i,17)=theta_T; %focal length

Quarray_out=Quarray;


%save saveddata.mat Quarray -v7.3;
%{
matObject = matfile('saveddata.mat','Writable',true);
[nrowsB,ncolsB] = size(matObject,'Quarray');
matObject.Quarray((nrowsB+1),:) = Quarray(1,:);
%}
end



