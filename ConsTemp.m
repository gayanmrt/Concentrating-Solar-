function [UL,T_cover] = ConsTemp(T_r_k,T_a,g,D_co,PL,ste_co,ep_c,ep_p,Dr,D_ci)
%UNTITLED2 Summary of this function goes here
%calculate UL for fixed temperature
%   Detailed explanation goes here



T_r=T_r_k; %receier temparature(K) fixed

%T_ci=zeros(1,iter); %glass cover tube inner temparature
iter=30; %number of iterations to find the receiver temperature
T_c=zeros(1,iter); %glass cover temparature; initial guess
T_c(1)=T_a; %initial guess 

%{
Quarray=zeros(1,17);

%i=1;%index
%for h= 1:94:30 % receiver height 
    %for nmax=2:8:498 % mirror number
        %for wr=0.1:0.2:1.2 % receiver width

Ar=pi()*Dr*PL; %receiver area            
            
totAp=PL*total_apature(h,wm,sm,sm1,nmax);

mirrorgen_v2(h,nmax,Dr,sm,wm); %generate mirror field

system('"C:\SolTrace\2012.7.9\SolTrace.exe" -s "C:\Users\user\Desktop\DR_PC backup\3rd year\corona\20201020\20200716_CPC_v1.lk"'); % run soltrace

%Opt_eff=opticaleff_v1(h,nmax,wr);
%read from file
filename = 'Optical_effi.txt'; 
fileID = fopen(filename);
textin=textscan(fileID,'%f');
Opt_eff=cell2mat(textin);
fclose(fileID);

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

%for %k=1 :(iter-1)%tepmerature iteration for T_p

%}
  

for j=1 :(iter-1)%tepmerature iteration

%p refers to absorber plate while c refers to glass cover

%%%%%%%%%%%%%%%%
%for outer surface of the glass cover tube; convection
con=airprop((T_c(j)+T_a)/2); % getting air propertis for glass cover outer surface and ambient air
k_a=0.01*con(5);
Pr_co=con(6);
vi_k=10^(-5)*con(7);%(10-5 m2/s)
beta= 1/((T_c(j)+T_a)/2); %coeficient of volume expansion 

Ra_co=g*beta*abs((T_c(j)-T_a))*(D_co^3)*Pr_co/(vi_k^2);% Raleigh number, absolute value of Grashof number  

if (Ra_co>10^12 )
fprintf('Ra_co( %d ) is in out of range, 10^12.\n',Ra_co);
else
Nu_co=(0.6 + 0.387*((Ra_co)^(1/6))/((1+(0.559/Pr_co)^(9/16))^(8/27)))^2; %10^3<Re<2.6*10^5
end

h_co=Nu_co*k_a/D_co ; % heat loss from the bottom of the glass cover 

%Q_co_conv=pi()*D_co*PL*h_co*(T_co(j)-T_a);


%%%%%%%%%%%%%%%%
%for outer surface of the glass cover tube; radiation

%Q_co_rad=pi()*D_co*PL*ste_co*ep_c*((T_co(j))^4-T_a^4);


%%%%%%%%%%%%%%%%
%Q_loss=Q_co_conv + Q_co_rad;

R_ca= 1/((pi()*D_co*PL*h_co) + (pi()*D_co*PL*ste_co*ep_c*((T_c(j))^2+T_a^2)*(T_c(j)+T_a))); %thermal resistance between outer glass surface and ambient

%%%%%%%%%%%%%%%%
%for glass tube conduction
%Q_c_cond=2*pi()*k_gc*PL*(T_ci(j)-T_co(j))/log(D_co/D_ci);

%R_cico=log(D_co/D_ci)/(2*pi()*k_gc*PL); %thermal resistance

%calculate glass cover inner temperature
%T_ci(j)= Q_loss*R_cico + T_co(j);



%%%%%%%%%%%%%%%%
%for absorber tube and glass cover inner surface; radiation
R_pc=((1/ep_p) + (1-ep_c)/ep_c*Dr/D_ci)/(pi()*Dr*PL*ste_co*((T_r)^2+T_c(j)^2)*(T_r+T_c(j))); %thermal resistance

%Q_po_rad=(T_po(k)-T_ci(j))/R_pci;


%calculate glass cover outer temperature 
T_c(j+1)= (R_ca*T_r + R_pc*T_a)/(R_ca+R_pc); %T_ci(j) - Q_po_rad*log(D_co/D_ci)/(2*pi()*k_gc*PL);

R_T=R_pc + R_ca; % total thermal resistance

UL=1/(pi()*Dr*PL*R_T);
T_cover=T_c(j+1);
end

