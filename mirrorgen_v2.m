
function  mirrorgen_v2(h,nmax,Dr,sm,wm,theta_T)

%this code generates the mirror field and CPC receiver data set for soltrace

%cd('C:\Users\sirimanm\Desktop\DR\3rd year\codes\manual\20200226_paper\s10');    

%selected files
%wrmin=0.1;wrgap=0.1;wrmax=1.0;% receiver width, 0.5 for FRESDEMO (m)
%wmmin=0.1;wmgap=0.1;wmmax=1.0;% mirror width, 0.6 for FRESDEMO (m)
%smmin=0.10;smgap=0.05;smmax=0.10;%mirror spacing, 0.25 for FRESDEMO (m)
%sm1=1.1; % gap between first two mirrors, wm+2*sm for check with odd number of mirrors
%nmax=80; % maximum number of mirrors expected
%hmin=4;hgap=1;hmax=16; %maximum receiver height expected

wr=Dr;
d_cp=0.12; %considering a fixed cavity depth,gap between the absorber plate and glass cover
t_g=0.003; %3mm cover glass
%AsR=3; %ratio between receiver width/gap
theta_T=deg2rad(theta_T);

mirror_array=zeros(nmax,11);%for mirror field data
collecor_array=zeros(5,11);%for collector data
extra_array=zeros(1,11);%for other additional data
extra_array(1,1)=theta_T; %passing transveerse angle 

%for wr= wrmin:wrgap:wrmax % generate receiver width set
  %  for sm= smmin:smgap:smmax % generate mirrror spacing set
   %     for wm= wmmin:wmgap:wmmax % generate mirror width set
    %        for h= hmin:hgap:hmax % generate receiver height set 
                
                %define mirror field
                mirror_array(:,2)=0; %Y
                mirror_array(:,3)=0; %Z
                mirror_array(:,5)=0; %Yaimpoint
                mirror_array(:,6)=h; %Zaimpoint
                mirror_array(:,7)=0; %Zrotation
                mirror_array(:,8)=wm; %Apature 
                mirror_array(:,9)=1; %Apature/plant width
                mirror_array(:,10)=0; %flat surface1
                mirror_array(:,11)=0; %flat surface
                
                %define collector 
                collecor_array(:,2)=0; %Y
                collecor_array(:,5)=0; %Yaimpoint
                collecor_array(:,6)=h; %Zaimpoint
                collecor_array(:,7)=0; %Zrotation
                collecor_array(:,9)=1; %Apature/plant width
                collecor_array(:,10)=0; %flat surface1
                collecor_array(:,11)=0; %flat surface
                
                for n=2:2:nmax %generate mirror number set, code written for an even number of mirrors                
                sm1=sm;    
                
                %define mirror field
                hmfw= sm1/2 + (n/2-1)*(wm+sm) + wm/2; %half mirror field width, distance between the center points of the first and the last mirror
                beta_L=(theta_T-atan(hmfw/h))/2; %negative/west side of mirror field
                beta_R=(theta_T+atan(hmfw/h))/2; %negative/west side of mirror field
                
                Xaimpt_L=hmfw+h*tan(beta_L);
                Xaimpt_R=hmfw-h*tan(beta_R);

                mirror_array((nmax-n+2)/2,1)=-hmfw; %X negative x
                mirror_array((nmax+n)/2,1)=hmfw; %X possitive x
                mirror_array((nmax-n+2)/2,4)=-Xaimpt_L; %Xaimpoint negative x
                mirror_array((nmax+n)/2,4)=Xaimpt_R; %Xaimpoint possitive x
                %mirror_array((nmax-n+2)/2,8)=wm*cos((atan(hmfw/h))/2); %apature negative x
                %mirror_array((nmax+n)/2,8)=wm*cos((atan(hmfw/h))/2); %apature possitive x
               
                n_mirror_array=mirror_array((nmax-n+2)/2:(nmax+n)/2,:); %selecting mirror array
                
                %define collector  
                
                %d_cp=wr/AsR;% d_cp, gap between the absorber plate and glass cover
                theta_1=asin(hmfw/(sqrt(h^2+hmfw^2)));
                theta_2=asin((wm/2)/(sqrt(h^2+hmfw^2)));
                phi= (pi/2) - (theta_1+theta_2); %cavity wall inclination; slope 
                
                collecor_array(1,1)= 0; %X absorber
                collecor_array(2,1)= -(wr/2 + d_cp/2/tan(phi)); %X negative x, left wall
                collecor_array(3,1)= (wr/2 + d_cp/2/tan(phi)); %X possitive x, right wall
                collecor_array(4,1)= 0; %X glass cover top
                collecor_array(5,1)= 0; %X glass cover bottom/mirror side
                
                collecor_array(1,3)= h; %Z absorber
                collecor_array(2,3)= h-d_cp/2; %Z left wall
                collecor_array(3,3)= h-d_cp/2; %Z right wall
                collecor_array(4,3)= h-d_cp; %Z glass cover top
                collecor_array(5,3)= h-d_cp-t_g; %Z glass cover bottom/mirror side
                                
                collecor_array(1,4)= 0; %Xaimpoint absorber
                collecor_array(2,4)= -(wr/2 + d_cp/2/sin(phi)/cos(phi)); %Xaimpoint left wall
                collecor_array(3,4)= (wr/2 + d_cp/2/sin(phi)/cos(phi)); %Xaimpoint right wall
                collecor_array(4,4)= 0; %Xaimpoint glass cover top
                collecor_array(5,4)= 0; %Xaimpoint glass cover bottom/mirror side
                
                collecor_array(1,6)=(h+1); %Zaimpoint; to avoid Z and Zaimpoint being coincided
                
                collecor_array(1,8)= wr; %apature absorber
                collecor_array(2,8)= (d_cp/sin(phi)); %apature left wall
                collecor_array(3,8)= (d_cp/sin(phi)); %apature right wall
                collecor_array(4,8)= (wr + 2*d_cp/tan(phi)); %apature glass cover top
                collecor_array(5,8)= (wr + 2*d_cp/tan(phi)); %apature glass cover bottom/mirror side

                out_array=round([n_mirror_array;collecor_array],3);
                
                %smw=(100*sm); %priting without decimals
                %wmw=(100*wm);
                %wrw=(100*wr);
                
                end
                
                out_array=round([out_array;extra_array],3);
                
dlmwrite('mfield.txt', out_array,'Delimiter','\t','newline','pc');%'S%0.2f_W%0.2f_H%d_M%d.txt'
            %end
        %end
    %end
%end
 
%cd('C:\Users\sirimanm\Desktop\DR\3rd year\codes\manual\20200226_paper');
