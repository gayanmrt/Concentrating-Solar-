function [ totAp ] = total_apature( h,wm,sm,sm1,nmax )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%for even number of mirrors

totAp=0;
for n=2:2:nmax %generate mirror number set
    %if (n==1)
        %apature=wm;
    %else
        
        hmfw= sm1/2 + (n/2-1)*(wm+sm) + wm/2;%hmfw=(n-1)/2*(wm+sm); %half mirror field width, distance between the center points of the first and the last mirror
        apatureH=wm*cos((atan(hmfw/h))/2); 
        apature=2*apatureH; %apature form both sides
        
    %end
    
    totAp=totAp + apature;
    
end


end

