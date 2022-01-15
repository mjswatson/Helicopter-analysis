%User Guide
%Input the specifics of the rotor to be tested in the Rotor Specification
%section. This programme only works currently for a linearly twisted, NACA 0012 aerofoil.
%Then input the range of values for solidity and Theta to be tested and the
%step between the values. Finally specify the density of air at the
%required altitude along with the speed of sound. 
%The choice of output is also given. If a value of 1 is inputed an excel 
%file will be created with the results. The name and save location can be 
%specified below the output choice. Any other value will result in a graph 
%being produced by matlab. The axis of the graph produced can also be
%limited in the user section.
%No further user input is required. 
%To run the programme press run or type in hover_code to the command window.
%Any assumptions made are explained in the code next to where they were
%applied.

clear
%Rotor specification 
%give either RPM or Mtip data
R=9.82;         %Radius of Rotor
b=4;            %number of blades
x0=0.15;        %cut out at root
RPM='n/a';      %input either RPM or Mtip. One left blank fill in 'N/A'
Mtip=0.69;      %Mach number at the tip of the blade
n=20;           %number of blade segments
theta1=-8.67;   %twist rate

%Range of collective pitch looked at
thetastart=10;  %collective pitch start
thetaend=30;    %collective pitch end
thetastep=1;    %collecive pitch step
%Range of solidities to be checked 
solstart=0.00; %first solidity value to be run
solend=0.3;    %Difference between solidity value
solstep=0.005; %Final Solidity value to be run
%whole number of curves must be inputted

%conditions specification
p=1.2256;    %local air density
a=340;       %local speed of sound

%Inputs end
%Output choice
output=2;   %output=1 will produce an excel file with the data any other value will produce a matlab graph


if output==1 %Dont change
    cd 'C:\Users\micha\OneDrive\Documents\Uni stuff\Aerospace\Rotary Wing'; %input file location
    filename='hover data.xlsx'; %Inset filename and type
end

%axis limits. Type inf at the end for no limits
xstart=0;
xend=inf;
ystart=0;
yend=inf;

%end user input


if RPM=='n/a'
    RPM=(Mtip*a*60*180)/(360*pi*R);
end


%creating matrices needed
i=1:1:(n+1);
theta=(thetastart:thetastep:thetaend);
sol=(solstart:solstep:solend);

%setting up counts
I1=0; %theta for loop counter
I2=0; % sigma for loop counter

%loops to find a set of values for all solidity values at all theta values
for theta_use=theta 
   %I1 counts number of theta loops run
    I1=I1+1;
    I2=0;
    for sol_use=sol
       %I2 counts the number of solidity loops run
        I2=I2+1;
       c=(sol_use*R*pi)/b;       %finding c for each solidity. c= blade chord
        A=pi*R^2;               %finding the area of the rotor. A= rotor area
       O=RPM/60*360*pi/180;     %finding the omega. O=rotation speed in rad/s
       r=(R*x0)+(((R*(1-x0))/n)*(i-1)); %blade element radii. r= radius position of the element
       I3=0;                    %setting the counter for the blade element for loop
       for ruse=r
           I3=I3+1;
       x=ruse/R;           %x position of the blade elements
       M=x*R*O/a;       %mach number of the elements. M=mach number
       if M>0.725       %finding the blade lift slope
              adeg=(0.677-0.744*M); %adeg is the lift curve slope per degree
       else
           adeg=((0.1/((1-(M^2))^0.5))-0.01*M);
       end
       arad=adeg*(180/pi); %lift curve slope per radians
       thetadeg=theta_use+x*theta1; %the collective pitch in degrees
       thetarad=thetadeg*(pi/180); %thetarad=collective pitch in radians 
       % induced velocity, v1= induced velocity
       v1=((O*arad*b*c)/(16*pi))*(-1+sqrt(1+((32*pi*thetarad*ruse)/(arad*b*c)))); 
       phi=atan(v1/(ruse*O));   %phi=inflow angle 
       alpharad=thetarad-phi;   %angle of attack alpharad=angle of attack in radians
       alphadeg=alpharad*(180/pi);       %alphadeg=angle of attack in degrees
       if M>0.725           %finding the lift coefficient checking for sonic and stall effects
           alphaL=3.4;      %alphaL=stall angle
           Cl=(0.677-0.744*M)*alphadeg; %Sonic region    Cl=lift coefficient 
       else
             alphaL=15-16*M; 
           if alphadeg>alphaL   %in stall
             K1=0.0233+0.342*M^7.15;
             K2=2.05-(0.95*M);
             Cl=adeg*alphadeg-K1*(alphadeg-alphaL)^K2;
           else 
               Cl=((0.1/((1-(M^2))^0.5))-0.01*M)*alphadeg;  %Normal equation
           end
       end
      
       % drag coefficients CDincomp= incompressible drag coefficient 
       CDincomp=(0.0081+(-350*alphadeg+396*alphadeg^2-63.3*alphadeg.^3+3.66*alphadeg^4)*10^-6); 
       if M>0.725   %Checks if sonic
           alphaDdeg=0; %angle of drag in degrees
           alphaDrad=alphaDdeg*(pi/180);    %alphaDrad= angle of drag in radians
           da=alphadeg-alphaDdeg;
           CD=CDincomp+(0.00035*(da)^2.54+21*(M-0.725)^3.2); %Sonic equation  CD=Drag coefficient
       elseif M<0.1
           CD=CDincomp; 
       else
           alphaDdeg=17-23.4*M;
            alphaDrad=alphaDdeg*(pi/180);
            da=alphadeg-alphaDdeg;
           if da>0     %in this region angle of attack should be greater than angle where it converges. 
              da=da;  %If this isnt the case there is no divergence to account for so CD=CDincomp.
          else
             da=0;
          end
           CD=CDincomp+(0.00067*((da)^2.54));   %Subsonic compressible equation
       end
       
       CTdx=(Cl*b*c*x^2)/(2*pi*R);  %CTdx=elemental thrust coefficient
       
       CTdxa(:,I3)=CTdx;    %CTdxa=all the elemental thrust coefficient placed into a matrix
       xa(:,I3)=x;          %xa=matrix of all x values
       alphaLa(:,I3)=alphaL;%alphaLa=matrix of all stall angles for each element
       Ma(:,I3)=M;          %Ma=Matrix of all Mach numbers for each element
       Cla(:,I3)=Cl;        %Cla=Matrix of all lift coefficients for each element
       alphadega(:,I3)=alphadeg; %alphadega=Matrix of all elemtnts anlges of attack
       end
%trapezium rule done manually as trapz produced the wrong value
CTntl=((xa(end)-xa(1))/(2*(numel(xa)-1)))*(CTdxa(1)+2*(sum(CTdxa)-CTdxa(1)-CTdxa(end))+CTdxa(end)); %CTntl= Thrust coefficient without tip losses
       if CTntl>0.006
    B=1-(((2.27*CTntl-0.01)^0.5)/b);    %B=x position past which is ignored to account for tip loss
         else
    B=1-(0.06/b);
       end
    CTdx2=[CTdxa,interp1(xa,CTdxa,B)]; %CTdx2= Matrix of CTdx value including the one interpolated for B
    x2=[xa,B];  %x2=matrix for x positions including B
   %trapezium rule done manually 
    CT=CTntl+((CTdx2(end)+CTdx2(end-1))/2*(x2(end)-x2(end-1))); %CT= thrust coefficient 
    CQ0dx=sol_use/2.*CD.*xa.^3;  %CQ0dx=elemntal torque profile coefficent for each element
    %trapezium rule done manually 
    CQ0=((xa(end)-xa(1))/(2*(numel(xa)-1)))*(CQ0dx(1)+2*(sum(CQ0dx)-CQ0dx(1)-CQ0dx(end))+CQ0dx(end)); %CQ0=profile torque coefficient
    CQidx=(b.*c.*Cl.*sin(phi).*xa.^3)./(2*pi*R);    %CQidx=elemental induced torque coeffcient matrix
    CQidx2=[CQidx,interp1(xa,CQidx,B)]; %CQidx2 elemental induced torque with the interpolated value for B added
    %trapezium rule done manually Cqi=induced torque coefficient 
    Cqi=(((xa(end)-xa(1))/(2*(numel(xa)-1)))*(CQidx2(1)+2*(sum(CQidx2)-CQidx2(1)-CQidx2(end-1)-CQidx2(end))+CQidx2(end-1)))+((CQidx2(end)+CQidx2(end-1))/2*(x2(end)-1)); % give wrong value so check
    lim=sqrt(2*CT);     %lim=limit for induced power integral
    x3=[lim,xa];        %x3=x value with induced power limit included at the start
    dP=((1-sqrt(1-((2*CT)./(xa.^2)))).^2.*(xa.^3).*(1/CT)); %dP=ratio of swirl induced power for each element
    dP=[(lim^3/CT),dP];%creates first value of dP without the chance of imaginary number due to numerical error
    dP(imag(dP)~=0)=0; %removes imaginary numbers from dP
    %trapezium rule done manually Pratio=swirl induced power ratio
    Pratio=((xa(end)-xa(1))/(2*(numel(xa)-1)))*(dP(2)+2*(sum(dP)-dP(2)-dP(end)-dP(1))+dP(end))+((dP(1)+dP(2))/2*(xa(1)-lim));
    deltaCQi=Cqi*Pratio;    %deltaCQi=torque coefficient due to the wake
    DL=CT*p*(R*O)^2;        %DL=disk loading
    CTo=CT/sol_use;          %CTo=CT over the solidity
    DLCTo=DL*CT/sol_use;     %DLCTo=Disk loading and thrust coefficient factor
    %CalcvmeasureP= Calculated Power/ measured power
    CalcvmeasurP=2.39633*10^-10*DLCTo.^5-4.12582*10^-8*DLCTo.^4+3.32009*10^-6*DLCTo.^3-1.81173*10^-4*DLCTo.^2+7.074*10^-3.*DLCTo+0.941289;
    CQ=(deltaCQi+Cqi+CQ0)*CalcvmeasurP; %CQ=total torque coefficient
    CQo=CQ/sol_use;                      %CQo=total torque coefficient over the solidity
    T=(p*A*(O*R)^2*CT)/1000;            %T=Thrust
    Q=(p*A*R*(O*R)^2*CQ)/1000;          %Q=torque
    P=(p*A*(O*R)^3*CQ)/1000;            %P=Power
  
    %saves values for one solidity run into an array
   CTosol(:,I2)=CTo;    %CTosol=Matrix of all CTo value for each solidity run
   CQosol(:,I2)=CQo;    %CQosol=Matrix of all CQo values for each solidity 
   Psol(:,I2)=P;        %Psol= Power value for each solidity
    end
    
    %saves the values for the solidities for each theta into an array
  CToall(I1,:)=CTosol;  %CToall=All CTo solidity matrices saved for all theta values
  CQoall(I1,:)=CQosol;  %CQoall=All CQo solidity matrices saved for all theta values
  Pall(I1,:)=Psol;      %Pall= All power solidity matrices saved for all theta values
end

% CQoall(imag(CQoall)~=0)=0; %removes complex unfeasible designs


if output==1
   thetad=transpose(theta);
    xlswrite(filename,sol,'CToall','B1');
   xlswrite(filename,thetad,'CToall','A2');
   xlswrite(filename,CToall,'CToall','B2');
    xlswrite(filename,sol,'CQoall','B1');
   xlswrite(filename,thetad,'CQoall','A2')
   xlswrite(filename,CQoall,'CQoall','B2');
else
%plotting
%plotting constant solidity curves
plot(CQoall,CToall);
hold on
%plotting constant theta curves
CQoallt=transpose(CQoall);
CToallt=transpose(CToall);
plot(CQoallt,CToallt);
hold off
xlabel('CQ/\sigma')
ylabel('CT/\sigma')
grid on
xlim([xstart xend])
ylim([ystart yend])
end

    