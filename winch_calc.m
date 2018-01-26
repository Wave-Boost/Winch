%% Parameters
% Glider
S 	= 17.95	; % glider wing area [m^2]
Speed 	= 28	; % glider towing speed [m/s]
M 	= 550	; % glider laden mass [kg]
Cdg 	= 0.038	; % glider drag coefficient
Clg 	= 1.19	; % glider lift coefficient
% Cable
t	= 0.004	; % cable thickness [m]
Cdc 	= 1	; % cable drag coefficient
beta	= 0	; % angle between cable drag vector and horizontal [degree]
mu 	= 0.0093; % cable linear density [kg/m]
% Earth surface
Rho 	= 1.255	; % air density [kg/m^3]
g 	= 9.81	; % gravitational constant [m/s^2] 

%% Initialise variables
figure		;
i	= 1	;
Delta	= []	;
x	= []	;
y	= []	;
pace_l	= 500	;
l	= [5000:pace_l:13000];
compt	= 0	;
pace	= 0.1	;
X	= []	;
Yg	= []	;
%% Forces acting on the system
W 	= M*g			; % glider weight [N]
D 	= 0.5*Rho*S*Speed*Speed*Cdg	; % glider drag [N]
L 	= 0.5*Rho*S*Speed*Speed*Clg	; % glider lift [N]
d 	= 0.5*Rho*t*Speed*Speed*Cdc	; % cable linear drag  [N/m]
Yg(1)	= 0.585*l(1)		; % initialise Yg [m] 
%% Computation of the shape
for k=1:length(l)	% shape for different lengths
     Wc 	= l(k)*mu*g		; % cable weight [N]
     Th 	= D+d*Yg(i)		; % horizontal tension component [N]
     a 	= Th/(mu*g)		; % shape factor equation
    if (L-W-Wc)>0
        K1 	= asinh((L-W-Wc)/Th)	; % constant 1
        K2 	= -a*cosh(K1)		; % constant 2
        Xg 	= a*(asinh((L-W)/D)-K1)	; % x-value of cable top end
        y(i) 	= a*cosh(Xg)		;	
 
        fun	= @(x) sqrt(1+(sinh(x/a+K1)).^2)	; % length_curve(a,K1)
        long	= quadgk(fun,0,0)		; % call integral solver

        while long<l(i)
            compt   =compt+pace		;
            long    = quadgk(fun,0,compt)	;
        end
        
        u=1;
        for j	= 1:pace:(compt-pace)
            x(u)	= j			; 
            y(u)	= ((a*cosh(x(u)/a+K1)+K2))	;
            u=u+1;
        end
    else
        x = 0	;
        y = 0	; 
    end

    hold on			;
    plot(x,y)			;
    Yg(i+1)=y(end)		;
    i=i+1			;
    xlabel('Distance in m')		;
    ylabel('Height reached in m')	;
 end


for i=1:length(l)
        Delta(i)=Yg(i)/l(i)	;
        X(i)=l(1)+(i-1)*pace_l	;       
end
 
figure			;
plot(X, Delta)		;
xlabel('Cable length')		;
ylabel('Percentage of cable length converted to height')