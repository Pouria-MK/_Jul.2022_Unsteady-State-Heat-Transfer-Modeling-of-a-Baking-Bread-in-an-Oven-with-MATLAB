%Pouria Motahhari (99171099)
clc
clear
%% global parameters (SI Unit)
k=0.223;
cp=3600;
rho=720;
L=14;
alpha_=k/(cp*rho);

bread_l=0.5;
bread_w=0.1;
bread_h=0.005;
bread_v=bread_h*bread_w*bread_l;

q_dot=1350/bread_v;
q_w=20000;
h=2.41;
beta_=h/(cp*rho);
T_amb=inline('(-250/7)*abs(x-7)+300'); %linear relation with x
%T_amb=inline('(-250/49)*(x-7)^2+300'); %parabolic relation with x (we can use both)

%Const. Temp. in width thus we use 2d mesh: 
m=25;
n=11;
delta_x=bread_l/m;
delta_y=bread_h/n;

%---radiation constant: ( hr=sigma*epsilon*(Tw^2+Tsur^2)*(Tw+Tsur) )
%assumptions: (sigma=Stefan–Boltzmann constant=5.699*10^-8)
%in hr: Tw=200C or Tw=473K and Tsur=25C or 298K
%the bread is a dark matter (epsilon = 1)
hr=5.699*10^-8*1*(473^2+298^2)*(473+298);
gamma_=hr/(cp*rho);

%%------------------------------------filling matrix in eq.:(C*t_new=A):

u=0.5; %assumed
taw=1; %time step
t_first=25; %temp. of bread in time=0
T_minimum=25;

%we use these constants to make our matrix smaller and better to read:
rx=alpha_*taw/delta_x^2;
ry=alpha_*taw/delta_y^2;

sx=beta_*taw/delta_x;
sy=beta_*taw/delta_y;

%%to calculate the temperatures with radiations, uncomment the following terms:
gx=0;%gamma_*taw/delta_x;  
gy=0;%gamma_*taw/delta_y;

q_dot_mesh= q_dot*taw/(rho*cp); %*in presentation*
q_w_mesh=q_w*taw/(rho*cp*delta_y); %*in presentation*


C=zeros(m*n); %-------------------------constants matrix
A=zeros(m*n,1);%------------------------knowns matrix

while T_minimum<=200 || T_minimum>=201
%     if T_minimum<=200 %method 1
%         u=u-0.01;
%     else
%         u=u+0.01;
%     end
    u=sqrt(T_minimum/200)*u; %method 2
    %both methods work and can be used (depends on the first assumption)
   
    t_new= zeros(m*n,1);%--------unknowns matrix
    
    for counter=1:L/(u*taw)
        
    
        for M=1:m
            for N=1:n
     
                k=(M-1)*n+N;%its T(m,n) in a non square matrix (*in presentation*)
                x=u*taw; %position x difference in one timestep (left side of the bread)
    
                if N==1
                    if M==1 %top left corner (edge of bread) node of mesh
                        A(k)= t_new(k)+ (rx/2)*(t_new(M*n+N)-t_new(k))+ (ry/2)*(t_new(k+1)-t_new(k))+ q_dot_mesh+ (sx/4)*( T_amb(counter*x)-T_amb((counter+1)*x)-t_new(k) )+ (sy/4)*( T_amb(counter*x)-t_new(k)-T_amb((counter+1)*x) )- (gx/4)*t_new(k)- (gy/4)*t_new(k);
                        C(k,k)=1+rx/2+ry/2-sy/4-sx/4+gy/4+gx/4;
                        C(k,k+1)= -ry/2;
                        C(k,M*n+N)= -rx/2;
    
                    elseif M==m %top right corner (edge of bread) node
                        A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)-t_new(k))+ (ry/2)*(t_new(k+1)-t_new(k))+ q_dot_mesh+ (sx/4)*( T_amb(counter*x+bread_l)-T_amb((counter+1)*x+bread_l)-t_new(k) )+ (sy/4)*( T_amb(counter*x+bread_l)-T_amb((counter+1)*x+bread_l)-t_new(k) )- (gx/4)*t_new(k)- (gy/4)*t_new(k);
                        C(k,k)=1+rx/2+ry/2-sy/4-sx/4+gy/4+gx/4;
                        C(k,k+1)= -ry/2;
                        C(k,(M-2)*n+N)= -rx/2;
    
                    else %top side (surface of bread) nodes of mesh
                        A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)+t_new(M*n+N)-2*t_new(k))+ (ry/2)*(t_new(k+1)-t_new(k))+ q_dot_mesh+ (sy/2)*( T_amb(counter*x+M*delta_x)-T_amb((counter+1)*x+M*delta_x)-t_new(k) )- (gy/2)*t_new(k); %*in presentation*
                        C(k,k)=1+rx+ry/2-sy/2+gy/2;
                        C(k,k+1)= -ry/2;
                        C(k,M*n+N)= -rx/2;
                        C(k,(M-2)*n+N)= -rx/2;
                        
                    end
    
                elseif N==n
                    if M==1 %bottom left corner (edge of bread) node
                        A(k)= t_new(k)+ (rx/2)*(t_new(M*n+N)-t_new(k))+ (ry/2)*(t_new(k-1)-t_new(k))+ q_dot_mesh+ q_w_mesh/2+ (sx/4)*( T_amb(counter*x)-T_amb((counter+1)*x)-t_new(k) )- (gx/4)*t_new(k); %t_counter*x is the x (range(x)=0:14) on the oven
                        C(k,k)=1+rx/2+ry/2-sx/4+gx/4;
                        C(k,k-1)= -ry/2;
                        C(k,M*n+N)= -rx/4;
    
                    elseif M==m %bottom right corner (edge of bread) node
                        A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)-t_new(k))+ (ry/2)*(t_new(k-1)-t_new(k))+ q_dot_mesh+ q_w_mesh/2+ (sx/4)*( T_amb(counter*x+bread_l)-T_amb((counter+1)*x+bread_l)-t_new(k) )- (gx/4)*t_new(k); %M=m therefore the node is at the end of the bread
                        C(k,k)=1+rx/2+ry/2-sx/4+gx/4;
                        C(k,k-1)= -ry/2;
                        C(k,(M-2)*n+N)= -rx/2;
    
                    else %bottom side (surface of bread that is on the oven) nodes
                        A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)+t_new(M*n+N)-2*t_new(k))+ (ry/2)*(t_new(k-1)-t_new(k))+ q_dot_mesh+ q_w_mesh;
                        C(k,k)=1+rx+ry/2;
                        C(k,k-1)= -ry/2;
                        C(k,M*n+N)= -rx/2;
                        C(k,(M-2)*n+N)= -rx/2;
    
                    end
    
                elseif M==1 %left side (vertical surface of bread) nodes
                    A(k)= t_new(k)+ (rx/2)*(t_new(M*n+N)-t_new(k))+ (ry/2)*(t_new(k+1)+t_new(k-1)-2*t_new(k))+ q_dot_mesh+ (sx/2)*( T_amb(counter*x)-T_amb((counter+1)*x)-t_new(k) )- (gx/2)*t_new(k);
                    C(k,k)=1+rx/2+ry-sx/2+gx/2;
                    C(k,k-1)= -ry/2;
                    C(k,k+1)= -ry/2;
                    C(k,M*n+N)= -rx/2;
    
                elseif M==m %right side nodes
                    A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)-t_new(k))+ (ry/2)*(t_new(k+1)+t_new(k-1)-2*t_new(k))+ q_dot_mesh+ (sx/2)*( T_amb(counter*x+bread_l)-T_amb((counter+1)*x+bread_l)-t_new(k) )- (gx/2)*t_new(k);
                    C(k,k)=1+rx/2+ry-sx/2+gx/2;
                    C(k,k-1)= -ry/2;
                    C(k,k+1)= -ry/2;
                    C(k,(M-2)*n+N)= -rx/2;
    
                else %center nodes of mesh
                    A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)+t_new(M*n+N)-2*t_new(k))+ (ry/2)*(t_new(k+1)+t_new(k-1)-2*t_new(k))+ q_dot_mesh;
                    C(k,k)=1+rx+ry;
                    C(k,k-1)= -ry/2;
                    C(k,k+1)= -ry/2;
                    C(k,M*n+N)= -rx/2;
                    C(k,(M-2)*n+N)= -rx/2;
    
                end 
            end
        end
    
        %Gauss-Seidel method-----------------------------------------
        si=size(t_new,1);
        normVal=Inf;
        nmax=1000; %number of maximum iterations which can be reached
        tol=0.001; %Tolerence 
        iter=0;

        while normVal>tol && iter<nmax
            t_old=t_new;
            for i=1:si

                guess=0;
                for j=1:i-1
                    guess=guess+C(i,j)*t_new(j);
                end

                for j=i+1:si
                    guess=guess+C(i,j)*t_old(j);
                end

                t_new(i)=(1/C(i,i))*(A(i)-guess);
            end
            iter=iter+1;
            normVal=norm(t_old-t_new);
        end %end of Gauss-Seidel method------------------------------
    
    end
    T_minimum=min(t_new);
    %disp(['at t= ' num2str(counter) ' : T_min= ' num2str(T_minimum)])

end
disp(['velocity(u) = ' num2str(u)])

bread_mesh=zeros(n,m); %temperature distribution on bread mesh (right view)
for j=1:m*n
    bread_mesh(j)=t_new(j);
end


%%third problem: t-x (at y=0.25cm and in the middle of furnace) plotting:

x_axis=0:bread_l/(m-1):bread_l; %x axis of the plot (x) *length issue will be presented*

%these are the same codes from line 74 to 153 except its meant to be used
%on L=7 the middle of the furnace:
t_new=t_first+ zeros(m*n,1);
for counter=1:L/(u*taw*2)
    

    for M=1:m
        for N=1:n
 
            k=(M-1)*n+N;
            x=u*taw; 

            if N==1
                if M==1 %top left corner (edge of bread) node of mesh
                    A(k)= t_new(k)+ (rx/2)*(t_new(M*n+N)-t_new(k))+ (ry/2)*(t_new(k+1)-t_new(k))+ q_dot_mesh+ (sx/4)*( T_amb(counter*x)-T_amb((counter+1)*x)-t_new(k) )+ (sy/4)*( T_amb(counter*x)-t_new(k)-T_amb((counter+1)*x) )- (gx/4)*t_new(k)- (gy/4)*t_new(k);
                    C(k,k)=1+rx/2+ry/2-sy/4-sx/4+gy/4+gx/4;
                    C(k,k+1)= -ry/2;
                    C(k,M*n+N)= -rx/2;

                elseif M==m %top right corner (edge of bread) node
                    A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)-t_new(k))+ (ry/2)*(t_new(k+1)-t_new(k))+ q_dot_mesh+ (sx/4)*( T_amb(counter*x+bread_l)-T_amb((counter+1)*x+bread_l)-t_new(k) )+ (sy/4)*( T_amb(counter*x+bread_l)-T_amb((counter+1)*x+bread_l)-t_new(k) )- (gx/4)*t_new(k)- (gy/4)*t_new(k);
                    C(k,k)=1+rx/2+ry/2-sy/4-sx/4+gy/4+gx/4;
                    C(k,k+1)= -ry/2;
                    C(k,(M-2)*n+N)= -rx/2;

                else %top side (surface of bread) nodes of mesh
                    A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)+t_new(M*n+N)-2*t_new(k))+ (ry/2)*(t_new(k+1)-t_new(k))+ q_dot_mesh+ (sy/2)*( T_amb(counter*x+M*delta_x)-T_amb((counter+1)*x+M*delta_x)-t_new(k) )- (gy/2)*t_new(k); %*in presentation*
                    C(k,k)=1+rx+ry/2-sy/2+gy/2;
                    C(k,k+1)= -ry/2;
                    C(k,M*n+N)= -rx/2;
                    C(k,(M-2)*n+N)= -rx/2;
                    
                end

            elseif N==n
                if M==1 %bottom left corner (edge of bread) node
                    A(k)= t_new(k)+ (rx/2)*(t_new(M*n+N)-t_new(k))+ (ry/2)*(t_new(k-1)-t_new(k))+ q_dot_mesh+ q_w_mesh/2+ (sx/4)*( T_amb(counter*x)-T_amb((counter+1)*x)-t_new(k) )- (gx/4)*t_new(k); %t_counter*x is the x (range(x)=0:14) on the oven
                    C(k,k)=1+rx/2+ry/2-sx/4+gx/4;
                    C(k,k-1)= -ry/2;
                    C(k,M*n+N)= -rx/4;

                elseif M==m %bottom right corner (edge of bread) node
                    A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)-t_new(k))+ (ry/2)*(t_new(k-1)-t_new(k))+ q_dot_mesh+ q_w_mesh/2+ (sx/4)*( T_amb(counter*x+bread_l)-T_amb((counter+1)*x+bread_l)-t_new(k) )- (gx/4)*t_new(k); %M=m therefore the node is at the end of the bread
                    C(k,k)=1+rx/2+ry/2-sx/4+gx/4;
                    C(k,k-1)= -ry/2;
                    C(k,(M-2)*n+N)= -rx/2;

                else %bottom side (surface of bread that is on the oven) nodes
                    A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)+t_new(M*n+N)-2*t_new(k))+ (ry/2)*(t_new(k-1)-t_new(k))+ q_dot_mesh+ q_w_mesh;
                    C(k,k)=1+rx+ry/2;
                    C(k,k-1)= -ry/2;
                    C(k,M*n+N)= -rx/2;
                    C(k,(M-2)*n+N)= -rx/2;

                end

            elseif M==1 %left side (vertical surface of bread) nodes
                A(k)= t_new(k)+ (rx/2)*(t_new(M*n+N)-t_new(k))+ (ry/2)*(t_new(k+1)+t_new(k-1)-2*t_new(k))+ q_dot_mesh+ (sx/2)*( T_amb(counter*x)-T_amb((counter+1)*x)-t_new(k) )- (gx/2)*t_new(k);
                C(k,k)=1+rx/2+ry-sx/2+gx/2;
                C(k,k-1)= -ry/2;
                C(k,k+1)= -ry/2;
                C(k,M*n+N)= -rx/2;

            elseif M==m %right side nodes
                A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)-t_new(k))+ (ry/2)*(t_new(k+1)+t_new(k-1)-2*t_new(k))+ q_dot_mesh+ (sx/2)*( T_amb(counter*x+bread_l)-T_amb((counter+1)*x+bread_l)-t_new(k) )- (gx/2)*t_new(k);
                C(k,k)=1+rx/2+ry-sx/2+gx/2;
                C(k,k-1)= -ry/2;
                C(k,k+1)= -ry/2;
                C(k,(M-2)*n+N)= -rx/2;

            else %center nodes of mesh
                A(k)= t_new(k)+ (rx/2)*(t_new((M-2)*n+N)+t_new(M*n+N)-2*t_new(k))+ (ry/2)*(t_new(k+1)+t_new(k-1)-2*t_new(k))+ q_dot_mesh;
                C(k,k)=1+rx+ry;
                C(k,k-1)= -ry/2;
                C(k,k+1)= -ry/2;
                C(k,M*n+N)= -rx/2;
                C(k,(M-2)*n+N)= -rx/2;

            end 
        end
    end


    %Gauss-Seidel method-----------------------------------------
        si=size(t_new,1);
        normVal=Inf;
        nmax=1000; %number of maximum iterations which can be reached
        tol=0.001; %Tolerence 
        iter=0;

        while normVal>tol && iter<nmax
            t_old=t_new;
            for i=1:si

                guess=0;
                for j=1:i-1
                    guess=guess+C(i,j)*t_new(j);
                end

                for j=i+1:si
                    guess=guess+C(i,j)*t_old(j);
                end

                t_new(i)=(1/C(i,i))*(A(i)-guess);
            end
            iter=iter+1;
            normVal=norm(t_old-t_new);
        end %end of Gauss-Seidel method------------------------------

end


bread_mesh_half_furnace=zeros(n,m); %temperature distribution on bread mesh in x=L/2=7 (right view)
for j=1:m*n
    bread_mesh_half_furnace(j)=t_new(j);
end

y_axis=bread_mesh_half_furnace((n+1)/2,:); %y axis of the plot (Temperature) ((n+1)/2: odd numbers of vertical nodes)


plot(x_axis,y_axis,'b','LineWidth',2)
grid on
xlabel('x (m)')
ylabel('Temperature (°C)')
title('t-x diagram (at y=0.25cm and in the middle of furnace)')
