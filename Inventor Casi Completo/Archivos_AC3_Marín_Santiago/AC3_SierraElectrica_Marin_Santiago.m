format long
t_ini=0; %s
t_fin=2; %s
n=481;
dt=(t_fin-t_ini)/(n-1);
l1=45; %longitud bancada mm
r2=75;%longitud manivela mm
r3=170;% longitud acoplador mm 
t2=30*(pi/180); %posición angular inicial para la manivela
res_mat=zeros(n,11);%tiempo (col1),theta1(col2),theta2(col3),theta3(col4),r1(col5),dtheta1_dt(col6),dtheta2_dt(col7),dtheta3_dt(col8),dr1_dt(col9),Vb(col10),Vhojax(col11)
res_mat(:,1)=[t_ini:dt:t_fin]'; %columna de tiempo
t_camb=((2*pi)*2)/(0+(70*(2*pi/60)));
dtheta1dt2=((70*(2*pi/60))-0)/(t_camb);
%Ecuaciones de Gobierno Cinemático
for i=1:n
    if res_mat(i,1)<t_camb
        res_mat(i,7)=0+dtheta1dt2*res_mat(i,1); %dtheta2/dt [rad/s]
        res_mat(i,3)=t2+0*res_mat(i,1)+0.5*(dtheta1dt2)*(res_mat(i,1))^2;%theta2 [rad]
    else
        res_mat(i,7)=(70*(2*pi/60)); %dtheta2/dt [rad/s]
        res_mat(i,3)=res_mat(i,7)*res_mat(i,1)+t2-2*pi; %theta2 [rad]
    end
end
%Análisis de Posición
res_mat(:,4)=pi-asin((l1-r2*sin(res_mat(:,3)))/r3);%theta3 [rad]
res_mat(:,2)=pi+acot((r2*cos(res_mat(:,3))+r3*cos(res_mat(:,4)))/l1);%theta1 [rad]
res_mat(:,5)=l1./(sin(res_mat(:,2)));%r1 mm
for i=1:n
    %Análisis de Velocidad
    A_v=[cos(res_mat(i,2)) -res_mat(i,5)*sin(res_mat(i,2)) r3*sin(res_mat(i,4));sin(res_mat(i,2)) res_mat(i,5)*cos(res_mat(i,2)) -r3*cos(res_mat(i,4));sin(res_mat(i,2)) res_mat(i,5)*cos(res_mat(i,2)) 0];
    B_v=[-r2*res_mat(i,7)*sin(res_mat(i,3));r2*res_mat(i,7)*cos(res_mat(i,3));0];
    x_v=linsolve(A_v,B_v); %x=[dr1_dt;dtheta1_dt;dtheta3_dt]
    res_mat(i,9)=x_v(1); %actualizando dr1_dt en la matriz de resultados mm/s
    res_mat(i,6)=x_v(2);%actualizando dtheta1_dt en la matriz de resultados [rad/s]
    res_mat(i,8)=x_v(3);%actualizando dtheta3_dt en la matriz de resultados [rad/s]  
end
%Velocidad de la Hoja 
res_mat(:,10)=(((r2*res_mat(:,7).*cos(res_mat(:,3)+pi/2)+r3*res_mat(:,8).*cos(res_mat(:,4)+pi/2)).^2+(r2*res_mat(:,7).*sin(res_mat(:,3)+pi/2)+r3*res_mat(:,8).*sin(res_mat(:,4)+pi/2)).^2).^0.5)/1000; %Velocidad Abs de la Hoja en m/s
res_mat(:,11)=(r2*res_mat(:,7).*cos(res_mat(:,3)+pi/2)+r3*res_mat(:,8).*cos(res_mat(:,4)+pi/2))/1000; %Velocidad en x de la Hoja m/s

%Conversiones para la Comparación
res_mat(:,2:4)=rad2deg(res_mat(:,2:4)); %conversión a grados de las posiciones theta1,theta2,theta3
res_mat(:,6:8)=rad2deg(res_mat(:,6:8)); %conversión a grados de las velocidades dtheta1_dt,dtheta2_dt,dtheta3_dt

%Importar resultados de Excel conteniendo la soulción numérica obtenida en
%Inventor
res_mat_inventor=readmatrix('Datos_Inventor.xlsx','Sheet','Data','Range','A2:F482');

%gráfica comparativa de la velocidad angular del acoplador
f1=figure;
figure(f1);
plot(res_mat(:,1),res_mat(:,8),'b-o',res_mat(:,1),res_mat_inventor(:,5),'g-->')
title('Comparación de resultados de velocidad angular del acoplador')
xlabel('t (s)')
ylabel('Velocidad Angular del Acoplador (grados/s)')
grid on

%gráfica comparativa de la velocidad de la hoja
f2=figure;
figure(f2);
plot(res_mat(:,1),res_mat(:,11),'b-o',res_mat(:,1),res_mat_inventor(:,6),'r--x')
title('Comparación de resultados de velocidad de la hoja')
xlabel('t (s)')
ylabel('Velocidad de la Hoja (m/s)')
grid on