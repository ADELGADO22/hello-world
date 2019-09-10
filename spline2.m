clear
clc
close all
n = input ('Enter n positions, n: ')
T=zeros(1,n)
X=zeros(1,n)
Y=zeros(1,n)
Z=zeros(1,n)

for i=1:n
    fprintf('Enter time t (in seconds) and (x,y,z) position (in meters) values on separate lines: \n', i, i, i)
    T(i)= input ('')
    X(i)= input ('')
    Y(i)= input ('')
    Z(i)= input ('')
end

fig1=figure()
set(fig1,'color','white')
set(gca,'fontsize',18)
p1=plot(T,X,'b*')
hold on
p2=plot(T,Y,'r*')
hold on
p3=plot(T,Z,'g*')
hold on
grid on
hold on


Nx= length(X)
Mx=Nx-1
Coefx=4*Mx

%%create X vector
Bx=zeros(2*(Nx-2)+2+(Nx-2)+(Nx-2)+2,1)

Nrowx=2*(Nx-2)+2
hx0=zeros(Nrowx,Coefx)


%%first part of the matrix
for i=1:(Nrowx/2)
    colx=((i-1)*4+1)
    rowx=(i-1)*2+1
    for j=colx:(colx+3)
        hx0(rowx,j)=T(i).^(3-j+colx)
        hx0(rowx+1,j)=T(i+1).^(3-j+colx)
    end
    Bx(rowx,1)=X(i)
    Bx(rowx+1,1)=X(i+1)
end

%%second part of the matrix using derivatives

hx1=zeros(Nx-2,Coefx)
for i=2:(Nx-1)
    colx=1+(i-2)*4
    hx1(i-1,colx)=3*T(i)^2
    hx1(i-1,colx+1)=2*T(i)
    hx1(i-1,colx+2)=1
    hx1(i-1,colx+3)=0
    hx1(i-1,colx+4)=-3*T(i)^2
    hx1(i-1,colx+5)=-2*T(i)
    hx1(i-1,colx+6)=-1
    hx1(i-1,colx+7)=0
end

%%third part of the matrix boundary conditions
hx2=0*hx1
for i=2:(Nx-1)
    colx=1+(i-2)*4
    hx2(i-1,colx)=6*T(i)
    hx2(i-1,colx+1)=2
    hx2(i-1,colx+2)=0
    hx2(i-1,colx+3)=0
    hx2(i-1,colx+4)=-6*T(i)
    hx2(i-1,colx+5)=-2
    hx2(i-1,colx+6)=0
    hx2(i-1,colx+7)=0
end

hxepts=zeros(2,Coefx)
hxepts(1,1)=6*T(1)
hxepts(1,2)=2*T(1)
hxepts(1,3)=1
hxepts(2, end)=1
hxepts(2, end-1)=2*T(end)
hxepts(2, end-2)=6*T(end)

hx=[hx0;hx1;hx2;hxepts]

%%solve
ABCsx=inv(hx)*Bx

for i=1:Mx
   rowx=1+(i-1)*4
   ax=ABCsx(rowx)
   bx=ABCsx(rowx+1)
   cx=ABCsx(rowx+2)
   dx=ABCsx(rowx+3)
   t=linspace(T(i),T(i+1),100)
   x=ax*t.^3+bx*t.^2+cx*t+dx
   px=plot(t,x,'b')
end


%%%%Y
Ny= length(Y)
My=Ny-1
Coefy=4*My

%%create Y vector
By=zeros(2*(Ny-2)+2+(Ny-2)+(Ny-2)+2,1)

Nrowy=2*(Ny-2)+2
hy0=zeros(Nrowy,Coefy)


%%first part of the matrix
for i=1:(Nrowy/2)
    coly=((i-1)*4+1)
    rowy=(i-1)*2+1
    for j=coly:(coly+3)
        hy0(rowy,j)=T(i).^(3-j+coly)
        hy0(rowy+1,j)=T(i+1).^(3-j+coly)
    end
    By(rowy,1)=Y(i)
    By(rowy+1,1)=Y(i+1)
end

%%second part of the matrix using derivatives

hy1=zeros(Ny-2,Coefy)
for i=2:(Ny-1)
    coly=1+(i-2)*4
    hy1(i-1,coly)=3*T(i)^2
    hy1(i-1,coly+1)=2*T(i)
    hy1(i-1,coly+2)=1
    hy1(i-1,coly+3)=0
    hy1(i-1,coly+4)=-3*T(i)^2
    hy1(i-1,coly+5)=-2*T(i)
    hy1(i-1,coly+6)=-1
    hy1(i-1,coly+7)=0
end

%%third part of the matrix boundary conditions
hy2=0*hy1
for i=2:(Ny-1)
    coly=1+(i-2)*4
    hy2(i-1,coly)=6*T(i)
    hy2(i-1,coly+1)=2
    hy2(i-1,coly+2)=0
    hy2(i-1,coly+3)=0
    hy2(i-1,coly+4)=-6*T(i)
    hy2(i-1,coly+5)=-2
    hy2(i-1,coly+6)=0
    hy2(i-1,coly+7)=0
end

hyepts=zeros(2,Coefy)
hyepts(1,1)=6*T(1)
hyepts(1,2)=2*T(1)
hyepts(1,3)=1
hyepts(2, end)=1
hyepts(2, end-1)=2*T(end)
hyepts(2, end-2)=6*T(end)

hy=[hy0;hy1;hy2;hyepts]

%%solve
ABCsy=inv(hy)*By

for i=1:My
   rowy=1+(i-1)*4
   ay=ABCsy(rowy)
   by=ABCsy(rowy+1)
   cy=ABCsy(rowy+2)
   dy=ABCsy(rowy+3)
   t=linspace(T(i),T(i+1),100)
   y=ay*t.^3+by*t.^2+cy*t+dy
   py=plot(t,y,'r')
end


%%% Z

Nz= length(Z)
Mz=Nz-1
Coefz=4*Mz

%%create Z vector
Bz=zeros(2*(Nz-2)+2+(Nz-2)+(Nz-2)+2,1)

Nrowz=2*(Nz-2)+2
hz0=zeros(Nrowz,Coefz)


%%first part of the matrix
for i=1:(Nrowz/2)
    colz=((i-1)*4+1)
    rowz=(i-1)*2+1
    for j=colz:(colz+3)
        hz0(rowz,j)=T(i).^(3-j+colz)
        hz0(rowz+1,j)=T(i+1).^(3-j+colz)
    end
    Bz(rowz,1)=Z(i)
    Bz(rowz+1,1)=Z(i+1)
end

%%second part of the matrix using derivatives

hz1=zeros(Nz-2,Coefz)
for i=2:(Nz-1)
    colz=1+(i-2)*4
    hz1(i-1,colz)=3*T(i)^2
    hz1(i-1,colz+1)=2*T(i)
    hz1(i-1,colz+2)=1
    hz1(i-1,colz+3)=0
    hz1(i-1,colz+4)=-3*T(i)^2
    hz1(i-1,colz+5)=-2*T(i)
    hz1(i-1,colz+6)=-1
    hz1(i-1,colz+7)=0
end

%%third part of the matrix boundary conditions
hz2=0*hz1
for i=2:(Nz-1)
    colz=1+(i-2)*4
    hz2(i-1,colz)=6*T(i)
    hz2(i-1,colz+1)=2
    hz2(i-1,colz+2)=0
    hz2(i-1,colz+3)=0
    hz2(i-1,colz+4)=-6*T(i)
    hz2(i-1,colz+5)=-2
    hz2(i-1,colz+6)=0
    hz2(i-1,colz+7)=0
end

hzepts=zeros(2,Coefz)
hzepts(1,1)=6*T(1)
hzepts(1,2)=2*T(1)
hzepts(1,3)=1
hzepts(2, end)=1
hzepts(2, end-1)=2*T(end)
hzepts(2, end-2)=6*T(end)

hz=[hz0;hz1;hz2;hzepts]

%%solve
ABCsz=inv(hz)*Bz

for i=1:Mz
   rowz=1+(i-1)*4
   az=ABCsz(rowz)
   bz=ABCsz(rowz+1)
   cz=ABCsz(rowz+2)
   dz=ABCsz(rowz+3)
   t=linspace(T(i),T(i+1),100)
   z=az*t.^3+bz*t.^2+cz*t+dz
   pz=plot(t,z,'g')
end


fig2=figure()
set(fig2,'color','white')
set(gca,'fontsize',18)
p4=plot3(X,Y,Z,'b*')
grid on 
hold on

for i=1:Mz
    t=T(i):T(i+1)
    rowx=1+(i-1)*4
   ax=ABCsx(rowx)
   bx=ABCsx(rowx+1)
   cx=ABCsx(rowx+2)
   dx=ABCsx(rowx+3)
   
   x=ax*t.^3+bx*t.^2+cx*t+dx
    
    rowy=1+(i-1)*4
   ay=ABCsy(rowy)
   by=ABCsy(rowy+1)
   cy=ABCsy(rowy+2)
   dy=ABCsy(rowy+3)
  
   y=ay*t.^3+by*t.^2+cy*t+dy
    
    rowz=1+(i-1)*4
   az=ABCsz(rowz)
   bz=ABCsz(rowz+1)
   cz=ABCsz(rowz+2)
   dz=ABCsz(rowz+3)
  
   z=az*t.^3+bz*t.^2+cz*t+dz
   
p5=plot3(x,y,z,'b')
end

grid on
hold on

