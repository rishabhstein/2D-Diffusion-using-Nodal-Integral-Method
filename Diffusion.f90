!COPYRIGHT
!
!This solver is based on Nodal Integral Method for the 2-Dimesional Steady State Diffusion equation. 
!Copyright (C) 2017 Rishabh Prakash Sharma. 
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.!
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>
!\***************************************************************\
 

Program Diffusion
!Do not use very high no. of grids.
!Very high no. of grids leads to high aspect ratio and divergence of scheme.
parameter(n=20, m=20, pt=7)		  
	double precision  Sx(m,n),Sy(m,n),Tx(m,n),Ty(m,n),SxxSyy(m,n),TxxTyy(m,n),Sx2Sy2(m,n),Tx2Ty2(m,n),SxTxSyTy(m,n),SxTx(m,n),SyTy(m,n)
	double precision  A2(m,n),A3(m,n),A4(m,n),A5(m,n),A6(m,n),A7(m,n),A8(m,n),A9(m,n)
	double precision  B2(m,n),B4(m,n),B3(m,n),B5(m,n),B6(m,n),B7(m,n),B8(m,n),B9(m,n)
	double precision  F11(m,n),F12(m,n),F13(m,n),F14(m,n),F15(m,n),F16(m,n),F17(m,n),F18(m,n),F29(m,n),F30(m,n)
	double precision  F21(m,n),F22(m,n),F23(m,n),F24(m,n),F25(m,n),F26(m,n),F28(m,n),Sm1(m,n),Sm2(m,n),Sm3(m,n),Sm4(m,n)
	double precision  dr,dtheta,err5(m+1,n),err6(m+1,n),f(m,n),Tm1(m,n),Tm2(m,n),Tm3(m,n),Tm4(m,n)
	double precision  Ts(m+1,n),Tt(m,n+1),Tc(m+1,n+1),Tst(m,n),q(m,n),rs(m,n),Thetat(m,n)
	double precision  tol5,tol6,temp5(m,n+1),temp6(m+1,n),temp7,temp8,temp9,temp10,temp12(m,n),temp13(m,n)
	double precision  temptol5,temptol6,xdiv,ydiv,tol,s(7),t(7),ws(7),wt(7)
	double precision  x1(m,n),x2(m,n),x3(m,n),x4(m,n),x5(m,n),x6(m,n),x7(m,n),x8(m,n),x9(m,n)
	double precision  y1(m,n),y2(m,n),y3(m,n),y4(m,n),y5(m,n),y6(m,n),y7(m,n),y8(m,n),y9(m,n)
	double precision  xs(m+1,n),ys(m+1,n),xt(m,n+1),yt(m,n+1),FSx(7,7),FSy(7,7),FTx(7,7),FTy(7,7)
	double precision  FSx2Sy2(7,7),FTx2Ty2(7,7),FSxTxSyTy(7,7),FSxTx(7,7),FSyTy(7,7)
	double precision  l1(m,n),l2(m,n),l3(m,n),l4(m,n),m1x(m,n),m1y(m,n),m2x(m,n),m2y(m,n),m3x(m,n),m3y(m,n),m4x(m,n),m4y(m,n)
	
	real	T1, T2,r1,r2,angle
	integer i,j,ks,kt,var
	 
	CALL CPU_TIME(T1)
	
!=================================Parameters============================!	
	r1=1.0;		!Inner radius
	r2=2.0;		!Outer radius
	xdiv=n;
	ydiv=m;
	tol=1.0e-08;	!tolerance

	angle=2*ATAN(1.0);	!Angle=Pi/2
	dr=(r2-r1)/n;		
	dtheta=angle/m;
!-----------------------------------------------------------------------!

!============================Grid generation=============================!
!-----------Co-ordinates for plotting-------------!
do i=1,m+1
    do j=1,n
        xs(i,j)=(r1+(dr/2.0+(j-1)*dr))*(cos((i-1)*dtheta));
        ys(i,j)=(r1+(dr/2.0+(j-1)*dr))*(sin((i-1)*dtheta));
	rs(i,j)=sqrt(xs(i,j)**2+ys(i,j)**2)
    enddo
enddo
do i=1,m
    do j=1,n+1
        xt(i,j)=(r1+(j-1)*dr)*(cos(dtheta/2.0+(i-1)*dtheta));
        yt(i,j)=(r1+(j-1)*dr)*(sin(dtheta/2.0+(i-1)*dtheta));
	Thetat(i,j)=ATAN(yt(i,j)/xt(i,j))
    enddo
enddo
    
!-------Discretization of Inner domain using linear elements--------------!
if (n.EQ.2)  then
        var=2; 
    else
        var=n-1; 
    endif

do i=1,m
    do j=2,var
        x1(i,j)=(r1+(j)*dr)*(cos((i)*dtheta));
        y1(i,j)=(r1+(j)*dr)*(sin((i)*dtheta));
        x2(i,j)=(r1+(j-1)*dr)*(cos((i)*dtheta));
        y2(i,j)=(r1+(j-1)*dr)*(sin((i)*dtheta));
        x3(i,j)=(r1+(j-1)*dr)*(cos((i-1)*dtheta));
        y3(i,j)=(r1+(j-1)*dr)*(sin((i-1)*dtheta));
        x4(i,j)=(r1+(j)*dr)*(cos((i-1)*dtheta));
        y4(i,j)=(r1+(j)*dr)*(sin((i-1)*dtheta));
        x5(i,j)=(x1(i,j)+x2(i,j))/2;
        y5(i,j)=(y1(i,j)+y2(i,j))/2;
        x6(i,j)=(x3(i,j)+x2(i,j))/2;
        y6(i,j)=(y3(i,j)+y2(i,j))/2;
        x7(i,j)=(x3(i,j)+x4(i,j))/2;
        y7(i,j)=(y3(i,j)+y4(i,j))/2;
        x8(i,j)=(x1(i,j)+x4(i,j))/2;
        y8(i,j)=(y1(i,j)+y4(i,j))/2;
        x9(i,j)=(x1(i,j)+x2(i,j)+x3(i,j)+x4(i,j))/4;
        y9(i,j)=(y1(i,j)+y2(i,j)+y3(i,j)+y4(i,j))/4;
    enddo
enddo
!-------------------------------------------------------------------------!

!------------Discretization of curved boundaries using 9-noded quadratic elements---------------!
! For inner radius(r1=1)
do i=1,m
    	j=1;
        x1(i,j)=x2(i,j+1);
        y1(i,j)=y2(i,j+1);
        x2(i,j)=(r1+(j-1)*dr)*(cos((i)*dtheta));
        y2(i,j)=(r1+(j-1)*dr)*(sin((i)*dtheta));
        x3(i,j)=(r1+(j-1)*dr)*(cos((i-1)*dtheta));
        y3(i,j)=(r1+(j-1)*dr)*(sin((i-1)*dtheta));
        x4(i,j)=x3(i,j+1);
        y4(i,j)=y3(i,j+1);
        x6(i,j)=(r1+(j-1)*dr)*(cos(dtheta/2.0+(i-1)*dtheta));
        y6(i,j)=(r1+(j-1)*dr)*(sin(dtheta/2.0+(i-1)*dtheta));
        x7(i,j)=(r1+((dr/2.0)+(j-1)*dr))*(cos((i-1)*dtheta));
        y7(i,j)=(r1+((dr/2.0)+(j-1)*dr))*(sin((i-1)*dtheta));
        x5(i,j)=(r1+((dr/2)+(j-1)*dr))*(cos((i)*dtheta));
        y5(i,j)=(r1+((dr/2)+(j-1)*dr))*(sin((i)*dtheta));
        x8(i,j)=(x1(i,j)+x4(i,j))/2;
        y8(i,j)=(y1(i,j)+y4(i,j))/2;
        x9(i,j)=(x1(i,j)+x2(i,j)+x3(i,j)+x4(i,j)+x5(i,j)+x6(i,j)+x7(i,j)+x8(i,j))/8;
        y9(i,j)=(y1(i,j)+y2(i,j)+y3(i,j)+y4(i,j)+y5(i,j)+y6(i,j)+y7(i,j)+y8(i,j))/8;
    
enddo
!For Outer Radius(r2=2)
do i=1,m
    j=n;
        x1(i,j)=(r1+(j)*dr)*(cos((i)*dtheta));
        y1(i,j)=(r1+(j)*dr)*(sin((i)*dtheta));
        x2(i,j)=x1(i,j-1);
        y2(i,j)=y1(i,j-1);
        x3(i,j)=x4(i,j-1);
        y3(i,j)=y4(i,j-1);
        x4(i,j)=(r1+(j)*dr)*(cos((i-1)*dtheta));
        y4(i,j)=(r1+(j)*dr)*(sin((i-1)*dtheta));
        x6(i,j)=(x2(i,j)+x3(i,j))/2;
        y6(i,j)=(y2(i,j)+y3(i,j))/2;
        x7(i,j)=(r1+((dr/2.0)+(j-1)*dr))*(cos((i-1)*dtheta));
        y7(i,j)=(r1+((dr/2.0)+(j-1)*dr))*(sin((i-1)*dtheta));
        x5(i,j)=(r1+((dr/2)+(j-1)*dr))*(cos((i)*dtheta));
        y5(i,j)=(r1+((dr/2)+(j-1)*dr))*(sin((i)*dtheta));
        x8(i,j)=(r1+(j)*dr)*(cos(dtheta/2.0+(i-1)*dtheta));
        y8(i,j)=(r1+(j)*dr)*(sin(dtheta/2.0+(i-1)*dtheta));
        x9(i,j)=(x1(i,j)+x2(i,j)+x3(i,j)+x4(i,j)+x5(i,j)+x6(i,j)+x7(i,j)+x8(i,j))/8;
        y9(i,j)=(y1(i,j)+y2(i,j)+y3(i,j)+y4(i,j)+y5(i,j)+y6(i,j)+y7(i,j)+y8(i,j))/8;
enddo
!----------------------------------------------------------------------------------------!

!=========================Transformation Coefficicents=======================!
do i=1,m
    do j=1,n
        A2(i,j)=(x8(i,j)-x6(i,j))/2.0;
        A3(i,j)=(x5(i,j)-x7(i,j))/2.0;
        A4(i,j)=(x6(i,j)+x8(i,j))/2.0-x9(i,j);
        A5(i,j)=(x5(i,j)+x7(i,j))/2.0-x9(i,j);
        A6(i,j)=(x1(i,j)-x2(i,j)+x3(i,j)-x4(i,j))/4.0;
        A7(i,j)=(x1(i,j)+x2(i,j)-x3(i,j)-x4(i,j))/4.0-(-x7(i,j)+x5(i,j))/2.0;
        A8(i,j)=(x1(i,j)-x2(i,j)-x3(i,j)+x4(i,j))/4.0-(-x6(i,j)+x8(i,j))/2.0;
        A9(i,j)=(x1(i,j)+x2(i,j)+x3(i,j)+x4(i,j))/4.0-(x5(i,j)+x6(i,j)+x7(i,j)+x8(i,j))/2.0+x9(i,j);
        B2(i,j)=(y8(i,j)-y6(i,j))/2.0;
        B3(i,j)=(y5(i,j)-y7(i,j))/2.0;
        B4(i,j)=(y6(i,j)+y8(i,j))/2.0-y9(i,j);
        B5(i,j)=(y5(i,j)+y7(i,j))/2.0-y9(i,j);
        B6(i,j)=(y1(i,j)-y2(i,j)+y3(i,j)-y4(i,j))/4.0;
        B7(i,j)=(y1(i,j)+y2(i,j)-y3(i,j)-y4(i,j))/4.0-(-y7(i,j)+y5(i,j))/2.0;
        B8(i,j)=(y1(i,j)-y2(i,j)-y3(i,j)+y4(i,j))/4.0-(-y6(i,j)+y8(i,j))/2.0;
        B9(i,j)=(y1(i,j)+y2(i,j)+y3(i,j)+y4(i,j))/4.0-(y5(i,j)+y6(i,j)+y7(i,j)+y8(i,j))/2.0+y9(i,j);
    enddo
enddo
!----------------------------------------------------------------------------!

!==============Length of each cell==========!
do i=1,m
    do j=1,n
        l1(i,j)=sqrt((x1(i,j)-x4(i,j))**2+(y1(i,j)-y4(i,j))**2);
        l2(i,j)=sqrt((x2(i,j)-x1(i,j))**2+(y2(i,j)-y1(i,j))**2);
        l3(i,j)=sqrt((x3(i,j)-x2(i,j))**2+(y3(i,j)-y2(i,j))**2);
        l4(i,j)=sqrt((x4(i,j)-x3(i,j))**2+(y4(i,j)-y3(i,j))**2);
    enddo
enddo
!-------------------------------------------!

!=================Normal unit vectors at each surface of a cell==================!
do i=1,m
    do j=1,n
        m1x(i,j)=((y1(i,j)-y4(i,j)))/l1(i,j);
        m1y(i,j)=-((x1(i,j)-x4(i,j)))/l1(i,j);
        
        m2x(i,j)=((y2(i,j)-y1(i,j)))/l2(i,j);
        m2y(i,j)=-((x2(i,j)-x1(i,j)))/l2(i,j);
        
        m3x(i,j)=((y3(i,j)-y2(i,j)))/l3(i,j);
        m3y(i,j)=-((x3(i,j)-x2(i,j)))/l3(i,j);
        
        m4x(i,j)=((y4(i,j)-y3(i,j)))/l4(i,j);
        m4y(i,j)=-((x4(i,j)-x3(i,j)))/l4(i,j);
    enddo
enddo
!---------------------------------------------------------------------------------!

!=============Derivatives used in transfomation of Equation from Global to local coordinate system=======!
Sx(m,n)=0.0;
Sy(m,n)=0.0;
Tx(m,n)=0.0;
Ty(m,n)=0.0;
Sx2Sy2(m,n)=0.0;
Tx2Ty2(m,n)=0.0;
SxxSyy(m,n)=0.0;
TxxTyy(m,n)=0.0;
SxTxSyTy(m,n)=0.0;
SxTx(m,n)=0.0;
SyTy(m,n)=0.0;

!-----------------Gauss-Quadrature points and corresponding weights for Numerical integration----------------!
! 7 point gauss quadrature is used

s(1)=-0.9491079123  ;
s(2)=-0.7415311856 ;
s(3)=-0.4058451514 ;
s(4)=0.0;
s(5)=0.4058451514 ;
s(6)=0.7415311856 ;
s(7)=0.9491079123 ;
t(1)=-0.9491079123;
t(2)=-0.7415311856 ;
t(3)=-0.4058451514 ;
t(4)=0.0;
t(5)=0.4058451514 ;
t(6)=0.7415311856 ;
t(7)=0.9491079123 ;
ws(1)=0.1294849662 ;
ws(2)=0.2797053915 ;
ws(3)=0.3818300505 ;
ws(4)=0.4179591837 ;
ws(5)=0.3818300505 ;
ws(6)=0.2797053915 ;
ws(7)=0.1294849662; 
wt(1)=0.1294849662 ;
wt(2)=0.2797053915 ;
wt(3)=0.3818300505 ;
wt(4)=0.4179591837 ;
wt(5)=0.3818300505 ;
wt(6)=0.2797053915 ;
wt(7)=0.1294849662 ;
!-----------------------------------------------------------------------------------!
do i=1,m
     do j=1,n
        
		do ks=1,pt
			do kt=1,pt
			FSx(ks,kt)=(B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))/&
                           ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt))))-& 
                           (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))));
                	
			FSy(ks,kt)=(-A3(i,j) - s(ks)*(A6(i,j) + A7(i,j)*s(ks)) - 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))/&
                           ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))) -&
                           (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))));
                	
			FTx(ks,kt)=(-B2(i,j) - 2*B4(i,j)*s(ks) - t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt))))/&
                           ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))) -&
                           (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))));

                	FTy(ks,kt)= (A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt))))/&
                           ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))) -& 
                           (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))));

                 	FSxTx(ks,kt)=((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(-B2(i,j) - 2*B4(i,j)*s(ks) - t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))))/&
                              ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))) -& 
                              (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))))**2; 

                 	FSyTy(ks,kt)=((-A3(i,j) - s(ks)*(A6(i,j) + A7(i,j)*s(ks)) - 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))))/&
                              ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))) -&
                              (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))))**2;         
             enddo
		enddo
		do ks=1,pt
			do kt=1,pt
				temp7=ws(ks)*wt(kt)*FSx(ks,kt);
				Sx(i,j)=Sx(i,j)+temp7;
				temp8=ws(ks)*wt(kt)*FSy(ks,kt);
				Sy(i,j)=Sy(i,j)+temp8;
				temp9=ws(ks)*wt(kt)*FTx(ks,kt);
				Tx(i,j)=Tx(i,j)+temp9;
				temp10=ws(ks)*wt(kt)*FTy(ks,kt);
				Ty(i,j)=Ty(i,j)+temp10;
				temp12(i,j)=ws(ks)*wt(kt)*FSxTx(ks,kt);
				SxTx(i,j)=SxTx(i,j)+temp12(i,j);
				temp13(i,j)=ws(ks)*wt(kt)*FSyTy(ks,kt);
				SyTy(i,j)=SyTy(i,j)+temp13(i,j);
			enddo
		enddo
     
	 
	 enddo
enddo


Sx=Sx/4.0;
Sy=Sy/4.0;
Tx=Tx/4.0;
Ty=Ty/4.0;
SxxSyy=SxxSyy/4.0;
TxxTyy=TxxTyy/4.0

do i=1,m
     do j=1,n
        do ks=1,pt
		do  kt=1,pt

                FSx2Sy2(ks,kt)=((A3(i,j)+s(ks)*(A6(i,j)+A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))**2+ (B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))**2)/&
                               ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))) -& 
                               (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))))**2;
                FTx2Ty2(ks,kt)=((A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt))))**2+ (B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt))))**2)/&
                               ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))) -& 
                               (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))))**2;
                 FSxTxSyTy(ks,kt)=(-((A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt))))) -& 
                                 (B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))))/&
                                 ((B3(i,j) + s(ks)*(B6(i,j) + B7(i,j)*s(ks)) + 2*(B5(i,j) + s(ks)*(B8(i,j) + B9(i,j)*s(ks)))*t(kt))*(A2(i,j) + 2*A4(i,j)*s(ks) + t(kt)*(A6(i,j) + A8(i,j)*t(kt) + 2*s(ks)*(A7(i,j) + A9(i,j)*t(kt)))) -& 
                                 (A3(i,j) + s(ks)*(A6(i,j) + A7(i,j)*s(ks)) + 2*(A5(i,j) + s(ks)*(A8(i,j) + A9(i,j)*s(ks)))*t(kt))*(B2(i,j) + 2*B4(i,j)*s(ks) + t(kt)*(B6(i,j) + B8(i,j)*t(kt) + 2*s(ks)*(B7(i,j) + B9(i,j)*t(kt)))))**2;
             enddo
        enddo
        do ks=1,pt
            do kt=1,pt
                temp7=ws(ks)*wt(kt)*FSx2Sy2(ks,kt);
                Sx2Sy2(i,j)=Sx2Sy2(i,j)+temp7;
                temp8=ws(ks)*wt(kt)*FTx2Ty2(ks,kt);
                Tx2Ty2(i,j)=Tx2Ty2(i,j)+temp8;
                temp11=ws(ks)*wt(kt)*FSxTxSyTy(ks,kt);
                SxTxSyTy(i,j)=SxTxSyTy(i,j)+temp11;
            enddo
        enddo
     enddo
enddo
Sx2Sy2=Sx2Sy2/4.0;
Tx2Ty2=Tx2Ty2/4.0;
SxTxSyTy=SxTxSyTy/4.0;
!------------------------------------------------------------------------------------------------------!


!================Direction Cosine for each surface of a cell=================!
do i=1,m
    do j=1,n
        Sm1(i,j)=(m1x(i,j)*Sx(i,j)+m1y(i,j)*Sy(i,j)); ! For edge 1
        Sm2(i,j)=(m2x(i,j)*Sx(i,j)+m2y(i,j)*Sy(i,j)); ! For edge 2
        Sm3(i,j)=-(m3x(i,j)*Sx(i,j)+m3y(i,j)*Sy(i,j)); ! For edge3
        Sm4(i,j)=-(m4x(i,j)*Sx(i,j)+m4y(i,j)*Sy(i,j)); ! For edge4
        Tm1(i,j)=(m1x(i,j)*Tx(i,j)+m1y(i,j)*Ty(i,j)); ! For edge 1
        Tm2(i,j)=(m2x(i,j)*Tx(i,j)+m2y(i,j)*Ty(i,j)); ! For edge 2
        Tm3(i,j)=-(m3x(i,j)*Tx(i,j)+m3y(i,j)*Ty(i,j)); ! For edge3
        Tm4(i,j)=-(m4x(i,j)*Tx(i,j)+m4y(i,j)*Ty(i,j)); ! For edge4
    enddo
enddo
!-----------------------------------------------------------------------------!

!=================Coefficients of final discrete equation=================!
do i=1,m-1
    do j=1,n
        F12(i,j)=(1*(Tm2(i,j))/2+1*(Tm4(i+1,j))/2+(3*(Tm2(i,j))/2)*(Sx2Sy2(i,j))/(((Sx2Sy2(i,j))+(Tx2Ty2(i,j))))+(3*(Tm4(i+1,j))/2)*(Sx2Sy2(i+1,j))/(((Sx2Sy2(i+1,j))+Tx2Ty2(i+1,j))));
    enddo
enddo
do i=1,m
    do j=1,n
        F11(i,j)=(1*(Tm4(i,j))/2-(3*Tm4(i,j)/2)*(Sx2Sy2(i,j))/((Sx2Sy2(i,j)+Tx2Ty2(i,j))));
        F13(i,j)=(1*(Tm2(i,j))/2-(3*Tm2(i,j)/2)*(Sx2Sy2(i,j))/((Sx2Sy2(i,j)+Tx2Ty2(i,j))));
        F14(i,j)=(3*(Tm2(i,j))/2)*(Sx2Sy2(i,j))/(((Sx2Sy2(i,j))+(Tx2Ty2(i,j))));
        F15(i,j)=(3*(Tm4(i,j))/2)*(Sx2Sy2(i,j))/(((Sx2Sy2(i,j))+Tx2Ty2(i,j)));
    enddo
enddo
do i=1,m
    do j=1,n
        F16(i,j)=1*(Tm2(i,j))/(((Sx2Sy2(i,j))+(Tx2Ty2(i,j))));
        F26(i,j)=1*(Tm4(i,j))/(((Sx2Sy2(i,j))+(Tx2Ty2(i,j))));
        F17(i,j)=2*(SxTxSyTy(i,j));
        F18(i,j)=(1/2+(3/2)*(Sx2Sy2(i,j))/((Sx2Sy2(i,j))+(Tx2Ty2(i,j))));
    enddo
enddo
do i=1,m
    do j=1,n-1
        F22(i,j)=(1*(Sm1(i,j))/2+1*(Sm3(i,j+1))/2+(3*(Sm1(i,j))/2)*(Tx2Ty2(i,j))/(((Sx2Sy2(i,j))+(Tx2Ty2(i,j))))+(3*(Sm3(i,j+1))/2)*(Tx2Ty2(i,j+1))/(((Sx2Sy2(i,j+1))+Tx2Ty2(i,j+1))));
    enddo
enddo
do i=1,m
    do j=1,n
        F21(i,j)=(1*(Sm3(i,j))/2-(3*Sm3(i,j))/2*(Tx2Ty2(i,j))/(((Sx2Sy2(i,j))+Tx2Ty2(i,j))));
        F23(i,j)=(1*(Sm1(i,j))/2-(3*Sm1(i,j))/2*(Tx2Ty2(i,j))/(((Sx2Sy2(i,j))+Tx2Ty2(i,j))));
        F24(i,j)=(3*(Sm1(i,j))/2)*(Tx2Ty2(i,j))/(((Sx2Sy2(i,j))+Tx2Ty2(i,j)));
        F25(i,j)=(3*(Sm3(i,j))/2)*(Tx2Ty2(i,j))/(((Sx2Sy2(i,j))+Tx2Ty2(i,j)));
        F28(i,j)=(1/2+(3/2)*(Tx2Ty2(i,j))/((Sx2Sy2(i,j))+(Tx2Ty2(i,j))));
        F29(i,j)=(1/2-(3/2)*(Tx2Ty2(i,j))/((Sx2Sy2(i,j))+(Tx2Ty2(i,j))));
        F30(i,j)=((3/2)*(Tx2Ty2(i,j)))/(((Sx2Sy2(i,j))+Tx2Ty2(i,j)));
    enddo
enddo
!---------------------------------------------------------------------------------!

!=================Initialization of Main loop==================!
tol5=1 
tol6=1

do while ((tol5.GT.tol).OR.(tol6.GT.tol))
        

!----------------calculation do cross derivatives-----!
            Tc(1,1)=(Tt(1,1)+Ts(1,1))/2.0;
            Tc(1,n+1)=(Ts(1,n)+Tt(1,n+1))/2.0;
            Tc(m+1,1)=(Tt(m,1)+Ts(m+1,1))/2.0;
            Tc(m+1,n+1)=(Tt(m,n+1)+Ts(m+1,n))/2.0;
			
			do j=2,n
				Tc(1,j)=(Ts(1,j-1)+Ts(1,j))/2.0;
				Tc(m+1,j)=(Ts(m+1,j-1)+Ts(m+1,j))/2.0;
			enddo
			do i=2,m
				Tc(i,1)=(Tt(i-1,1)+Tt(i,1))/2.0;
				Tc(i,n+1)=(Tt(i-1,n+1)+Tt(i,n+1))/2.0;
			enddo
            
            		do i=2,m
                		do j=2,n
                    			Tc(i,j)=(Ts(i,j-1)+Ts(i,j)+Tt(i-1,j)+Tt(i,j))/4.0;
                		enddo
            		enddo
            		do i=1,m
                		do j=1,n
                    			Tst(i,j)=((Tc(i+1,j+1)-Tc(i,j+1)-Tc(i+1,j)+Tc(i,j))/4.0);
                    			f(i,j)=q(i,j)+(F17(i,j)*Tst(i,j));
                		enddo
            		enddo
!-----------------------------------------------------!		

!============================Boundary Conditions=========================!
		!Neumann boundary condition at Inner radius of the domain
    		do i=1,m
			Tt(i,1)=(F29(i,1)*Tt(i,2)+F30(i,1)*(Ts(i,1)+Ts(i+1,1))+0.5*(Tc(i+1,1)-Tc(i,1))*(Tm3(i,1)/Sm3(i,1)))/F28(i,1);
		enddo
		
		!Sinosuidal function is applied at outer radius 
		!Frequency of the applied sinosuidal boundary condition can be changed
		!by changing the "feq" parameter given below
		feq=2
		do i=1,m
			Tt(i,n+1)=COS(feq*(4*ATAN(1.0))*ATAN(yt(i,n+1)/xt(i,n+1))/angle);
		enddo
	
		!Dirichlet B.C with "zero" value at Theta=0deg 
        	do j=1,n
			Ts(1,j)=0;
        	enddo
		!Dirichlet B.C with "zero" value at Theta=90deg 
        	do j=1,n
        		Ts(m+1,j)=0;
        	enddo
!-------------------------------------------------------------------------!
!=============Loop for temperature(both s-averaged(Ts) and t-averaged(Tt))================!
         do i=1,m-1
         	do j=1,n-1
                Ts(i+1,j)=(F11(i1+1,j)*Ts(i+2,j)+F13(i,j)*Ts(i,j)+F14(i,j)*(Tt(i,j)+Tt(i,j+1))+F15(i+1,j)*(Tt(i+1,j)+Tt(i+1,j+1))-(f(i,j)*F16(i,j)+f(i+1,j)*F16(i+1,j))+((Tc(i+1,j+1)-Tc(i+1,j))/2)*(Sm4(i+1,j)-Sm2(i,j)))/F12(i,j);
            	Ts(i+1,n)=(F11(i+1,n)*Ts(i+2,n)+F13(i,n)*Ts(i,n)+F14(i,n)*(Tt(i,n)+Tt(i,n+1))+F15(i+1,n)*(Tt(i+1,n)+Tt(i+1,n+1))-(f(i,n)*F16(i,n)+f(i+1,n)*F16(i+1,n))+((Tc(i+1,n+1)-Tc(i+1,n))/2)*(Sm4(i+1,n)-Sm2(i,n)))/F12(i,n);
           	Tt(i,j+1)=(F21(i,j+1)*Tt(i,j+2)+F23(i,j)*Tt(i,j)+F24(i,j)*(Ts(i,j)+Ts(i+1,j))+F25(i,j+1)*(Ts(i,j+1)+Ts(i+1,j+1))-(f(i,j)*F26(i,j)+f(i,j+1)*F26(i,j+1))+((Tc(i+1,j+1)-Tc(i,j+1))/2)*(Tm3(i,j+1)-Tm1(i,j)))/F22(i,j);
           	Tt(m,j+1)=(F21(m,j+1)*Tt(m,j+2)+F23(m,j)*Tt(m,j)+F24(m,j)*(Ts(m,j)+Ts(m+1,j))+F25(m,j+1)*(Ts(m,j+1)+Ts(m+1,j+1))-(f(m,j)*F26(m,j)+f(m,j+1)*F26(m,j+1))+((Tc(m+1,j+1)-Tc(m,j+1))/2)*(Tm3(m,j+1)-Tm1(m,j)))/F22(m,j);
        	enddo
        enddo
!-----------------------------------------------------------------------------------------!

!=======================Residual==========================!
			do i=1,m
				do j=1,n+1
					err5(i,j)=(Tt(i,j)-temp5(i,j));
				enddo
			enddo
			do i=1,m+1
				do j=1,n
					err6(i,j)=(Ts(i,j)-temp6(i,j));
				enddo
			enddo
			do i=1,m
				do j=1,n+1
					temp5(i,j)=Tt(i,j);
				enddo
			enddo
			do i=1,m+1
				do j=1,n
					temp6(i,j)=Ts(i,j);
				enddo
			enddo
				
			temptol5=0.0
			temptol6=0.0
			do i=1,m
				do j=1,n+1
					temptol5=temptol5+err5(i,j)**2;
				enddo
			enddo
			tol5=sqrt(temptol5)/sqrt(xdiv*ydiv);
			do i=1,m+1
				do j=1,n
					temptol6=temptol6+err6(i,j)**2;
				enddo
			enddo
			tol6=sqrt(temptol6)/sqrt(xdiv*ydiv);

			print *,'Error in Tt is ', tol5
			print *,'Error in Ts is ', tol6
!-------------------------------------------------------------!
	   
	enddo !end of while loop
	CALL CPU_TIME(T2)
!----------------------------End of Main loop--------------------------------------!
!================Output Files===============!
!S-averged surface
open(1,file='Xs')
do i=1,m+1
write(1,*) xs(i,:)
enddo
close(1)

open(2,file='Ys')
do i=1,m+1
write(2,*) ys(i,:)
enddo
close(2)

open(3,file='Ts')
do i=1,m+1
write(3,*) Ts(i,:)
enddo
close(3)

!T-averaged surface
open(4,file='Xt')
do i=1,m
write(4,*) xt(i,:)
enddo
close(4)

open(5,file='Yt')
do i=1,m
write(5,*) yt(i,:)
enddo
close(5)

open(6,file='Tt')
do i=1,m
write(6,*) Tt(i,:)
enddo
close(6)

!Temperature profile at Theta=Pi/4
open(7,file='T_P1')
do j=1,n
write(7,*) rs((m/2)+1,j),Ts((m/2)+1,j)
enddo
close(7)

!Temperature profile at r=1.5
open(8,file='T_P2')
do i=1,m
write(8,*) Thetat(i,(n/2)+1),Tt(i,(n/2)+1)
enddo
close(8)

!-------------------------------------------!
 end program       
