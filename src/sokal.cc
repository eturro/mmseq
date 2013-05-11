/*
 *  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
 *  Copyright (c) 2003       David Hastie  <david.hastie@bristol.ac.uk>
 *  Copyright (c) 2003       Peter Green   <p.j.green@bristol.ac.uk>
 *
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 2, as
 *  published by the Free Software Foundation.
 *
 *  BGX is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

#include <math.h>
#include <iostream>

#include "sokal.hh"

int fft(double* xreal, double* ximag, int length);

// Estimate integrated autocorrelation time using the method of Sokal.
// Note that this is actually the sum from -infty to +infty of the acf, hence 
// twice Sokal's definition.  The implementation follows that of 
// Green and Han (1990).
// WARNING: The input array x is destroyed in the process.
int sokal(int *n_, double *x, double *var, double *tau, int* m)
{
  int n=*n_;
  if(n>2<<20)
    {
      std::cerr << "Auto-correlation length exceeded" << std::endl;
      return 100;
    }

  // Declare and initialize the imaginary part of x
  double *x_im = new double[n];
  for(int i=0; i<n; ++i)
    x_im[i]=0;

  // Use FFT to compute autocorrelations (returned in x), and 
  // variances (returned in var).
  int test=fft(x,x_im,n);
  if(test!=0) return test;

  for(int i=0; i<n; ++i)
    {
      x[i] = x[i]*x[i] + x_im[i]*x_im[i];
      x_im[i] = 0;
    }
  *x=0;

  test=fft(x,x_im,n);
  delete [] x_im;
  if(test!=0) return test;

  *var = *x/((double)n*(n-1));

  double c = 1.0/ *x;
  for(int i=0; i<n; ++i)
    x[i] *= c;

  // Use Sokal's adaptive truncated periodogram method to estimate 
  // integrated autocorrelation time (returned in tau).
  // The number of terms used is returned in m.  If m=n+1 then the 
  // algorithm did not converge before all terms of the acf were used up.
  double sum = -0.333333333333333333333;
  *m=n+1;
  for(int i=0; i<n; ++i)
    {
      sum += x[i] - 0.166666666666666666666;
      if(sum<0)
	{
	  *m=i+1;
	  break;
	}
    }
  *tau = 2*(sum+(*m-1.0)/6.0);

  return 0;
}

// An fft routine I (Graeme Ambler) nabbed from David Hastie.  He in turn 
// got it from Peter Green, who got it from Erik Renshaw (I think) ... 
// There were no copyright notices on David or Peter's code, so I have added 
// David and Peter to the copyright list at the top of this file to get 
// the attribution as correct as I can.  Peter's code was in f77.  David 
// translated that to C and I added C++ iostreams support, a couple of 
// tweaks and a bugfix or two.
int fft(double* xreal, double* ximag, int length)
{
  double z,bcos,bsin,temp,cw1,sw1;
  double cw2=0; double cw3=0; double sw2=0; double sw3=0;
  double xs0,xs1,xs2,xs3,ys0,ys1,ys2,ys3,x1,x2,x3,y1,y2,y3;

  int i0,i1,i2,i3,ul[20],n,counta,countb,time,test;
  int indic=0;
  int j0,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17,j18,j19;

  n = length<0 ? -length : length;

  if(n<4){
    std::cerr << "Fast Fourier length too short" << std::endl;
    return 200;
  }
  
  if(length<0){
    for(i1=0;i1<n;i1++){
      ximag[i1]=-ximag[i1];
    }
  }

  test=n;
  while(test>1){
    if(test&1){
      std::cerr << "Fast Fourier length must be power of 2" << std::endl;
      return 201;
    }
    test>>=1;
  }
    
  counta=n/4;
  time=0;
  
  while(counta>0){
    countb=counta*4;
    time+=2;
    
    z=M_PI/countb;
    double sin_z=sin(z);
    bcos=-2.0*sin_z*sin_z;
    bsin=sin(2.0*z);
    
    cw1=1.0;
    sw1=0.0;
    for(j1=0;j1<counta;j1++){
      for(j2=j1;j2<n;j2+=countb){
	i0=j2;
	i1=i0+counta;
	i2=i1+counta;
	i3=i2+counta;
	xs0=xreal[i0]+xreal[i2];
	xs1=xreal[i0]-xreal[i2];
	ys0=ximag[i0]+ximag[i2];
	ys1=ximag[i0]-ximag[i2];
	xs2=xreal[i1]+xreal[i3];
	xs3=xreal[i1]-xreal[i3];
	ys2=ximag[i1]+ximag[i3];
	ys3=ximag[i1]-ximag[i3];

	xreal[i0]=xs0+xs2;
	ximag[i0]=ys0+ys2;
	
	x1=xs1+ys3;
	y1=ys1-xs3;
	x2=xs0-xs2;
	y2=ys0-ys2;
	x3=xs1-ys3;
	y3=ys1+xs3;

	if(j1==0){
	  xreal[i2]=x1;
	  ximag[i2]=y1;
	  xreal[i1]=x2;
	  ximag[i1]=y2;
	  xreal[i3]=x3;
	  ximag[i3]=y3;
	}
	else{
	  xreal[i2]=x1*cw1+y1*sw1;
	  ximag[i2]=y1*cw1-x1*sw1;
	  xreal[i1]=x2*cw2+y2*sw2;
	  ximag[i1]=y2*cw2-x2*sw2;
	  xreal[i3]=x3*cw3+y3*sw3;
	  ximag[i3]=y3*cw3-x3*sw3;
	}
      }
      if(j1<(counta-1)){
	z=cw1*bcos-sw1*bsin+cw1;
	sw1=bcos*sw1+bsin*cw1+sw1;
	temp=1.5-0.5*(z*z+sw1*sw1);
	cw1=z*temp;
	sw1=sw1*temp;
	cw2=cw1*cw1-sw1*sw1;
	sw2=2.0*cw1*sw1;
	cw3=cw1*cw2-sw1*sw2;
	sw3=cw1*sw2+cw2*sw1;
      }
    }

    indic=0;
    if(counta>1){
      counta/=4;
      indic=1;
    }
    else{
      counta=0;
    }
  }
  if(indic){
    for(j1=0;j1<n;j1+=2){
      temp=xreal[j1]+xreal[j1+1];
      xreal[j1+1]=xreal[j1]-xreal[j1+1];
      xreal[j1]=temp;
      temp=ximag[j1]+ximag[j1+1];
      ximag[j1+1]=ximag[j1]-ximag[j1+1];
      ximag[j1]=temp;
    }
    time++;
  }
  
  if(length<0){
    for(j1=0;j1<n;j1++){
      ximag[j1]=-ximag[j1];
    }
  }

  if(time>20){
    std::cerr << "Fast Fourier length too long" << std::endl;
    return 202;
  }

  i1=20-time;
  for(j1=0;j1<i1;j1++){
    ul[j1]=1;
  }
  if(i1==0){
    ul[0]=2;
    i1++;
  }
  for(j1=i1;j1<20;j1++){
    ul[j1]=ul[j1-1]<<1;
  }

  i0=0;
  for(j0=0;j0<ul[0];j0++){
    for(j1=j0;j1<ul[1];j1+=ul[0]){
      for(j2=j1;j2<ul[2];j2+=ul[1]){
	for(j3=j2;j3<ul[3];j3+=ul[2]){
	  for(j4=j3;j4<ul[4];j4+=ul[3]){
	    for(j5=j4;j5<ul[5];j5+=ul[4]){
	      for(j6=j5;j6<ul[6];j6+=ul[5]){
		for(j7=j6;j7<ul[7];j7+=ul[6]){
		  for(j8=j7;j8<ul[8];j8+=ul[7]){
		    for(j9=j8;j9<ul[9];j9+=ul[8]){
		      for(j10=j9;j10<ul[10];j10+=ul[9]){
			for(j11=j10;j11<ul[11];j11+=ul[10]){
			  for(j12=j11;j12<ul[12];j12+=ul[11]){
			    for(j13=j12;j13<ul[13];j13+=ul[12]){
			      for(j14=j13;j14<ul[14];j14+=ul[13]){
				for(j15=j14;j15<ul[15];j15+=ul[14]){
				  for(j16=j15;j16<ul[16];j16+=ul[15]){
				    for(j17=j16;j17<ul[17];j17+=ul[16]){
				      for(j18=j17;j18<ul[18];j18+=ul[17]){
					for(j19=j18;j19<ul[19];j19+=ul[18]){
					  if((i0-j19)<0){
					    temp=xreal[i0];
					    xreal[i0]=xreal[j19];
					    xreal[j19]=temp;
					    temp=ximag[i0];
					    ximag[i0]=ximag[j19];
					    ximag[j19]=temp;
					  }
					  i0++;
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return 0;
}
