package mainpackage;

import static java.util.concurrent.ThreadLocalRandom.current;
import org.jtransforms.fft.*;
import java.util.ArrayList;

public class fastlomb {
	public double[] px;
	public double[] f;

public void FastLomb(double[] x){
	double sum=0;
	double[] t=new double[x.length];
	for (int i=0;i<x.length;i++){
		sum=sum+x[i];
		t[i]=sum/1000;
	}
	double[]spectrum = null;
	int ofac  = 4;
    int hifac = 1;
    double T = t[t.length-1]- t[0];
    
    int n = t.length;
    double Ts=T/(n-1);
  //  System.out.println(Ts);
   int nout = (int) Math.round(0.5*ofac*hifac*n);
   //System.out.println(nout);
   double[] fr=new double[nout];
   for (int i =0; i<nout; i++){
        fr[i]=   i/(n*Ts*ofac);
   }
   double fMaxChnl= max(fr);
   double meanx=getAverage(x);
   for(int i=0;i<n;i++){
	   x[i]=x[i]-meanx;
	  // System.out.println(x[i]);
   }
   int p = nextpow2(x.length)-1;
   
   int MACC=4;
   //System.out.println(nextpow2(ofac*hifac*n*MACC));
   int nfreq  = nextpow2(ofac*hifac*n*MACC); 
   int fndim = 2*nfreq;

   //% Initialize
   double[] wk11 = new double[fndim];
   double[] wk22 = new double[fndim];
   double fac = fndim/(n*Ts*ofac);
  
   double[] nfac = {1,1,2,6,24,120,720,5040,40320,362880};
   double nden = nfac[MACC-1];
   
   double tmin=t[0];
   //System.out.println(tmin);
  
   for(int j=0;j<n;j++){
	    double ck  = 1 + ((t[j]-tmin)*fac)%fndim;
	   
	    //System.out.println(x[j]);
	    double ckk = 1 + (2*(ck-1))%fndim;
	    //System.out.println("ckk");
	    //System.out.println(ckk);
	    wk11 = spread(x[j],wk11,fndim,ck ,MACC,nden);
	    wk22 = spread(1,wk22,fndim,ckk,MACC,nden);
   }
   double[] wk1 = new double[fndim];
   double[] wk2 = new double[fndim];
   for(int j=0;j<fndim-1;j++){
	   wk1[j]=wk11[j+1];
	   //System.out.println(wk1[j]);
	   wk2[j]=wk22[j+1];
	   //System.out.println(wk2[j]);
   }
   //System.out.println(wk1.length);
   DoubleFFT_1D fft=new DoubleFFT_1D(wk1.length);
   double[] a=perparefft(wk1);
   fft.complexForward(a);
   double[] rwk1=new double[nout];
   double[] iwk1=new double[nout];
   //System.out.println(nout);
   //System.out.println(nfreq);
   for(int i=1;i<nout;i++){
	   rwk1[i-1]=a[2*i]; 
	   //System.out.println(i);
	
	   iwk1[i-1]=a[2*i+1];
	   //System.out.println(iwk1[i]);
   }
   DoubleFFT_1D fft2=new DoubleFFT_1D(wk2.length);
   double[] b=perparefft(wk2);
   fft2.complexForward(b);
   double[] rwk2=new double[nout];
   double[] iwk2=new double[nout];
   double[] hypo=new double[nout];
   for(int i=1;i<nout;i++){
	  // System.out.println("rwk2");
	   hypo[i-1]=Math.sqrt(Math.pow(b[2*i],2)+Math.pow(b[2*i+1],2));
	  // System.out.println(hypo[i-1]);
	   rwk2[i-1]=b[2*i]; 
	  // System.out.println(rwk2[i-1]);
	   
	   iwk2[i-1]=b[2*i+1];
	 // System.out.println( iwk2[i-1]);
   }
   double[] hc2wt=new double[nout];
   double[] hs2wt=new double[nout];
   double[] swt=new double[nout];
   double[] cwt=new double[nout];
   double[] den=new double[nout];
   double[]  cterm=new double[nout];
   double[] sterm=new double[nout];
   
   double[] wwk2=new double[nout-1];
   for(int i=0;i<nout-1;i++){
	   hc2wt[i]= 0.5*rwk2[i]/hypo[i];
	   //System.out.println( "swt");
	   //System.out.println( hc2wt[i]);
	   hs2wt[i] = 0.5*iwk2[i]/hypo[i];
	   cwt[i]   = Math.sqrt(0.5+hc2wt[i]);
	   swt[i]   = Math.signum(hs2wt[i])*(Math.sqrt(0.5-hc2wt[i]));
	   //System.out.println( swt[i]);
	   den[i]=0.5*n + hc2wt[i]*rwk2[i] + hs2wt[i]*iwk2[i];
	   cterm[i] = (cwt[i]*rwk1[i] + swt[i]*iwk1[i])*(cwt[i]*rwk1[i] + swt[i]*iwk1[i])/den[i];
	   sterm[i] = (cwt[i]*iwk1[i] - swt[i]*rwk1[i])*(cwt[i]*iwk1[i] - swt[i]*rwk1[i])/(n-den[i]);
	   wwk2[i]= cterm[i]+sterm[i];
	  // System.out.println( "ai");
	   //System.out.println( wwk2[i]);
   }
   double[] px=new double[nout];
   
  for (int i=0;i<nout-1;i++)
   {
   px[i]=wwk2[i]*Ts;
   //System.out.println(String.valueOf(i));
   //System.out.println(String.valueOf(fr[i]));
   //System.out.println(String.valueOf(px[i]));
   }
this.f=fr;
this.px=px;
}
public double[] getp(){
	return this.px;
}
public double[] getf(){
	return this.f;
}
public static double getAverage(double []array){
    int sum = 0;
    
    for(int i = 0;i < array.length;i++){
        sum += array[i];
    }
    return (double)(sum / array.length);
}
public static int nextpow2(int n) {
if(n<=0)return 1;
else
{

 if((n & (n-1)) == 0)
	return   n;
 else{
	 int i=0;
	 int t=1;
     while((t=(int) Math.pow(2,i))<=n){
                i++;}
    return t;
}
}
}

public static double max(double []arr){ //（这个函数getMax后面的形式参数为什么不 能换成 arr ）
double max = arr[0]; //int []arr已经定义了是数组，那么直接写arr和 再输入 int []arr 有何不同啊
for(int i=0; i< arr.length;i++){
if(arr[i]>max){
max=arr[i];
}
}
return max;
}
public static double min(double []arr){
double min = arr[0]; //int []arr已经定义了是数组，那么直接写arr和 再输入 int []arr 有何不同啊
for(int i=0; i< arr.length;i++){
if(arr[i]<min){
min=arr[i];
}
}
return min;
}
public static double prod(double []arr){
	double prod=1;
	for (int i=0;i<arr.length;i++){
		prod=prod*arr[i];
	}
	return prod;
}
public static double[] spread(double y, double[] yy,int n,double x, int m, double nden){
	if(x==Math.round(x))
		yy[(int) x]=yy[(int) x]+y;
	else
	{
		int ilo=(int)Math.min(Math.max(Math.floor(x-0.5*m+1), 1),n-m+1);
		int ihi = ilo+m-1;
		double[] inew = new double[m-1];
		for(int j=ilo+1;j<ihi+1;j++){
       	inew[j-ilo-1]=x-j;
		}
		double fac = (x-ilo)*prod(inew);
		yy[ihi-1] = yy[ihi-1] + y*fac/(nden*(x-ihi));
		 for (int j = ihi-1;j>=ilo;j--){
			        nden = (nden/(j+1-ilo))*(j-ihi);
			        yy[j-1] = yy[j-1] + y*fac/(nden*(x-j));			        
			        }
			  
	}
	return yy;
}
public static double[] perparefft(double[] x){
	double[] a=new double[2*x.length];
	for (int i=0;i<x.length;i++){
		a[2*i]=x[i];
		a[2*i+1]=0;
		
	}
	return a;
}

public static void main(String[] args){
	double [] rr={400,1700,1250,950,350,550,1600,900,900,1050,600,1300,700,850,1700,550,1200,300,500,350,250,600,650,650};
	//double [] t={0.400000000000000,2.10000000000000,3.35000000000000,4.30000000000000,4.65000000000000,5.20000000000000,6.80000000000000,7.70000000000000,8.60000000000000,9.65000000000000,10.2500000000000,11.5500000000000,12.2500000000000,13.1000000000000,14.8000000000000,15.3500000000000,16.5500000000000,16.8500000000000,17.3500000000000,17.7000000000000,17.9500000000000,18.5500000000000,19.2000000000000,19.8500000000000};
	fastlomb ff=new fastlomb(); 
	ff.FastLomb(rr);
	
double[] f=ff.getf();
double[] p=ff.getp();
	double pfp=0;

	for (int i=0;i<f.length;i++)
	{
		if(f[i]<=0.40&f[i]>0)
			pfp=pfp+p[i]*f[i];
			
		else{
			
		}
	}
	
	double lnpfp=Math.log(pfp);
	 //System.out.println(String.valueOf(lnpfp));
	}
}


