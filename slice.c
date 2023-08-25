#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#include "martin.h"
#define BSIZE 28800000
#define _XOPEN_SOURCE 700

/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;
XImage *x_image;
unsigned char *x_buffer;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();


void disp (unsigned char *,int,int);

void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}


int zoom_it (float *in, float *out, float zoom, float xc, float yc)
{
	int x,y,p,pr,pzi,xzi;
	float xz,yz,pz,qz,dy,dx;

	for (y=0;y<Y_SIZE;y++)
	{
		p=y*3*X_SIZE;
		pz=((float)y-yc)*zoom;
		pz+=yc;
		pzi=pz;
		dy=pz-pzi;
		if (pzi+1>=Y_SIZE || pzi<0){ continue;} 
		pzi*=(3*X_SIZE);
		for (x=0;x<X_SIZE;x++)
		{
			xz=((float)x-xc)*zoom;
			xz+=xc;
			xzi=xz;
			dx=xz-xzi;
			//printf ("dx %f dy %f \n",dx,dy);
			if (xzi+1>=X_SIZE || xzi< 0){ continue;} 
			int i;
			for (i=0;i<3;i++)
			{
				float a,b,c,d,x1,x2;
				a=in[pzi+(xzi*3)+i];
				b=in[pzi+((xzi+1)*3)+i];

				c=in[(pzi+(3*X_SIZE))+(xzi*3)+i];
				d=in[(pzi+(3*X_SIZE))+((xzi+1)*3)+i];


				x1=((1-dx)*a)+(dx*b);
				x2=((1-dx)*c)+(dx*d);

				out[p+(x*3)+i]=(x1*(1-dy))+(x2*dy);
			}
		}
	}
}

int add_wav_cur (short *wav, float *image, float xp, float yp, int slice,float max)
{
        int sx,sy,w_len,wp,i,xl,yl,xr,yr,size,totl,totr,l,r,thresh;
        float s,rad,span,xp1,yp1;

        s=(float)slice/18000;
        rad=1500;span=200;
        w_len=44100*2/30;
        thresh=w_len*100;
        wp=slice*w_len;

        size=1; totl=0; totr=0;
        for (i=0;i<w_len/2;i++)
        {
                l=wav[wp+(2*i)];
                r=wav[wp+(2*i)+1];
                if (l>0){ totl+=l;}else{totl-=l;}
                if (r>0){ totr+=r;}else{totr-=r;}
        }
        if (totl<thresh){totl=0;}
        if (totr<thresh){totr=0;}

	int circ;

	circ=(Y_SIZE/2);

        for (i=0;i<w_len/2;i++)
        {
                float left,right,xlf,xlt,ylt,ylf,xrf,xrt,yrf,yrt,phi;

		// half circle
		phi=(-M_PI/4)+((float)i*3*M_PI/(float)(w_len));

                left=wav[wp+(2*i)]*span/24768;
                right=wav[wp+(2*i)+1]*span/32768;

                xlf=left+(s*X_SIZE/2)+span;
                xlt=(X_SIZE/2)-((sin(phi)*(left+circ)));
		xl=((1-s)*xlf)+(s*xlt);

                ylf=i*Y_SIZE*2/w_len;
                ylt=(Y_SIZE/2)-((cos(phi)*(left+circ)));
		yl=((1-s)*ylf)+(s*ylt);

                xrf=X_SIZE+right-(s*X_SIZE/2)-span;
                xrt=(X_SIZE/2)+((sin(phi)*(right+circ)));

		xr=((1-s)*xrf)+(s*xrt);

                yrf=i*Y_SIZE*2/w_len;
                yrt=(Y_SIZE/2)-((cos(phi)*(right+circ)));

		yr=((1-s)*yrf)+(s*yrt);


                for (sx=0;sx<size;sx++)
                {
                        for (sy=0;sy<size;sy++)
                        {
                                if (totl) {plotf( image,xl+sx,yl+sy,255,255,255);}
                                if (totr) {plotf( image,xr+sx,yr+sy,255,255,255);}
                                //plotf( image,xl+sx,yl+sy,255,255,255);
                                //plotf( image,xr+sx,yr+sy,255,255,255);
                        }
                }
        }
}


int add_wav_lin (short *wav, float *image, float xp, float yp, int slice,float max)
{
	int sx,sy,w_len,wp,i,xl,yl,xr,yr,size,totl,totr,l,r,thresh;
	float s,rad,span,xp1,yp1;

	s=(float)slice/18000;
	rad=1500;span=250;
	w_len=44100*2/30;
	thresh=w_len*100;
	wp=slice*w_len;

	size=1; totl=0; totr=0;
	for (i=0;i<w_len/2;i++)
	{
		l=wav[wp+(2*i)];
		r=wav[wp+(2*i)+1];
		if (l>0){ totl+=l;}else{totl-=l;}
		if (r>0){ totr+=r;}else{totr-=r;}
	}
	if (totl<thresh){totl=0;}
	if (totr<thresh){totr=0;}
	
	for (i=0;i<w_len/2;i++)
	{
		float left,right;
		left=wav[wp+(2*i)]*span/28768;
		right=wav[wp+(2*i)+1]*span/32768;

		xl=left+(s*X_SIZE/2)+span;
		yl=i*Y_SIZE*2/w_len;

		xr=X_SIZE+right-(s*X_SIZE/2)-span;
		yr=i*Y_SIZE*2/w_len;
		
		for (sx=0;sx<size;sx++)
		{
			for (sy=0;sy<size;sy++)
			{
				if (totl) {plotf( image,xl+sx,yl+sy,255,255,255);}
				if (totr) {plotf( image,xr+sx,yr+sy,255,255,255);}
			}
		}
	}
}

int add_wav_rad (short *wav, float *image, float xp, float yp, int slice,float max)
{
	int sx,sy,w_len,wp,i,xl,yl,xr,yr,size;
	float rad,span,xp1,yp1;
	float ang,left,right,off,s;

	s=(float)slice/18000;

	off=30*M_PI*s;

	xp1=(X_SIZE)-xp;
	yp1=(Y_SIZE)-yp;
	//xp1=xp;
	//yp1=yp;
	//
	xp=(X_SIZE/2)+((1-s)*X_SIZE*sin(2*M_PI*12*s)/2);
	yp=(Y_SIZE/2)+((1-s)*Y_SIZE*cos(2*M_PI*12*s)/2);

	xp1=(X_SIZE/2)-((1-s)*(1-s)*X_SIZE*sin(2*M_PI*14*s)/2);
	yp1=(Y_SIZE/2)-((1-s)*(1-s)*Y_SIZE*cos(2*M_PI*14*s)/2);

	rad=30;span=120;

	w_len=44100*2/30;
	wp=slice*w_len;
	size=1;

	int leftc,rightc,leftp,rightp;
	int left_start,right_start;

	leftc=0;rightc=0;
	leftp=0;rightp=0;

	for (i=0;i<w_len/2;i++)
	{
		left=wav[wp+(2*i)];
		right=wav[wp+(2*i)+1];
		if (leftc==0 && left<0){ leftc=1;}
		if (rightc==0 && right<0){ rightc=1;}
		if (leftc==1 && left>0){ leftc=2; leftp=i;}
		if (rightc==1 && right>0){ rightc=2; rightp=i;}
	}


	for (i=0;i<w_len/2;i++)
	{
		int lefto,righto;
		float phi;
		lefto=i+leftp;if (lefto>=(w_len/2)){lefto-=(w_len/2);}
		righto=i+rightp; if (righto>=(w_len/2)){righto-=(w_len/2);}
		left=wav[wp+(2*lefto)]*span/32768;
		right=wav[wp+(2*righto)+1]*span/32768;
		ang=(float)i*4*M_PI/(float)w_len;

		xl=xp+((rad+left)*sin(ang+off));
		yl=yp+((rad+left)*cos(ang+off));

		xr=xp1+(((rad)+right)*sin(ang+off));
		yr=yp1+(((rad)+right)*cos(ang+off));
		float rr,gg,bb,r,g,b;
		rr=max*(1-s);
		gg=max*s;
		bb=max;

		/*r=2*(((float)i*bb)+((((float)w_len/2)-(float)i)*rr))/(float)w_len;
		g=2*(((float)i*gg)+((((float)w_len/2)-(float)i)*bb))/(float)w_len;
		b=2*(((float)i*rr)+((((float)w_len/2)-(float)i)*gg))/(float)w_len;*/

		phi=8*(float)i*M_PI/(float)w_len;

		r=(1+sin(phi))*rr/2;
		g=(1+sin(phi+(2*M_PI/3)))*gg/2;
		b=(1+sin(phi+(4*M_PI/3)))*bb/2;


		for (sx=0;sx<size;sx++)
		{
			for (sy=0;sy<size;sy++)
			{
				plotf( image,xl+sx,yl+sy,r,g,b);
				plotf( image,xr+sx,yr+sy,g,b,r);
			}
		}
	}	
}


int main(int argc,char *argv[])
{
	float *image1,*image2;
	unsigned char *image3;
	short *wav,*out,*p1,*p2;
	char *mask[3];
	char name[200];
	int dims[3],x_size,y_size,off,x,y,i;
	float    max,min;
	float left[18000];
	float right[18000];
	long length,lp1,lp2;
  	DIR *dp;
  	struct dirent *ep;     
	init_x();

        /*image1=(float *)malloc(sizeof (float)*X_SIZE*Y_SIZE*3); // disp buffer
        image2=(float *)malloc(sizeof (float)*X_SIZE*Y_SIZE*3); // disp buffer*/
        image1=(float *)malloc(sizeof (float)*X_SIZE*Y_SIZE*3); // disp buffer
        image2=(float *)malloc(sizeof (float)*X_SIZE*Y_SIZE*3); // disp buffer
        image3=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        mask[0]=(char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        mask[1]=(char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        mask[2]=(char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        wav=(short *)malloc(sizeof (short)*2*60*11*44100); // disp buffer
        out=(short *)malloc(sizeof (short)*2*60*11*44100); // disp buffer
        p1=(short *)malloc(sizeof (short)*2*60*11*44100); // disp buffer
        p2=(short *)malloc(sizeof (short)*2*60*11*44100); // disp buffer
							   //
							   //
	length=load_wav(wav,"./orch.wav");
	printf ("Loaded orch %ld \n",length);
	lp1=load_wav(p1,"./divided.wav");
	printf ("Loaded p1 %ld \n",lp1);
	//lp2=load_wav(p2,"./p2.wav");
	//printf ("Loaded p2 %ld \n",lp2);
	for (i=0;i<length;i++){ out[i]=0;}

	float maxl,maxr,minl,minr;
	maxl=0;maxr=0;
	minl=59070485; minr=59070485;
	int secs;
	secs=600;
	/*
	for (i=0;i<600;i++)
	{
		//44100 pers sec
		//30 frames per sec 1470 samps per frame;
		int tick,ll,rr,wp,ulp,urp;
		ll=0;rr=0;
		wp=i*(2*44100);
		for (tick=0;tick<44100;tick++)
		{
			ulp=wav[wp+(tick*2)];
			urp=wav[wp+(tick*2)+1];
			if (ulp<0){ ll-=ulp;}else{ll+=ulp;}
			if (urp<0){ rr-=urp;}else{rr+=urp;}
		}
		if (ll<minl){minl=ll;} if (rr<minr){minr=rr;}
		if (ll>maxl){maxl=ll;} if (rr>maxr){maxr=rr;}
		left[i]=ll;
		right[i]=rr;
	}
	*/
	/*
	for (i=0;i<600;i++)
	{
		left[i]=(left[i]-minl)/(maxl-minl);
		right[i]=(right[i]-minr)/(maxr-minr);
		printf("i %d l %f r %f maxl %f minl %f\n",i,left[i],right[i], maxl, minl);
	}

	int current_chunk,wanted_gap,along,chunk_size,per,pp1,pp2;
	chunk_size=10;
	per=10;
	along=0;
	pp1=0;pp2=0;
	for (i=0;i<length/2;i++)
	{
		out[i*2]=wav[i*2]+((1-left[i/44100])*p1[pp1]);
		out[(i*2)+1]=wav[(i*2)+1]+((1-right[i/44100])*p2[pp2]);
		pp1+=2;pp2+=2;
		if (pp1>lp1){pp1=0;}
		if (pp2>lp2){pp2=0;}
	}
	save_wav(out,"./out.wav", 44100, 2, length );*/

	maxl=0;maxr=0;
	minl=5970485; minr=5970485;
	for (i=0;i<18000;i++)
	{
		//44100 pers sec
		//30 frames per sec 1470 samps per frame;
		int tick,ll,rr,wp,ulp,urp;
		ll=0;rr=0;
		wp=i*(2*1460);
		for (tick=0;tick<1460;tick++)
		{
			ulp=wav[wp+(tick*2)];
			urp=wav[wp+(tick*2)+1];
			if (ulp<0){ ll-=ulp;}else{ll+=ulp;}
			if (urp<0){ rr-=urp;}else{rr+=urp;}
		}
		if (ll<minl){minl=ll;} if (rr<minr){minr=rr;}
		if (ll>maxl){maxl=ll;} if (rr>maxr){maxr=rr;}
		left[i]=ll;
		right[i]=rr;
	}
	for (i=0;i<18000;i++)
	{
		left[i]=(left[i]-minl)/(maxl-minl);
		right[i]=(right[i]-minr)/(maxr-minr);
		printf("i %d l %f r %f maxl %f minl %f\n",i,left[i],right[i], maxl, minl);
	}


	for (y=0;y<Y_SIZE;y++)
	{
		float r,g,b,xx,yy;
		int p;
		p=y*3*X_SIZE;
		for(x=0;x<X_SIZE;x++)
		{
			r=165+(90*(sin((float)x*(float)(2)*M_PI/(X_SIZE*3))));
			g=165+(90*(sin((float)y*(float)(2)*M_PI/(Y_SIZE*3))));
			b=165+(90*(cos((float)(y+x)*(float)(2)*M_PI/(Y_SIZE*4))));
			r=0;g=0;b=0;
			image1[p+(x*3)]=r; image1[p+(x*3)+1]=g; image1[p+(x*3)+2]=b;
			if (rand()%500==1){image1[p+(x*3)]=128; image1[p+(x*3)+1]=128; image1[p+(x*3)+2]=255;}
		       	mask[0][p+(x*3)]=1; mask[1][p+(x*3)+1]=1; mask[2][p+(x*3)+2]=1;
		}

	}

	float slice,size,size2,pixel,pixel2;
	float cvl,cvr,cvtot,cvplus;
	int col,id[3],t,ud,fo;
	col=0;
	cvl=0;cvr=0;cvtot=0;
	ud=1;
	fo=1;
	for (slice=0;slice<18000;slice+=1)
	{
		float theta,m,c,sheer,s,thetaa,m2,c2,q;
		float x1,x2,y1,y2,xx1,xx2,yy1,yy2;

		m=1;
		c=Y_SIZE;
		fo=1-fo;

		s=(float)slice/18000;
		q=(1-s);
		max=40+(215*s*s);
		//min=40-(19*((1-cos(5*M_PI*s))));
		min=40-(s*40);

		pixel=(cos(s*2*M_PI*20));
		pixel2=(sin(s*2*M_PI*19));

		cvl=left[(int)slice];
		cvr=right[(int)slice];
		cvtot+=(cvl+cvr);
		cvl=2*sqrt(cvl);cvr=2*sqrt(cvr);
		cvplus=(cvl+cvr);

		if (cvtot>800)
		{
			cvtot=0;
			ud=1-ud;
			for (y=0;y<Y_SIZE;y++)
			{
				float r,g,b,xx,yy;
				int p;
				p=y*3*X_SIZE;
				for(x=0;x<X_SIZE;x++)
				{
					//r=165+(90*(sin((double)x*(float)(2)*M_PI/(slice))));
					//g=165+(90*(sin((double)y*(float)(2)*M_PI/(slice))));
					//b=165+(90*(cos((double)(y+x)*(float)(2)*M_PI/(slice))));
					//image1[p+(x*3)]+=r; image1[p+(x*3)+1]+=g; image1[p+(x*3)+2]+=b;
					//image1[p+(x*3)]/=2; image1[p+(x*3)+1]/=2; image1[p+(x*3)+2]/=2;
		       			//mask[0][p+(x*3)]=1; mask[1][p+(x*3)+1]=1; mask[2][p+(x*3)+2]=1;
				}
			}
		}	

		theta=M_PI-(cvl*M_PI/2);
		thetaa=(cvr*M_PI/2);


		//xoffc=(Y_SIZE/2)*(cos(-thetaa));
		//yoffc=(Y_SIZE/2)*(sin(-thetaa));

		//xxoffc=s*(Y_SIZE/2)*(cos(-theta));
		//yyoffc=s*(Y_SIZE/2)*(sin(-theta));

		//xoffc=(X_SIZE/4); yoffc=s*100;
		//xxoffc=(-X_SIZE/4); yyoffc=-s*100;

		
		x1=(X_SIZE/2)+((100)*cos(-thetaa))-(s*X_SIZE);
		x2=(X_SIZE/2)+((100)*cos(-thetaa+(M_PI/6)))-(s*X_SIZE);
		y1=(Y_SIZE-(Y_SIZE*s))+(100*sin(-thetaa));
		y2=(Y_SIZE-(Y_SIZE*s))+(100*sin(-thetaa+(M_PI/6)));

		xx1=(X_SIZE/2)+((100)*cos(-theta))+(s*X_SIZE);
		xx2=(X_SIZE/2)+((100)*cos(-theta-(M_PI/6)))+(s*X_SIZE);
		yy1=(Y_SIZE-(Y_SIZE*s))+(100*sin(-theta));
		yy2=(Y_SIZE-(Y_SIZE*s))+(100*sin(-theta-(M_PI/6)));

		t=60+(s*120);
		printf("doing %f %f %f %f\n",slice,cvtot,max,min);
	
		if (x2==x1){ x1+=0.001;}
		if (xx2==xx1){ xx1+=0.001;}
		m=(y2-y1)/(x2-x1);
		m2=(yy2-yy1)/(xx2-xx1);
		//y=mx+c
		c=y1-(m*x1);
		c2=yy1-(m2*xx1);

	
	for (y=0;y<Y_SIZE;y++)
		{
			int p;
			float ypect,vib,wib,r,yp;
			yp=(cvl*40*sin((float)y*2*M_PI*pixel2*cvr/500));
			p=y*3*X_SIZE;
			for(x=0;x<X_SIZE;x++)
			{
				int i,point;
				if (fo) {ypect=m*(float)x+c;}else{
				ypect=m2*(float)x+c2;}

				ypect+=yp+(cvr*40*sin((float)x*2*M_PI*pixel*cvl/500));

				if ((ypect<y && ud ) || (ypect>y && !ud))
				{
					for (i=0;i<3;i++)
					{
						point=p+(x*3)+i;
						if (i==(col/t))
						{
							if (image1[point]+mask[i][point]>max){ mask[i][point]=0;}
							else if (image1[point]+mask[i][point]<min){ mask[i][point]=1;}
							if (mask[i][point])
							{
							image2[point]=image1[point]+cvplus;
							}else{
							image2[point]=image1[point]-cvplus;
							}
						} else {
				 			image2[point]=image1[point];
						}

					}
				} else{
					for (i=0;i<3;i++)
					{
						point=p+(x*3)+i;
						image2[point]=image1[point];
					}
				}
				if (rand()%5000==1){image2[p+(x*3)]=max; image2[p+(x*3)+1]=max; image2[p+(x*3)+2]=max;}
			}
		}


		col++;if (col>(3*t)){col=0;}

		float zoomx,zoomy;

		zoomx=(X_SIZE/2)+(((0.5-s)*Y_SIZE/2)*sin(2*M_PI*s*15));
	     	zoomy=(Y_SIZE/2)+(((0.5-s)*Y_SIZE/2)*cos(2*M_PI*s*15));

		add_wav_rad (wav,image2,zoomx,zoomy,slice,max);
		if (slice>=2700){add_wav_cur (p1,image2,zoomx,zoomy,slice-2700,max);}
		zoom_it (image2,image1, 0.9993-(0.0005*cvtot/300), zoomx, zoomy);

		for (i=0;i<3*X_SIZE*Y_SIZE;i++)
		{
			image3[i]=image2[i];
			//image1[i]=image2[i];
		}
		//memcpy(image1,image2,X_SIZE*Y_SIZE*3);
		//blurt(image2,2);
		disp (image3,slice,1);
	}
	disp (image3,18009,1);
	scanf("%c",name);	
	exit (0);


	close_x();

	exit(0);
}	

void disp (unsigned char *image2,int fram,int ab)
{
	int x,y;
	char *input;
	input=malloc(300);


       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*X_SIZE*3;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=image2[xpoint];
			x_buffer[X_POINT+1]=image2[xpoint+1];
			x_buffer[X_POINT]=image2[xpoint+2];
                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/kinal%05d.jpg",fram);
	if (ab){jayit(image2,X_SIZE, Y_SIZE, input);}
	free (input);
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};

