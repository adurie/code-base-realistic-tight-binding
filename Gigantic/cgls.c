#include <math.h>
/* #include "nrutil.h" */
#define ITMAX 200
#define EPS 1.0e-10
#define TOL 2.0e-4	//Tolerance passed to brent.
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
/*Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is
a small number that protects against trying to achieve fractional accuracy for a minimum that
happens to be exactly zero. */
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
/*
Here ITMAX is the maximum allowed number of iterations, while EPS is a small number to
rectify the special case of converging to exactly zero function value.
i*/
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
/* Here GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the
maximum magnification allowed for a parabolic-fit step. */

int ncom;		//Global variables communicate with f1dim.
float *pcom,*xicom,(*nrfunc)(float []);
void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []))
/*Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p. This is actually all accomplished by calling the
routines mnbrak and brent.*/
{
	float brent(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin);
	float f1dim(float x);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float (*func)(float));
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;		//Define the global variables.
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;		//Initial guess for brackets.
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {i	//Construct the vector results to return.
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}

float f1dim(float x)
{
	int j;
	float f,*xt;
	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
		f=(*nrfunc)(xt);
		free_vector(xt,1,ncom);
	return f;
}

float brent(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin)
/*Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
between ax and cx, and f(bx) is less than both f(ax) and f(cx) ), this routine isolates
the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
the minimum is returned as xmin , and the minimum function value is returned as brent , the
returned function value. */
{
	int iter;
	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;	// This will be the distance moved on the step before last.
	a=(ax < cx ? ax : cx); 	// a and b must be in ascending order,
	b=(ax > cx ? ax : cx);	// but input abscissas need not be.
	x=w=v=bx; 		// Initializations...
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {	// Main program loop.
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {	// Test for done here.
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {	// Construct a trial parabolic fit.
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
/* The above conditions determine the acceptability of the parabolic fit. Here we
take the golden section step into the larger of the two segments. */
			else {
				d=p/q;	// Take the parabolic step.
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} 
		else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		// This is the one function evaluation per iteration.
		if (fu <= fx) {	// Now decide what to do with our function evaluation.
			if (u >= x) 
				a=x; 
			else 
				b=x;
			SHFT(v,w,x,u)	// Housekeeping follows:
			SHFT(fv,fw,fx,fu)
		}
	       	else {
			if (u < x) 
				a=u;
		       	else 
				b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}		//Done with housekeeping. Back for another iteration.
	}
	nrerror("Too many iterations in brent");
	*xmin=x;		// Never get here.
	return fx;
}

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float (*func)(float))
/* Given a function func, and given distinct initial points ax and bx, this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa, fb, and fc. */
{
	float ulim,u,r,q,fu,dum;
	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {		// Switch roles of a and b so that we can go
		SHFT(dum,*ax,*bx,dum)	// downhill in the direction from a to b.
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);	// First guess for c.
	*fc=(*func)(*cx);
	while (*fb > *fc) {		// Keep returning here until we bracket.
		r=(*bx-*ax)*(*fb-*fc);	// Compute u by parabolic extrapolation from
		q=(*bx-*cx)*(*fb-*fa);	// a, b, c. TINY is used to prevent any possible division by zero.
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		// We won’t go farther than this. Test various possibilities:
		if ((*bx-u)*(u-*cx) > 0.0) {	// Parabolic u is between b and c: try it.
			fu=(*func)(u);
			if (fu < *fc) {		// Got a minimum between b and c.
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			}
		       	else if (fu > *fb) {	// Got a minimum between between a and u.
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);	// Parabolic fit was no use. Use default mag-
			fu=(*func)(u);		// nification.
		} else if ((*cx-u)*(u-ulim) > 0.0) {	// Parabolic fit is between c and its
			fu=(*func)(u);			// allowed limit.
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {	// Limit parabolic u to maximum
			u=ulim;					// allowed value.
			fu=(*func)(u);
		} else {	// Reject parabolic u, use default magnification
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)	// Eliminate oldest point and continue.
		SHFT(*fa,*fb,*fc,fu)
	}
}

void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
{
/* Given a starting point p[1..n] , Fletcher-Reeves-Polak-Ribiere minimization is performed on a
function func , using its gradient as calculated by a routine dfunc . The convergence tolerance
on the function value is input as ftol . Returned quantities are p (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine linmin is called to perform line minimizations. */

	void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []));
	int j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;
	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);			//Initializations.
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {	//Loop over iterations.
		*iter=its;
		linmin(p,xi,n,fret,func);	//Next statement is the normal return:
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp= *fret;
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			/* dgg += xi[j]*xi[j]; */ 	//This statement for Fletcher-Reeves.
			dgg += (xi[j]+g[j])*xi[j]; 	//This statement for Polak-Ribiere.
		}
		if (gg == 0.0) { 		// Unlikely. If gradient is exactly zero then
			FREEALL				//we are already done.
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
