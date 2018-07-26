#ifndef CGLS_H 
#define CGLS_H 
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
/*Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is
a small number that protects against trying to achieve fractional accuracy for a minimum that
happens to be exactly zero. 
Here ITMAX is the maximum allowed number of iterations, while EPS is a small number to
rectify the special case of converging to exactly zero function value.
 Here GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the
maximum magnification allowed for a parabolic-fit step. */

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Matrix<dcomp, 9, 9> m3;

template <typename func1>
void dxargs(VectorXd&, VectorXd&, func1&&, const vector<VectorXd>&,
	       	const m3&, const vector<Vector3d>&,
		const VectorXd&, const int, const int);

template <typename func1, typename... Args>
double f1dim(double x, int n, VectorXd &pcom, VectorXd &xicom, func1&& func, Args&&... params)
{
	double f;
	VectorXd xt(n);
	xt = pcom + x*xicom;
	f=forward<func1>(func)(xt, forward<Args>(params)...);
	return f;
}

/* template <typename func1, typename... Args> */
/* void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc, func1&& func, Args&&... params) */
template <typename... Args>
void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc, Args&&... params)
/* Given a function func, and given distinct initial points ax and bx, this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa, fb, and fc. */
{
	double GOLD = 1.618034;
	double TINY = 1e-20;
	double GLIMIT = 100.;
	double ulim,u,r,q,fu,dum;
	double tmp1, tmp2;
	fa=f1dim(ax, params...);
	fb=f1dim(bx, params...);
	if (fb > fa) {		// Switch roles of a and b so that we can go
		dum = ax;
		ax = bx;
		bx = dum;	// downhill in the direction from a to b.
		dum = fb;
		fb = fa;
		fa = dum;
	}
	cx=(bx)+GOLD*(bx-ax);	// First guess for c.
	fc=f1dim(cx, params...);
	while (fb > fc) {		// Keep returning here until we bracket.
		r=(bx-ax)*(fb-fc);	// Compute u by parabolic extrapolation from
		q=(bx-cx)*(fb-fa);	// a, b, c. TINY is used to prevent any possible division by zero.
		tmp1 = abs(q-r);
		tmp2 = ((tmp1 > TINY) ? tmp1 : TINY);
		tmp1 = ((q-r) >= 0.0 ? abs(tmp2) : -abs(tmp2));
		u=(bx)-((bx-cx)*q-(bx-ax)*r)/(2.0*tmp1);
		ulim=(bx)+GLIMIT*(cx-bx);
		// We won’t go farther than this. Test various possibilities:
		if ((bx-u)*(u-cx) > 0.0) {	// Parabolic u is between b and c: try it.
			fu=f1dim(u, params...);
			if (fu < fc) {		// Got a minimum between b and c.
				ax=bx;
				bx=u;
				fa=fb;
				fb=fu;
				return;
			}
		       	else if (fu > fb) {	// Got a minimum between between a and u.
				cx=u;
				fc=fu;
				return;
			}
			u=cx+GOLD*(cx-bx);	// Parabolic fit was no use. Use default magnification
			fu=f1dim(u, params...);
		} else if ((cx-u)*(u-ulim) > 0.0) {	// Parabolic fit is between c and its allowed limit.
			fu=f1dim(u, params...);
			if (fu < fc) {
				bx = cx;
				cx = u;
				u = cx + GOLD*(cx-bx);
				fb = fc;
				fc = fu;
				fu=f1dim(u, params...);
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) {	// Limit parabolic u to maximum
			u=ulim;					// allowed value.
			fu=f1dim(u, params...);
		} else {	// Reject parabolic u, use default magnification
			u=cx+GOLD*(cx-bx);
			fu=f1dim(u, params...);
		}
		ax = bx;	// Eliminate oldest point and continue.
		bx = cx;
		cx = u;
		fa = fb;
		fb = fc;
		fc = fu;
	}
}

/* template <typename func1, typename... Args> */
/* double brent(double ax, double bx, double cx, double tol, double &xmin, func1&& func, Args&&... params) */
template <typename... Args>
double brent(double ax, double bx, double cx, double tol, double &xmin, Args&&... params)
/*Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
between ax and cx, and f(bx) is less than both f(ax) and f(cx) ), this routine isolates
the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
the minimum is returned as xmin , and the minimum function value is returned as brent , the
returned function value. */
{
 	double CGOLD = 0.3819660;
	double ZEPS = 1e-10;
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;	// This will be the distance moved on the step before last.
	a=(ax < cx ? ax : cx); 	// a and b must be in ascending order,
	b=(ax > cx ? ax : cx);	// but input abscissas need not be.
	x=w=v=bx; 		// Initializations...
	fw=fv=fx=f1dim(x, params...);
	int ITMAX = 100;
	for (iter=1;iter<=ITMAX;iter++) {	// Main program loop.
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*abs(x)+ZEPS);
		if (abs(x-xm) <= (tol2-0.5*(b-a))) {	// Test for done here.
			xmin=x;
			return fx;
		}
		if (abs(e) > tol1) {	// Construct a trial parabolic fit.
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=abs(q);
			etemp=e;
			e=d;
			if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
/* The above conditions determine the acceptability of the parabolic fit. Here we
take the golden section step into the larger of the two segments. */
			else {
				d=p/q;	// Take the parabolic step.
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d = ((xm - x) >= 0.0 ? abs(tol1) : -abs(tol1));
			}
		} 
		else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(abs(d) >= tol1 ? x+d : x+(d >= 0.0 ? abs(tol1) : -abs(tol1)));
		fu = f1dim(u, params...);
		// This is the one function evaluation per iteration.
		if (fu <= fx) {	// Now decide what to do with our function evaluation.
			if (u >= x) 
				a=x; 
			else 
				b=x;
			// Housekeeping follows:
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
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
	cout<<"Too many iterations in brent"<<endl;
	xmin=x;		// Never get here.
	return fx;
}

template <typename... Args>
void linmin(VectorXd &p, VectorXd &xi, int n, double &fret, Args&&... params)
/*Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p. This is actually all accomplished by calling the
routines mnbrak and brent.*/
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	VectorXd pcom(n), xicom(n);
	pcom = p;
	xicom = xi;
	ax = 0.0;		//Initial guess for brackets.
	xx = 1.0;
	double TOL = 2e-4;
	mnbrak(ax,xx,bx,fa,fx,fb,n,pcom,xicom,params...);
	fret = brent(ax,xx,bx,TOL,xmin,n,pcom,xicom,params...);
	xi = xi*xmin;
	p = p + xi;
}

template <typename func1, typename... Args>
void frprmn(VectorXd &p, double ftol, int &iter, double &fret,
	func1&& func, Args&&... params)
{
/* Given a starting point p[1..n] , Fletcher-Reeves-Polak-Ribiere minimization is performed on a
function func , using its gradient as calculated by a routine dfunc . The convergence tolerance
on the function value is input as ftol . Returned quantities are p (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine linmin is called to perform line minimizations. */

	int j,its;
	int n = p.size();
	double gg,gam,fp,dgg;
	VectorXd g(n), h(n), xi(n);
	fp = forward<func1>(func)(p, forward<Args>(params)...);
	dxargs(p, xi, func, params...);
	g = -xi;
	h = g;
	xi = h;
	double EPS = 1e-10;
	int ITMAX = 500;
	for (its=1;its<=ITMAX;its++) {	//Loop over iterations.
		iter=its;
		linmin(p,xi,n,fret,func,params...);	//Next statement is the normal return:
		if (2.0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) 
			return;
		fp = fret;
		dxargs(p, xi, func, params...);
		dgg = gg = 0.0;
		gg = g.dot(g);
		/* dgg = xi.dot(g); 		//This statement for Fletcher-Reeves. */
		dgg = (xi + g).dot(xi);		//This statement for Polak-Ribiere.
		if (gg == 0.0)  		// Unlikely. If gradient is exactly zero then
			return;			//we are already done.
		gam = dgg/gg;
		g = -xi;
		h = g + gam*h;
		xi = h;
	}
	cout<<"Too many iterations in frprmn"<<endl;
}
#endif
