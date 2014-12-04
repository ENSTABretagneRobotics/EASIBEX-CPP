// Compatiblity layer with ibex for the simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#include "interval.h"

#if defined(_MSC_VER) || defined(__BORLANDC__) 
// Enable the use of isnan().
#include <float.h>
#ifndef isnan
#define isnan _isnan
#endif // isnan
// Used to define NAN (Not A Number).
#ifdef NAN_CONSTS_NEEDED
const unsigned long nan[2] = {0xffffffff, 0xffffffff};
const double nan_double = -*(double*)nan;
#endif // NAN_CONSTS_NEEDED
#endif // defined(_MSC_VER) || defined(__BORLANDC__) 

#ifdef NAI_CONST_NEEDED
// Used to define NAI (Not An Interval).
const interval nai = interval();
#endif // NAI_CONST_NEEDED

#include <cmath>

using namespace std;
using namespace ibex;

//----------------------------------------------------------------------
// Useful real-valued functions
//----------------------------------------------------------------------
double Min(vector<double>& x)
{
	double d = oo;
	for (unsigned int i = 0; i < x.size(); i++)
		d = min(x[i], d);
	return d;
}
//----------------------------------------------------------------------
double Max(vector<double>& x)
{
	double d = -oo;
	for (unsigned int i = 0; i < x.size(); i++)
		d = max(x[i], d);
	return d;
}
//----------------------------------------------------------------------
double Sign(const double x)
{
	if (x > 0) return (1);
	else return (0);
}
//----------------------------------------------------------------------
double Chi(const double a, const double b, const double c)
{
	if (a < 0) return (b);
	else return (c);
}
//----------------------------------------------------------------------
double Arccossin(const double x, const double y)
{
	if (y > 0) return (acos(x));
	else return (-acos(x));
}
//----------------------------------------------------------------------
double Arg(const double x, const double y)
{
	double r = sqrt(x*x + y*y);
	if (r == 0) return 0;
	return Arccossin(x / r, y / r);
}
//----------------------------------------------------------------------
double Det(double ux, double uy, double vx, double vy)
{
	return (ux*vy - vx*uy);
}
//----------------------------------------------------------------------
double DistanceDirSegment(double mx, double my, double theta, double ax, double ay, double bx, double by)
{      
	// Distance directionnelle du point m au segment [a,b]. La direction est donnee par theta
	double ma_x = ax - mx;
	double ma_y = ay - my;
	double mb_x = bx - mx;
	double mb_y = by - my;
	double ab_x = bx - ax;
	double ab_y = by - ay;
	double ux = cos(theta);
	double uy = sin(theta);
	double z1 = Det(ma_x, ma_y, ux, uy);
	double z2 = Det(ux, uy, mb_x, mb_y);
	double z3 = Det(ma_x, ma_y, ab_x, ab_y);
	double z4 = Det(ux, uy, ab_x, ab_y);
	double z5 = min(z1, min(z2, z3));
	double d1 = 0;
	if (z4 == 0) d1 = oo;
	else d1 = z3 / z4;
	return (Chi(z5, oo, d1));
}
//----------------------------------------------------------------------
void DistanceDirSegment(double& d,double& phi, double mx, double my, double theta, double ax, double ay, double bx, double by)
{      
	// Distance directionnelle du point m au segment [a,b].
	double ma_x=ax-mx;
	double ma_y=ay-my;
	double mb_x=bx-mx;
	double mb_y=by-my;
	double ab_x=bx-ax;
	double ab_y=by-ay;
	double ux=cos(theta);
	double uy=sin(theta);
	double z1=Det(ma_x,ma_y,ux,uy);
	double z2=Det(ux,uy,mb_x,mb_y);
	double z3=Det(ma_x,ma_y,ab_x,ab_y);
	double z4=Det(ux,uy,ab_x,ab_y);
	double z5=min(z1,min(z2,z3));
	double d1=z3/z4;
	d=Chi(z5,oo,d1);
	phi=atan2(-ab_x,ab_y); //phi is the angle of the normal vector of [a,b]
}
//----------------------------------------------------------------------
double DistanceDirSegments(double mx, double my, double theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// Distance directionnelle relativement a un polygone
	vector<double> dist(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++)
		dist[j] = DistanceDirSegment(mx, my, theta, ax[j], ay[j], bx[j], by[j]);
	double distmin = Min(dist);
	return (distmin);
}
//----------------------------------------------------------------------
void DistanceDirSegments(double& d, double& phi, double mx, double my, double theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{     
	// Distance directionnelle relativement a un polygone (le triedre m-a-b doit etre direct, sinon cette distance est infinie)
	d = oo;
	for (unsigned int j = 0; j < ax.size(); j++)
	{
		double dj,phij;
		DistanceDirSegment(dj,phij,mx,my,theta,ax[j],ay[j],bx[j],by[j]);
		if (dj < d) { d=dj;phi=phij; };
	}
}
//----------------------------------------------------------------------
double DistanceDirCircle(double mx, double my, double theta, double cx, double cy, double r)
{      
	// Distance directionnelle du point m au cercle de centre c et de rayon r.
	double ux, uy, alpha, beta, a, b, c, delta, px1, py1, px2, py2, d1, d2;
	ux = cos(theta);
	uy = sin(theta);
	if (fabs(uy) > 0.00)  //pour eviter la division par zero. Il conviendrait de traiter le cas uy=0
	{
		alpha = ux / uy;
		beta = mx - my*alpha;
		a = alpha*alpha + 1;
		b = 2 * alpha*(beta - cx) - 2 * cy;
		c = (beta - cx)*(beta - cx) + cy*cy - r*r;
		delta = b*b - 4 * a*c;
		if (delta < 0) return(oo);
		py1 = (-b - sqrt(delta)) / (2 * a);
		px1 = alpha*py1 + beta;
		py2 = (-b + sqrt(delta)) / (2 * a);
		px2 = alpha*py2 + beta;
		d1 = Chi((px1 - mx)*ux + (py1 - my)*uy, oo, sqrt((px1 - mx)*(px1 - mx) + (py1 - my)*(py1 - my)));
		d2 = Chi((px2 - mx)*ux + (py2 - my)*uy, oo, sqrt((px2 - mx)*(px2 - mx) + (py2 - my)*(py2 - my)));
		return min(d1, d2);
	}
	return oo;
}
//----------------------------------------------------------------------
void DistanceDirCircle(double& d, double& phi, double mx, double my, double theta, double cx, double cy, double r)
{      
	//  d is the directional distance from m to the circle of center c and radius r.
	//  phi is the angle of the normal vector of of the impact point
	double ux,uy,alpha,beta,a,b,c,delta,px1,py1,px2,py2,d1,d2,phi1,phi2;
	ux=cos(theta);
	uy=sin(theta);
	if (fabs(uy)<0.01)    // pour eviter une division par zero, on permutte x et y
	{ double aux; aux=mx; mx=my; my=aux;  aux=cx; cx=cy; cy=aux;    aux=ux; ux=uy; uy=aux; }
	alpha=ux/uy;
	beta=mx-my*alpha;
	a=alpha*alpha+1;
	b=2*alpha*(beta-cx)-2*cy;
	c=(beta-cx)*(beta-cx)+cy*cy-r*r;
	delta=b*b-4*a*c;
	if (delta<0) {d=oo;phi=-999;return;};
	py1=(-b-sqrt(delta))/(2*a);
	px1=alpha*py1+beta;
	phi1=atan2(cy-py1,cx-px1);
	py2=(-b+sqrt(delta))/(2*a);
	px2=alpha*py2+beta;
	phi2=atan2(cy-py2,cx-px2);

	d1=oo;d2=oo;
	if ((px1-mx)*ux+(py1-my)*uy>0)   d1=sqrt((px1-mx)*(px1-mx)+(py1-my)*(py1-my));
	if ((px2-mx)*ux+(py2-my)*uy>0)   d2=sqrt((px2-mx)*(px2-mx)+(py2-my)*(py2-my));
	if (d1<d2) { d=d1; phi=phi1; } else { d=d2; phi=phi2; }
}
//----------------------------------------------------------------------
double DistanceDirCircles(double mx, double my, double theta, vector<double> cx, vector<double> cy, vector<double> r)
{      
	// Distance directionnelle relativement a plusieurs cercles.
	vector<double> dist(cx.size());
	for (unsigned int j = 0; j < cx.size(); j++)
		dist[j] = DistanceDirCercle(mx, my, theta, cx[j], cy[j], r[j]);
	double distmin = Min(dist);
	return (distmin);
}
//----------------------------------------------------------------------
void DistanceDirCircles(double& d,double& phi, double mx, double my, double theta, vector<double> cx, vector<double> cy, vector<double> r)
{       
	d = oo;
	for (unsigned int j = 0; j < cx.size(); j++)
	{
		double dj, phij;
		DistanceDirCircle(dj,phij,mx,my,theta,cx[j],cy[j],r[j]);
		if (dj < d) { d=dj; phi=phij; };
	}
}
//----------------------------------------------------------------------
double DistanceDirSegmentsOrCircles(double mx, double my, double theta,
									vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
									vector<double> cx, vector<double> cy, vector<double> r)
{     
	double d1a, d1b;
	d1a = DistanceDirSegments(mx, my, theta, ax, ay, bx, by);
	d1b = DistanceDirCircles(mx, my, theta, cx, cy, r);
	return min(d1a, d1b);
}
//----------------------------------------------------------------------
void DistanceDirSegmentsOrCircles(double& d, double& phi, double mx, double my, double theta,
								  vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
								  vector<double> cx, vector<double> cy, vector<double> r)
{     
	// returns the distance and orientation collected by a laser rangefinder in a room made with segments and circles
	double phi1a, phi1b, d1a, d1b;
	DistanceDirSegments(d1a,phi1a,mx,my,theta,ax,ay,bx,by);
	DistanceDirCircles(d1b,phi1b,mx,my,theta,cx,cy,r);
	if (d1a < d1b) { d = d1a; phi = phi1a-theta; } else { d = d1b; phi = phi1b-theta; }
}
//----------------------------------------------------------------------
// Operators
//----------------------------------------------------------------------
//std::ostream& operator<<(std::ostream& os, const interval& a)
//{
//	if (a.is_empty()) os << "EmptyInterval";
//	else if (a.lb() != a.ub())
//	{ 
//		os << "[" << setprecision(4) << a.lb() << ", " << setprecision(4) << a.ub() << "] "; 
//	}
//	else os << a.lb();
//	return os;
//}
//----------------------------------------------------------------------
// Interval-valued functions
//----------------------------------------------------------------------
interval Min(const interval& x, const interval& y)
{
	return ibex::min(x, y);
}
//----------------------------------------------------------------------
interval Min(const interval& x, const interval& y, const interval& z)
{
	return ibex::min(ibex::min(x, y), z);
}
//----------------------------------------------------------------------
interval Max(const interval& x, const interval& y)
{
	return ibex::max(x, y);
}
//----------------------------------------------------------------------
interval Max(const interval& x, const interval& y, const interval& z)
{
	return ibex::max(ibex::max(x, y), z);
}
//----------------------------------------------------------------------
interval InvSqrt(interval& X)
{
	interval Y(-oo,oo);
	bwd_sqrt(Y, X);
	return Y;
}
//----------------------------------------------------------------------
interval PowFrac(const interval& x, int num, int den)
{
	// [x]^(num/den)
	if (x.is_empty()) return interval();
	if (num*den > 0)
	{
		double a1 = min(x.lb(), x.ub()), a2 = max(x.lb(), x.ub());
		double n = num, m = den;
		if (den % 2 == 0)
		{
			if ((a1 >= 0) || (a1*a2 <= 0))
				return interval(-pow(a2, n / m), pow(a2, n / m));
			else return interval();
		}
		else
		{
			if (a1*a2 <= 0)
				return interval(-pow(fabs(a1), n / m), pow(a2, n / m));
			if (a1 > 0) return interval(pow(a1, n / m), pow(a2, n / m));
			else return interval(-pow(fabs(a1), n / m), -pow(fabs(a2), n / m));
		}
	}
	else return interval(1.0) / PowFrac(x, abs(num), abs(den));
}
//----------------------------------------------------------------------
interval PowRoot(const interval& x, int num, int den)
{
	// [x]^(num/den)
	if (x.is_empty()) return interval();
	if (num*den > 0)
	{
		double a1 = min(x.lb(), x.ub()), a2 = max(x.lb(), x.ub()), n = num, m = den;
		if (den % 2 == 0)
		{
			if (a1 >= 0) return interval(pow(a1, n / m), pow(a2, n / m));
			if (a1*a2 <= 0) return interval(0, pow(a2, n / m));
			else return interval();
		}
		else
		{
			if (a1*a2 <= 0) return interval(-pow(fabs(a1), n / m), pow(a2, n / m));
			if (a1 > 0) return interval(pow(a1, n / m), pow(a2, n / m));
			else return interval(-pow(fabs(a1), n / m), -pow(fabs(a2), n / m));
		}
	}
	else return interval(1.0) / PowFrac(x, abs(num), abs(den));
}
//----------------------------------------------------------------------------
interval Det(interval& ux, interval& uy, interval& vx, interval& vy)
{
	return (ux*vy - vx*uy);
}
//----------------------------------------------------------------------
interval Det(interval& ux, interval& uy, double& vx, double& vy)
{
	return (vy*ux - vx*uy);
}
//----------------------------------------------------------------------
interval Inter(const interval& a, const interval& b)
{
	return a&b;
}
//----------------------------------------------------------------------
interval Inter(vector<interval> x)
{
	interval r = interval::EMPTY_SET; // r is empty
	for (unsigned int i = 0 ; i < x.size(); i++)
		r = r&x[i];
	return r;
}
//----------------------------------------------------------------------
interval Union(const interval& a, const interval& b)
{   
	return a|b;
}
//----------------------------------------------------------------------
interval Union(vector<interval> x)
{   
	interval r = interval::EMPTY_SET; // r is empty
	for (unsigned int i = 0 ; i < x.size(); i++)
		r = r|x[i];
	return r;
}
//----------------------------------------------------------------------
interval Inflate(const interval& a, double eps)
{
	interval r(a.lb() - eps, a.ub() + eps);
	if (a.is_empty()) r = interval::EMPTY_SET; // r is empty
	return interval(r);
}
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
double Inf(const interval& a)
{
	return a.lb();
}
//-----------------------------------------------------------------------
double Sup(const interval& a)
{
	return a.ub();
}
//-----------------------------------------------------------------------
double Center(const interval& a)
{
	return a.mid();
}
//-----------------------------------------------------------------------
double Width(const interval& a)
{
	return a.diam();
}
//------------------------------------------------------------------------------
double ToReal(const interval& a)
{
	if ((a.is_empty())||(a.ub() != a.lb())) return NAN;
	return a.lb();
}
//----------------------------------------------------------------------
bool Disjoint(const interval& a, const interval& b)
{
	if (a.is_empty() || b.is_empty()) return true;
	return ((a.ub() < b.lb()) || (b.ub() < a.lb()));
}
//----------------------------------------------------------------------
bool Subset(const interval& a, const interval& b)
{
	if (a.is_empty()) return true;
	if (b.is_empty()) return false;
	return ((a.lb() >= b.lb()) && (a.ub() <= b.ub()));
}
//----------------------------------------------------------------------
iboolean In(const interval& F, const interval& Y)
{
	if (Disjoint(F, Y)) return false;
	if (Subset(F, Y)) return true;
	return iboolean(iperhaps);
}
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
void Cadd(interval& Z, interval& X, interval& Y, int sens)
{
	if (sens != -1) Z &= X+Y; 
	if (sens != 1) bwd_add(Z,X,Y); 
}
//----------------------------------------------------------------------
void Cadd(interval& Z, double x, interval& Y, int sens)
{
	if (sens != -1) Z &= x+Y; 
	if (sens != 1) 
	{
		interval X(x);
		bwd_add(Z,X,Y); 
	}
}
//----------------------------------------------------------------------
void Cadd(interval& Z, interval& X, double y, int sens)
{
	if (sens != -1) Z &= X+y; 
	if (sens != 1) 
	{
		interval Y(y);
		bwd_add(Z,X,Y); 
	}
}
//----------------------------------------------------------------------
void Csub(interval& Z, interval& X, interval& Y, int sens)
{
	if (sens != -1) Z &= X-Y; 
	if (sens != 1) bwd_sub(Z,X,Y); 
}
//----------------------------------------------------------------------
void Csub(interval& Z, double x, interval& Y, int sens)
{
	if (sens != -1) Z &= x-Y; 
	if (sens != 1) 
	{
		interval X(x);
		bwd_sub(Z,X,Y); 
	}
}
//----------------------------------------------------------------------
void Csub(interval& Z, interval& X, double y, int sens)
{
	if (sens != -1) Z &= X-y; 
	if (sens != 1) 
	{
		interval Y(y);
		bwd_sub(Z,X,Y); 
	}
}
//----------------------------------------------------------------------
void Cmul(interval& Z, interval& X, interval& Y, int sens)
{
	if (sens != -1) Z &= X*Y; 
	if (sens != 1) bwd_mul(Z,X,Y); 
}
//----------------------------------------------------------------------
void Cmul(interval& Z, double x, interval& Y, int sens)
{
	if (sens != -1) Z &= x*Y; 
	if (sens != 1) 
	{
		interval X(x);
		bwd_mul(Z,X,Y); 
	}
}
//----------------------------------------------------------------------
void Cmul(interval& Z, interval& X, double y, int sens)
{
	if (sens != -1) Z &= X*y; 
	if (sens != 1) 
	{
		interval Y(y);
		bwd_mul(Z,X,Y); 
	}
}
//----------------------------------------------------------------------
void Cdiv(interval& Z, interval& X, interval& Y, int sens)
{    
	if (sens != -1) Z &= X/Y;
	if (sens != 1) bwd_div(Z,X,Y);
}
//----------------------------------------------------------------------
void Csign(interval& Y, interval& X)
{  
	Y &= sign(X);
	bwd_sign(Y,X);
}
//----------------------------------------------------------------------
void Csqr(interval &Y, interval &X, int sens)
{   
	if (sens != -1) Y &= sqr(X);
	if (sens != 1) bwd_sqr(Y,X);
}
//-----------------------------------------------------------------------
void Cexp(interval& Y, interval& X, int sens)
{   
	if (sens != -1) Y &= exp(X);
	if (sens != 1) bwd_exp(Y,X);
}
//----------------------------------------------------------------------
void Cpow(interval& Y, interval& X, int n)
{   
	static Variable x, y;
	static NumConstraint C(x,y, y=pow(x,n));
	static CtcHC4 ctc(C);
	IntervalVector P(2);  P[0]=X; P[1]=Y;
	try {ctc.contract(P);} catch(EmptyBoxException) {X=interval::EMPTY_SET;Y=X; return;}
	X=P[0]; Y=P[1];
}
//----------------------------------------------------------------------
void Ccos(interval &Y, interval &X, int sens)
{   
	if (sens != -1) Y &= cos(X);
	if (sens != 1) bwd_cos(Y,X);
}
//-----------------------------------------------------------------------
void Csin(interval &Y, interval &X, int sens)
{   
	if (sens != -1) Y &= sin(X);
	if (sens != 1) bwd_sin(Y,X);
}
//-----------------------------------------------------------------------
void Ctan(interval& Y, interval& X, int sens)
{
	if (sens == 1)  { Y &=tan(X); }
	if (sens == -1) 
	{ 
		bwd_tan(Y,X); 
		Y &= tan(X); //On est oblige de le rappeler car sinon par optimal // bug Ibex
	}
	//if (Y.is_empty()) qDebug()<<"aaa";
}
//----------------------------------------------------------------------
void Catan(interval& Y, interval& X, int sens)
{
	// Y=atan(X)                  =>  sens=1;
	// X=atan-1(Y)=tan(Y)         =>  sens=-1; 
	Ctan(X, Y, -sens);
}
//----------------------------------------------------------------------
void Cnorm(interval& N, interval& X, interval& Y)
{
	interval SqrX, SqrY, SqrN;
	SqrX = Sqr(X);
	SqrY = Sqr(Y);
	SqrN = Sqr(N);
	Cadd(SqrN, SqrX, SqrY, 1);
	Cadd(SqrN, SqrX, SqrY, -1);
	Csqr(SqrY, Y, -1);
	Csqr(SqrX, X, -1);
	Csqr(SqrN, N, -1);
}
//----------------------------------------------------------------------
void Cnorm(interval& N, interval& X, interval& Y, interval& Z, int sens)
{
	interval SqrX = Sqr(X), SqrY = Sqr(Y), SqrZ = Sqr(Z), SqrN = Sqr(N);
	if (sens != -1)
	{
		SqrN = Inter(SqrN, SqrX+SqrY+SqrZ);
		Csqr(SqrN, N, -1);
	}
	if (sens != 1)
	{
		SqrZ = Inter(SqrZ, SqrN-SqrX-SqrY);
		SqrY = Inter(SqrY, SqrN-SqrX-SqrZ);
		SqrX = Inter(SqrX, SqrN-SqrZ-SqrY);
		Csqr(SqrZ, Z, -1);
		Csqr(SqrY, Y, -1);
		Csqr(SqrX, X, -1);
	}
}
//----------------------------------------------------------------------
void Cscal(interval& s, interval& ux, interval& uy, interval& vx, interval& vy)
{
	interval z1 = ux*vx;
	interval z2 = uy*vy;
	Cadd(s, z1, z2);
	Cmul(z2, uy, vy, -1);
	Cmul(z1, ux, vx, -1);
}
//----------------------------------------------------------------------
void Cscal(interval& s, double& ux, double& uy, interval& vx, interval& vy)
{
	interval z1 = ux*vx;
	interval z2 = uy*vy;
	Cadd(s, z1, z2);
	Cmul(z2, uy, vy, -1);
	Cmul(z1, ux, vx, -1);
}
//----------------------------------------------------------------------
void Cdet(interval& det, interval& ux, interval& uy, interval& vx, interval& vy, int sens)
{
	interval z1 = ux*vy;
	interval z2 = vx*uy;
	Csub(det, z1, z2, 1);
	if (sens == 1) return;
	Csub(det, z1, z2, -1);
	Cmul(z2, vx, uy, -1);
	Cmul(z1, ux, vy, -1);
}
//----------------------------------------------------------------------
void Cdet(interval& det, interval& ux, interval& uy, double& vx, double& vy, int sens)
{
	interval z1 = vy*ux;
	interval z2 = vx*uy;
	Csub(det, z1, z2, -1);
	if (sens == 1) return;
	Csub(det, z1, z2, -1);
	Cmul(z2, vx, uy, -1);
	Cmul(z1, vy, ux, -1);
}
//----------------------------------------------------------------------
void Cdet(interval& det, double& ux, double& uy, interval& vx, interval& vy, int sens)
{
	interval z1 = ux*vy;
	interval z2 = uy*vx;
	Csub(det, z1, z2, -1);
	if (sens == 1) return;
	Csub(det, z1, z2, -1);
	Cmul(z2, uy, vx, -1);
	Cmul(z1, ux, vy, -1);
}
//----------------------------------------------------------------------
void Cstep(interval& Y, interval& X)
{   
	vector<double> ax,ay,bx,by;
	ax.push_back(-1E10); ay.push_back(0); bx.push_back(0);    by.push_back(0);
	ax.push_back(0);     ay.push_back(0); bx.push_back(0);    by.push_back(1);
	ax.push_back(0);     ay.push_back(1); bx.push_back(1E10); by.push_back(1);
	CPointInSegments(X,Y,ax,ay,bx,by);
}
//----------------------------------------------------------------------
void CPointInSegment(interval& mx, interval& my, double ax, double ay, double bx, double by)
{      
	// contracte relativement a la contrainte : "m appartient au segment [a,b]"
	mx = Inter(mx, Union(interval(ax, ax), interval(bx, bx)));
	my = Inter(my, Union(interval(ay, ay), interval(by, by)));
	interval ma_x = ax - mx;
	interval ma_y = ay - my;
	double ab_x = bx - ax;
	double ab_y = by - ay;
	interval z1 = interval(0, 0);
	Cdet(z1, ab_x, ab_y, ma_x, ma_y, -1);
	Csub(ma_y, ay, my, -1);
	Csub(ma_x, ax, mx, -1);
	if ((mx.is_empty())||(my.is_empty())) { mx = interval::EMPTY_SET; my = interval::EMPTY_SET; }
}
//----------------------------------------------------------------------
void CPointInSegments(interval& mx, interval& my, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// contracte relativement a la contrainte : "m appartient au polygone dont les segments sont les [ai,bi]"
	if (ax.size() == 0) return;
	vector<interval> Mx(ax.size());
	vector<interval> My(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++)
	{
		interval mx0 = mx;
		interval my0 = my;
		CPointInSegment(mx0, my0, ax[j], ay[j], bx[j], by[j]);
		Mx[j] = mx0;
		My[j] = my0;
	}
	mx = Inter(mx, Union(Mx));
	my = Inter(my, Union(My));
}
//----------------------------------------------------------------------
class CtcInCircle : public Ctc {
public:
	CtcInCircle(double _cx, double _cy,double _r) : Ctc(2),
		C(dist()),cx(_cx),cy(_cy),r(_r) {}
	void contract(IntervalVector& P) 
	{
		IntervalVector p(3); p[0]=cx;  p[1]=cy; p[2]=r;
		IntervalVector Q = cart_prod(P,p);
		C.contract(Q);
		P = Q.subvector(0,1);
	}
	static NumConstraint& dist();
	CtcFwdBwd C;
	double cx,cy,r;
};
//----------------------------------------------------------------------
NumConstraint& CtcInCircle::dist() 
{   
	static NumConstraint dist1("xa","ya","xb","yb","r","(xa-xb)^2+(ya-yb)^2=r^2");
	return dist1;
}
//----------------------------------------------------------------------
void CPointInCircle(interval& X, interval& Y, double cx, double cy, double r)
{   
	CtcInCircle C(cx,cy,r);
	IntervalVector P(2); P[0]=X; P[1]=Y;
	try {C.contract(P);} catch(EmptyBoxException&) {X.set_empty(); Y.set_empty();};
	X=P[0];  Y=P[1];
}
//----------------------------------------------------------------------
//void Cnotin(interval& X, const interval& Y)
//{   
//	if (Y.is_empty()) return;
//	interval X1 = Inter(X, interval(-oo,Y.lb()));
//    interval X2 = Inter(X, interval(Y.ub(),+oo));
//    X = Inter(X, Union(X1,X2));
//}
//----------------------------------------------------------------------
void Cnotin(interval& X, interval& Y)
{ 
	if ((X.is_empty())||(Y.is_empty())) return;
	else
	{
		if (X.lb() <= Y.lb())
		{
			if ((X.ub() <= Y.ub())&&(Y.lb() <= X.ub())) X = interval(X.lb(), Y.lb());
			else return;
		}
		else
		{
			if (X.ub() <= Y.ub()) X = interval();
			else
			{
				if (X.lb() <= Y.ub()) X = interval(Y.ub(), X.ub());
				else return;
			}
		}
	}
}
//----------------------------------------------------------------------
void C_q_in(interval& x, int q, vector<interval>& y)
{  
	//Array<IntervalVector> V(y.size());   // ibex : pourquoi array et non pas vector ?
	//for (unsigned int i=0;i<y.size();i++)        // ibex : pas implemente sur les intervalles ?
	//{
	//	IntervalVector Vi(1,y[i]);
	//	V[i]=Vi;
	//}
	//x = (qinter(V,q))[0];
	Array<IntervalVector> V(y.size());
	for (unsigned int i = 0 ;i < y.size(); i++)
	{
		IntervalVector* Vi = new IntervalVector(1, y[i]);
		V.set_ref(i, *Vi);
	}
	x = (qinter(V, q))[0];
	for (unsigned int i = 0 ;i < y.size(); i++)
	{
		delete &(V[i]);
	}
}
//----------------------------------------------------------------------
void SinRing(interval& X, interval& Y, double cx, double cy, interval R, bool outer)
{   
	static Variable x, y, r;
	static NumConstraint C(x,y,r, pow(x-cx,2)+pow(y-cy,2)=pow(r,2));
	static CtcHC4 ctc(C);
	if (outer==true)
	{  
		IntervalVector P(3);  P[0]=X; P[1]=Y; P[2]=R;
		try {ctc.contract(P);} catch(EmptyBoxException) {X=interval::EMPTY_SET;Y=X; return;}
		X=P[0]; Y=P[1];
		return;
	}
	interval Xa,Ya,Xb,Yb;
	IntervalVector Pa(3);  Pa[0]=X; Pa[1]=Y; Pa[2]=interval(-1,R.lb());
	try {ctc.contract(Pa);} catch(EmptyBoxException) {Xa=interval::EMPTY_SET;Ya=Xa;}
	Xa=Pa[0]; Ya=Pa[1];
	IntervalVector Pb(3);  Pb[0]=X; Pb[1]=Y; Pb[2]=interval(R.ub(),POS_INFINITY);
	try {ctc.contract(Pb);} catch(EmptyBoxException) {Xb=interval::EMPTY_SET;Yb=Xb;}
	Xb=Pb[0]; Yb=Pb[1];
	X=Xa|Xb;
	Y=Ya|Yb;
}
//----------------------------------------------------------------------
