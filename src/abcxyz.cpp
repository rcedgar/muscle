#include "myutils.h"
#include "abcxyz.h"

void LogVec(const string &Msg, const vector<double> &v)
	{
	Log("%s(%.2f, %.2f, %.2f)\n", Msg.c_str(), v[X], v[Y], v[Z]);
	}

void LogMx(const string &Msg, const vector<vector<double> > &Mx)
	{
	Log("\n");
	Log("%-10.10s    x         y         z\n", Msg.c_str());
	for (int Motif = 0; Motif < 3; ++Motif)
		{
		Log("[%c]  ", "ABC"[Motif]);
		for (int Axis = 0; Axis < 3; ++Axis)
			{
			double v = Mx[Motif][Axis];
			Log("  %8.2f", v);
			}
		Log("\n");
		}

	double AB = GetDist3D(
	  Mx[A][X], Mx[A][Y], Mx[A][Z], 
	  Mx[B][X], Mx[B][Y], Mx[B][Z]);

	double AC = GetDist3D(
	  Mx[A][X], Mx[A][Y], Mx[A][Z], 
	  Mx[C][X], Mx[C][Y], Mx[C][Z]);

	double BC = GetDist3D(
	  Mx[B][X], Mx[B][Y], Mx[B][Z], 
	  Mx[C][X], Mx[C][Y], Mx[C][Z]);
	Log("AB=%.2f, BC=%.2f, AC=%.2f\n", AB, BC, AC);
	}

void GetIdentityMx(vector<vector<double> > &Mx)
	{
	Resize3x3(Mx);

	Mx[0][0] = 1;
	Mx[0][1] = 0;
	Mx[0][2] = 0;

	Mx[1][0] = 0;
	Mx[1][1] = 1;
	Mx[1][2] = 0;

	Mx[2][0] = 0;
	Mx[2][1] = 0;
	Mx[2][2] = 1;
	}

void CrossProduct(
  const vector<double> &a,
  const vector<double> &b,
  vector<double> &Prod)
	{
	Resize3(Prod);
	Prod[X] = a[Y]*b[Z] - a[Z]*b[Y];
	Prod[Y] = a[Z]*b[X] - a[X]*b[Z];
	Prod[Z] = a[X]*b[Y] - a[Y]*b[X];
	}

void XFormMx(
  const vector<vector<double> > &Mx,
  const vector<double> &t,
  const vector<vector<double> > &R,
  vector<vector<double> > &XMx)
	{
	Resize3x3(XMx);
	for (uint i = 0; i < 3; ++i)
		XFormPt(Mx[i], t, R, XMx[i]);
	}

void XFormPt(
  const vector<double> &Pt,
  const vector<double> &t,
  const vector<vector<double> > &R,
  vector<double> &XPt)
	{
	Resize3(XPt);

	vector<double> Tmp;
	Resize3(Tmp);
	for (uint i = 0; i < 3; ++i)
		Tmp[i] = Pt[i] + t[i];

	for (uint i = 0; i < 3; ++i)
//		XPt[i] = R[i][0]*Tmp[0] + R[i][1]*Tmp[1] + R[i][2]*Tmp[2];
		XPt[i] = R[0][i]*Tmp[0] + R[1][i]*Tmp[1] + R[2][i]*Tmp[2];

	}

void XFormPt(
  const vector<float> &Pt,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<float> &XPt)
	{
	Resize3(XPt);

	vector<float> Tmp;
	Resize3(Tmp);
	for (uint i = 0; i < 3; ++i)
		Tmp[i] = Pt[i] + t[i];

	for (uint i = 0; i < 3; ++i)
//		XPt[i] = R[i][0]*Tmp[0] + R[i][1]*Tmp[1] + R[i][2]*Tmp[2];
		XPt[i] = R[0][i]*Tmp[0] + R[1][i]*Tmp[1] + R[2][i]*Tmp[2];

	}

void RotateMx(const vector<vector<double> > &Mx,
  uint Axis, double Theta, vector<vector<double> > &RotatedMx)
	{
	RotatedMx.clear();
	RotatedMx.resize(3);

	uint OtherAxis_i = GetOtherAxis_i(Axis);
	uint OtherAxis_j = GetOtherAxis_j(Axis);

	const double cos_Theta = cos(Theta);
	const double sin_Theta = sin(Theta);

	for (uint k = 0; k < 3; ++k)
		{
		RotatedMx[k].resize(3);

		RotatedMx[k][Axis] = Mx[k][Axis];

		double Coord_i = Mx[k][OtherAxis_i];
		double Coord_j = Mx[k][OtherAxis_j];

		RotatedMx[k][OtherAxis_i] = Coord_i*cos_Theta - Coord_j*sin_Theta;
		RotatedMx[k][OtherAxis_j] = Coord_i*sin_Theta + Coord_j*cos_Theta;
		}

	AssertSameLengths(Mx, RotatedMx);
	AssertSameAngles(Mx, RotatedMx);
	}

double dot(const double a[3], const double b[3])
	{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	}

void transform(const double t[3], const double u[3][3],
  const double x[3], double x_transformed[3])
	{
	x_transformed[0] = t[0] + dot(u[0], x);
	x_transformed[1] = t[1] + dot(u[1], x);
	x_transformed[2] = t[2] + dot(u[2], x);
	}
