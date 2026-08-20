#pragma once

double dot(const double a[3], const double b[3]);
void transform(const double t[3], const double u[3][3],
  const double x[3], double x_transformed[3]);

static const double PI = 3.1415926535;
static const int X = 0;
static const int Y = 1;
static const int Z = 2;

static const int A = 0;
static const int B = 1;
static const int C = 2;

const uint AL = 12;
const uint BL = 14;
const uint CL = 8;
const uint PALMCOREM = 150;

extern const uint MotifLVec[3];

void LogVec(const string &Msg, const vector<double> &v);
void LogMx(const string &Msg, const vector<vector<double> > &Mx);
double GetMxDeterminant(const vector<vector<double> > &Mx);
void InvertMx(const vector<vector<double> > &Mx,
  vector<vector<double> > &InvMx);
void MulMxVec(const vector<vector<double> > &Mx,
  const vector<double> &Vec,
  vector<double> &Result);
void MulMx(
  const vector<vector<double> > &A,
  const vector<vector<double> > &B,
  vector<vector<double> > &Prod);
void GetBasisR(const vector<vector<double> > &Basis,
  vector<vector<double> > &R);
void LogMx(const string &Msg, const vector<vector<double> > &Mx);
void RotateMx(const vector<vector<double> > &Mx,
  uint Axis, double Theta, vector<vector<double> > &RotatedMx);
void CrossProduct(
  const vector<double> &a,
  const vector<double> &b,
  vector<double> &Prod);
void GetTriForm(
  vector<vector<double> > &MotifCoords,
  vector<double> &t,
  vector<vector<double> > &R);

void GetTriangleCentroid(const vector<vector<double> > &MotifCoords,
  vector<double> &CentroidCoords);
void GetTriangleBasis(const vector<vector<double> > &MotifCoords,
  vector<double> &CentroidCoords, vector<vector<double> > &Basis);

static inline uint GetOtherAxis_i(uint Axis) { return (Axis+1)%3;  }
static inline uint GetOtherAxis_j(uint Axis) { return (Axis+2)%3;  }

template <typename T> void Resize3(vector<T> &v) { v.resize(3); }

template <typename T> void Resize3x3(vector<vector<T> > &Mx)
	{
	Mx.resize(3);
	for (uint i = 0; i < 3; ++i)
		Mx[i].resize(3);
	}

template <typename T> void MulVecScalar(const vector<T> &v,
  T s, vector<T> &r)
	{
	asserta(SIZE(v) == 3);
	r.resize(3);
	r[0] = v[0]*s;
	r[1] = v[1]*s;
	r[2] = v[2]*s;
	}

template <typename T> T GetDist2(
  const vector<T> &Pt1,
  const vector<T> &Pt2)
	{
	assert(Pt1.size() == 3);
	assert(Pt2.size() == 3);
	double dx = Pt1[X] - Pt2[X];
	double dy = Pt1[Y] - Pt2[Y];
	double dz = Pt1[Z] - Pt2[Z];
	double d2 = dx*dx + dy*dy + dz*dz;
	return d2;
	}

template <typename T> T GetDist(
  const vector<T> &Pt1,
  const vector<T> &Pt2)
	{
	assert(Pt1.size() == 3);
	assert(Pt2.size() == 3);
	double dx = Pt1[X] - Pt2[X];
	double dy = Pt1[Y] - Pt2[Y];
	double dz = Pt1[Z] - Pt2[Z];
	double d2 = dx*dx + dy*dy + dz*dz;
	double d = sqrt(d2);
	return d;
	}

static inline double GetDist3D(
  double x1, double y1, double z1,
  double x2, double y2, double z2)
	{
	double dx = x1 - x2;
	double dy = y1 - y2;
	double dz = z1 - z2;
	double d2 = dx*dx + dy*dy + dz*dz;
	double d = sqrt(d2);
	return d;
	}

static inline float GetDist3D(
  float x1, float y1, float z1,
  float x2, float y2, float z2)
	{
	float dx = x1 - x2;
	float dy = y1 - y2;
	float dz = z1 - z2;
	float d2 = dx*dx + dy*dy + dz*dz;
	float d = sqrt(d2);
	return d;
	}

static inline double GetDist2_Mxij(const vector<vector<double> > &Mx,
  uint i, uint j)
	{
	double xi = Mx[i][X];
	double yi = Mx[i][Y];
	double zi = Mx[i][Z];

	double xj = Mx[j][X];
	double yj = Mx[j][Y];
	double zj = Mx[j][Z];

	double dx = xi - xj;
	double dy = yi - yj;
	double dz = zi - zj;

	double d2 = dx*dx + dy*dy + dz*dz;
	return d2;
	}

static inline double GetDist_Mxij(const vector<vector<double> > &Mx,
  uint i, uint j)
	{
	double xi = Mx[i][X];
	double yi = Mx[i][Y];
	double zi = Mx[i][Z];

	double xj = Mx[j][X];
	double yj = Mx[j][Y];
	double zj = Mx[j][Z];

	double dx = xi - xj;
	double dy = yi - yj;
	double dz = zi - zj;

	double d2 = dx*dx + dy*dy + dz*dz;
	double d = sqrt(d2);
	return d;
	}

template <typename T> T GetMod_xyz(T x, T y, T z)
	{
	T d2 = x*x + y*y + z*z;
	T Mod = sqrt(d2);
	return Mod;
	}

template <typename T> T GetMod_Vec(const vector<T> &v)
	{
	T x = v[X];
	T y = v[Y];
	T z = v[Z];
	T Mod = GetMod_xyz(x, y, z);
	return Mod;
	}

template <typename T> void NormalizeVec(vector<T> &v)
	{
	T Mod = GetMod_Vec(v);
	assert(Mod > 0);
	v[X] /= Mod;
	v[Y] /= Mod;
	v[Z] /= Mod;
	assert(feq(GetMod_Vec(v), 1));
	}

template <typename T> T GetMod_Mxi(const vector<vector<T> > &Mx,
  uint i)
	{
	T x = Mx[i][X];
	T y = Mx[i][Y];
	T z = Mx[i][Z];
	T Mod = GetMod_xyz(x, y, z);
	return Mod;
	}

template <typename T> T GetTheta3D(
  T xi, T yi, T zi,
  T xj, T yj, T zj)
	{
	T DotProd = xi*xj + yi*yj + zi*zj;
	T Modi = GetMod_xyz(xi, yi, zi);
	T Modj = GetMod_xyz(xj, yj, zj);
	if (abs(Modi*Modj) < 1e-6)
		return 0;
	T cos_theta = DotProd/(Modi*Modj);
	asserta(cos_theta >= -1.02 && cos_theta <= 1.02);
	if (cos_theta < -1)
		cos_theta = -1;
	else if (cos_theta > 1)
		cos_theta = 1;
	T theta = acos(cos_theta);
	return theta;
	}

// angle a-b-c where b is the vertex
template <typename T> T GetTheta3D_3pts(
	T xa, T ya, T za,
	T xb, T yb, T zb,
	T xc, T yc, T zc)
	{
    return GetTheta3D(
        xa - xb, ya - yb, za - zb,
        xc - xb, yc - yb, zc - zb);
	}

template <typename T> T GetTheta3D(
  const vector<T> &vi,
  const vector<T> &vj)
	{
	asserta(SIZE(vi) == 3);
	asserta(SIZE(vj) == 3);
	return GetTheta3D(vi[0], vi[1], vi[2], vj[0], vj[1], vj[2]);
	}

template <typename T> void Sub_Vecs(const vector<T> &vi,
  const vector<T> &vj, vector<T> &Diff)
	{
	Resize3(Diff);
	Diff[X] = vi[X] - vj[X];
	Diff[Y] = vi[Y] - vj[Y];
	Diff[Z] = vi[Z] - vj[Z];
	}

template <typename T> void Add_Vecs(const vector<T> &vi,
  const vector<T> &vj, vector<T> &Sum)
	{
	Resize3(Sum);
	Sum[X] = vi[X] + vj[X];
	Sum[Y] = vi[Y] + vj[Y];
	Sum[Z] = vi[Z] + vj[Z];
	}

template <typename T> T GetTheta_Vecs(const vector<T> &vi,
  const vector<T> &vj)
	{
	T xi = vi[X];
	T yi = vi[Y];
	T zi = vi[Z];

	T xj = vj[X];
	T yj = vj[Y];
	T zj = vj[Z];

	T theta = GetTheta3D(xi, yi, zi, xj, yj, zj);
	return theta;
	}

static inline double GetTheta_Mxij(const vector<vector<double> > &Mx,
  uint i, uint j)
	{
	double theta = GetTheta_Vecs(Mx[i], Mx[j]);
	return theta;
	}

template <typename T> T degrees(T Radians) { return Radians*180/T(PI); }

template <typename T> T degrees_0_to_360(T Radians)
	{
	T Deg = Radians*180.0/PI;
	Deg = fmod(Deg, 360);
	if (Deg < 0)
		Deg += 360;
	assert(Deg >= 0 && Deg < 360);
	return Deg;
	}

void GetIdentityMx(vector<vector<double> > &Mx);

void XFormPt(
  const vector<double> &Pt,
  const vector<double> &t,
  const vector<vector<double> > &R,
  vector<double> &XPt);

void XFormMx(
  const vector<vector<double> > &Mx,
  const vector<double> &t,
  const vector<vector<double> > &R,
  vector<vector<double> > &XMx);

void XFormPt(
  const vector<float> &Pt,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<float> &XPt);

void XFormMx(
  const vector<vector<float> > &Mx,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<vector<float> > &XMx);

#if DEBUG
static void AssertCanonicalUnitBasis(const vector<vector<double> > &Basis)
	{
	assert(SIZE(Basis) == 3);
	assert(SIZE(Basis[0]) == 3);
	assert(SIZE(Basis[1]) == 3);
	assert(SIZE(Basis[2]) == 3);

	assert(feq(Basis[0][X], 1));
	assert(feq(Basis[0][Y], 0));
	assert(feq(Basis[0][Z], 0));

	assert(feq(Basis[1][X], 0));
	assert(feq(Basis[1][Y], 1));
	assert(feq(Basis[1][Z], 0));

	assert(feq(Basis[2][X], 0));
	assert(feq(Basis[2][Y], 0));
	assert(feq(Basis[2][Z], 1));
	}
#else
#define AssertCanonicalUnitBasis(x)	0
#endif // DEBUG

static void AssertUnitBasisA(const vector<vector<double> > &Basis)
	{
// check length of basis vectors is 1 and angles are PI/2
	for (uint k = 0; k < 3; ++k)
		{
		double Modk = GetMod_Mxi(Basis, k);
		assert(Modk >= 0.98 && Modk <= 1.02);

		uint Axis_i = GetOtherAxis_i(k);
		uint Axis_j = GetOtherAxis_j(k);

		double theta = GetTheta_Mxij(Basis, Axis_i, Axis_j);
		assert(feq(theta, PI/2));
		}
	}

#if DEBUG
static void AssertUnitBasis(const vector<vector<double> > &Basis)
	{
// check length of basis vectors is 1 and angles are PI/2
	for (uint k = 0; k < 3; ++k)
		{
		double Modk = GetMod_Mxi(Basis, k);
		assert(Modk >= 0.98 && Modk <= 1.02);

		uint Axis_i = GetOtherAxis_i(k);
		uint Axis_j = GetOtherAxis_j(k);

		double theta = GetTheta_Mxij(Basis, Axis_i, Axis_j);
		assert(feq(theta, PI/2));
		}
	}
#else
#define AssertUnitBasis(x)	0
#endif // DEBUG

#if DEBUG
static void AssertSameLengths(const vector<vector<double> > &Mx1,
  const vector<vector<double> > &Mx2)
	{
	for (uint k = 0; k < 3; ++k)
		{
		double Modk1 = GetMod_Mxi(Mx1, k);
		double Modk2 = GetMod_Mxi(Mx2, k);
		assert(feq(Modk1, Modk2));
		}
	}
#else
#define AssertSameLengths(x, y)	0
#endif // DEBUG

#if DEBUG
static void AssertSameAngles(const vector<vector<double> > &Mx1,
  const vector<vector<double> > &Mx2)
	{
	for (uint i = 0; i < 3; ++i)
		{
		for (uint j = i + 1; j < 3; ++j)
			{
			double Theta1 = GetTheta_Mxij(Mx1, i, j);
			double Theta2 = GetTheta_Mxij(Mx2, i, j);
			asserta(feq(Theta1, Theta2));
			}
		}
	}
#else
#define AssertSameAngles(x, y)	0
#endif // DEBUG

#if DEBUG
static void AssertMx3D(const vector<vector<double> > &Mx)
	{
	asserta(SIZE(Mx) == 3);
	for (uint i = 0; i < 3; ++i)
		asserta(SIZE(Mx[i]) == 3);
	}
#else
#define AssertMx3D(x)	0
#endif
