#pragma once

/***
Coordinates are stored as uint16_t "ic_t".
	L rows x 3 columns
	Same format as BCA
	x0,y0,z0 x1,y1,z1 ... xL-1,yL-1,zL-1

Distances are stored as uint16_t "sid_t"
	sid	= (ic*ic)/16
		= (d*d*100)/16

Convert sid_t to float in Angstroms
	d = sqrt(16*sid/100)

Convert sid_t to ic_t in 1/10th Angstroms
	d_ic = int(10*sqrt(16*sid/100) + 0.5)

    (maxic*maxic)/16 = 65535
	maxic = sqrt(65535*16) = 1024
	maximum distance = 102.4 Angstroms
	resolution ~0.4 Angstroms distance

Squared distance sd(i,j) stored for all 1 < |i-j| <= M
    M pairs (j values) for typical i
    <M pairs close to the ends

Flat distamx matrix layout:
    distmx[M*i + j - i - 1] where j = i+1, i+2 ... i+M
***/

using ic_t = uint16_t;	// 1/10th Angstrom units
using sid_t = uint16_t;	// (ic_t*ic_t*)/8

using p_uint8_t = uint8_t * __restrict;
using cp_uint8_t = const uint8_t * __restrict;
using p_ic_t = uint16_t * __restrict;
using cp_ic_t = const uint16_t * __restrict;
using p_sid_t = sid_t * __restrict;
using cp_sid_t = const sid_t * __restrict;
using p_uint16_t = uint16_t * __restrict;
using cp_uint16_t = const uint16_t * __restrict;

static inline ic_t coord2ic(float x) { return ic_t((x + 1000)*10 + 0.5); }
static inline float ic2coord(ic_t ic) { return float(ic/10.0f) - 1000; }

static inline sid_t dist2sid(float d)
	{
	return sid_t(d*d*(100.0f/16) + 0.5f);
	}

static inline float sid2dist(sid_t sid)
	{
	return sqrtf(16.0f*sid)/10.0f;
	}

static inline float sid2dist2(sid_t sid)
	{
	return 16.0f*sid/100.0f;
	}

static inline sid_t icxyzpair2sid(
	ic_t x1, ic_t y1, ic_t z1,
	ic_t x2, ic_t y2, ic_t z2)
	{
	int32_t dx = int32_t(x1) - int32_t(x2);
	int32_t dy = int32_t(y1) - int32_t(y2);
	int32_t dz = int32_t(z1) - int32_t(z2);
	sid_t sid = (dx*dx + dy*dy + dz*dz)/16;
	return sid; 
	}

extern const ic_t sid2ic[65536];
