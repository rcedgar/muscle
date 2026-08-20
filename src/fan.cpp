#include "myutils.h"
#include "fan.h"

const char *FAN2str(FAN fan)
	{
	switch (fan)
		{
#define f(x)	case FAN_##x: return #x;
#include "flat_featlist.h"
		}
	return "FAN_invalid";
	}

FAN str2FAN(const char *s)
	{
	if (0) ;
#define f(x)	else if (strcmp(s, #x) == 0) return FAN_##x;
#include "flat_featlist.h"
	Die("str2FA(%s)",s);
	return FAN_COUNT;
	}

FAN str2FAN(const string &s)
	{
	return str2FAN(s.c_str());
	}

bool is_quantized(FAN fan)
	{
	switch (fan)
		{
	case FAN_aa:
	case FAN_pm:
	case FAN_sec:
	case FAN_nensec:
	case FAN_rensec:
	case FAN_pensec:
	case FAN_mensec:
	case FAN_nu:
	case FAN_kappa:
		return false;

	case FAN_nendist:
	case FAN_rendist:
	case FAN_pendist:
	case FAN_mendist:
	case FAN_fendist:
	case FAN_pmdd:
	case FAN_pmdiff:
	case FAN_pack:
	case FAN_ppack:
	case FAN_mpack:
	case FAN_angle:
	case FAN_turnd:
		return true;
		}

	Die("is_quantized(%u=%s)", uint(fan), FAN2str(fan));
	return false;
	}
