#pragma once

typedef void OnPathGlobal_fn(const string &Path);
typedef void OnPathLocal_fn(uint PosA, uint PosB, const string &Path);
void EnumPathsLocal(uint LA, uint LB, OnPathLocal_fn OnPath);