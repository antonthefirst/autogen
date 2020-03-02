#pragma once

#define CELL_COUNT_MAX 8
#define SITE_COUNT_MAX 4
#define BODY_TYPE_COUNT_MAX 32 // hard limit, because of the "types inside" bitmask variable.
#define NO_TYPE u8(-1)
#define MOL_RAD 1.f
#define BIND_POS_THRESH 0.9f

enum { GENDER_BOTH, GENDER_FEMALE, GENDER_MALE };

#include "core/type_magic.h"
#define TYPE() SHELL_STATUS
#define VALUES(X) \
X(empty) \
X(capsid) \
X(inert) \
X(dormant) \
X(fertile) \
X(ripe) 
DECLARE_ENUM(TYPE(), VALUES)
#undef VALUES
#undef TYPE

#define NO_IDX (-1)

#define PATH "data/autogen/"

#include "struct_pre.inl"
#include "autogen.typ"
#include "struct_post.inl"

#include "serialize_pre.inl"
#include "autogen.typ"
#include "serialize_post.inl"

#include "gui_pre.inl"
#include "autogen.typ"
#include "gui_post.inl"
