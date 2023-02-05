
#include "SurfaceFactory.h"

SurfaceFactory& surfaceFactory()
{
	static SurfaceFactory reg;
	return reg;
}
