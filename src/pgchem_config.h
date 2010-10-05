/************************************************************************
 * pgchem_config.h compile-time configuration settings
 *
 * Copyright (c) 2004 by Ernst-G. Schmid
 * 
 * This file is part of the xchem::tigress project.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * lesser GNU General Public License for more details.
 ************************************************************************/

#ifdef WIN32
#define PGCHEM_VERSION "pgchem_tigress_1.3 GiST WIN32"
//#define PGCHEM_SDF_SEPARATOR "\r\n$$$$\r\n"
//#define TDFILE_PATH NULL
#else
#define PGCHEM_VERSION "pgchem_tigress_1.3 GiST POSIX"
//#define PGCHEM_INJECTOR "echo \""
//#define PGCHEM_SDF_SEPARATOR "\n\\$\\$\\$\\$\n"
#endif
//#define MAX_MOL_SIZE 262143   //Max. supported molfile size in byte
//#define HAS_INTEL_SSE         //Processor is Intel SSE capable
