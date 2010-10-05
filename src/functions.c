/************************************************************************
 * functions.c native chemistry handling functions
 *  
 * Copyright (c) 2004,2009 by Ernst-G. Schmid
 * Copyright (c) 2004,2005 by Bayer Business Services GmbH
 * for explicitly marked functions
 *   
 * This file is part of the xchem::tigress project.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the lesser GNU General Public License as published by
 * the Free Software Foundation version 2.1 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * lesser GNU General Public License for more details.
 ************************************************************************/

#ifndef WIN32
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <postgres.h>
#include <utils/varbit.h>
#include <fmgr.h>
#include "libpq/md5.h"
#include "functions.h"
#include "obwrapper.h"
#include "barsoi/barsoi.h"
#include "reaction/reaction.h"

PG_MODULE_MAGIC;

/*static void
bytesToHex(uint8 b[16], char *s)
{
	static const char *hex = "0123456789abcdef";
	int			q,
				w;

	for (q = 0, w = 0; q < 16; q++)
	{
		s[w++] = hex[(b[q] >> 4) & 0x0F];
		s[w++] = hex[b[q] & 0x0F];
	}
	s[w] = '\0';
}*/

/*static void binary_to_01(int bin, char *str)
{
    unsigned int mask;      // used to check each individual bit, 
                            //    unsigned to alleviate sign 
                            //    extension problems

    mask = 0x80000000;      // Set only the high-end bit
    while (mask)            // Loop until MASK is empty
    {
        if (bin & mask)     // test the masked bit
              *str = '1';   // if true, value is 1
          else 
              *str = '0';   // if false, value is 0
        str++;              // next character
        mask >>= 1;         // shift the mask 1 bit
    }
    *str = 0;               // add the trailing null 
}*/

//char molfile_buf[MAX_MOL_SIZE + 1] = "";

//char current_molfile[MAX_MOL_SIZE + 1] = "";

//pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

/*#define CSA(h,l, a,b,c) \
   {unsigned int u = a ^ b; unsigned int v = c; \
      h = (a & b) | (u & v); l = u ^ v;}

inline static int
countonbits (unsigned x)
{
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x + (x >> 4)) & 0x0F0F0F0F;
  x = x + (x >> 8);
  x = x + (x >> 16);
  return x & 0x0000003F;
}

inline static int
countonbitsinarray (unsigned A[], int n)
{

  int tot, i;
  unsigned int ones, twos, twosA, twosB, fours, foursA, foursB, eights;

  tot = 0;			// Initialize.
  fours = twos = ones = 0;

  for (i = 0; i <= n - 8; i = i + 8)
    {
      CSA (twosA, ones, ones, A[i], A[i + 1])
	CSA (twosB, ones, ones, A[i + 2], A[i + 3])
	CSA (foursA, twos, twos, twosA, twosB)
	CSA (twosA, ones, ones, A[i + 4], A[i + 5])
	CSA (twosB, ones, ones, A[i + 6], A[i + 7])
	CSA (foursB, twos, twos, twosA, twosB)
	CSA (eights, fours, fours, foursA, foursB)
	tot = tot + countonbits (eights);
    }
  tot =
    8 * tot + 4 * countonbits (fours) + 2 * countonbits (twos) +
    countonbits (ones);

  for (i = i; i < n; i++)	// Simply add in the last
    tot = tot + countonbits (A[i]);	// 0 to 7 elements.
  return tot;
}*/


//inline static bool
//match_exact (bytea * molfile_needle, bytea * molfile_haystack,
//           bool e_z_global, bool r_s_global)
//{
//
//  FILE *process;
//  char *rvbuffer;
//  bool retval = false;
//  int eno;
//  int modepos = 3;
//  cow_orker_info *cw_info;
//  const int rvbuffer_len = 5;
//  char match_mode[6];
//#ifndef WIN32
//  char *sdfile;
//#endif
//
//  int32 datalen_molfile_needle = VARSIZE (molfile_needle) - VARHDRSZ;
//  int32 datalen_molfile_haystack = VARSIZE (molfile_haystack) - VARHDRSZ;
//
//  match_mode[0] = '-';
//  match_mode[1] = 'x';
//  match_mode[2] = 's';
//
//  if (e_z_global)
//    {
//      match_mode[modepos++] = 'g';
//    }
//
//  if (r_s_global)
//    {
//      match_mode[modepos++] = 'G';
//    }
//
//  match_mode[modepos] = '\0';
//
//#ifdef WIN32
//
//  cw_info =
//    spawn_cow_orker ("matchmol", match_mode, VARDATA (molfile_needle),
//                   datalen_molfile_needle, VARDATA (molfile_haystack),
//                   datalen_molfile_haystack);
//
//  process = cw_info->process;
//#else
//  sdfile =
//    (char *) palloc (datalen_molfile_needle + strlen (PGCHEM_SDF_SEPARATOR) +
//                   datalen_molfile_haystack + 1);
//
//  memcpy (sdfile, VARDATA (molfile_needle), datalen_molfile_needle);
//  sdfile[datalen_molfile_needle] = '\0';      //terminate properly
//  strcat (sdfile, PGCHEM_SDF_SEPARATOR);
//  strncat (sdfile, VARDATA (molfile_haystack), datalen_molfile_haystack);
//
//  cw_info =
//    spawn_cow_orker ("matchmol", match_mode, (const char *) sdfile,
//                   strlen (sdfile));
//
//  process = cw_info->process;
//
//  //pfree (sdfile);
//#endif
//
//  if (process != NULL)
//    {
//      rvbuffer = (char *) palloc (rvbuffer_len);
//
//      fgets (rvbuffer, rvbuffer_len, process);
//      retval = (bool) (strstr (rvbuffer, "1:T") != NULL);
//
//      //pfree (rvbuffer);
//    }
//  else
//    {
//      ereport (ERROR, (errmsg ("matchmol start failed")));
//    }
//
//  eno = cow_orker_cleanup (cw_info);
//
//  if (eno == 127)
//    {
//      ereport (ERROR, (errmsg ("matchmol start failed: 127")));
//    }
//
//  /* return as psql type by reference */
//  return retval;
//}
//
//
//inline static bool
//match_substruct (bytea * molfile_needle, bytea * molfile_haystack,
//               bool strict_typing, bool e_z_global, bool r_s_global)
//{
//
//  FILE *process;
//  char *rvbuffer;
//  char match_mode[5];
//  bool retval = false;
//  int eno;
//  int modepos = 0;
//  cow_orker_info *cw_info;
//  const int rvbuffer_len = 5;
//#ifndef WIN32
//  char *sdfile;
//#endif
//
//  int32 datalen_molfile_needle = VARSIZE (molfile_needle) - VARHDRSZ;
//  int32 datalen_molfile_haystack = VARSIZE (molfile_haystack) - VARHDRSZ;
//
//  if (!(strict_typing && e_z_global && r_s_global))
//    {
//      match_mode[modepos] = '\0';
//    }
//  else
//    {
//      match_mode[modepos++] = '-';
//      if (strict_typing)
//      {
//        match_mode[modepos++] = 's';
//      }
//
//      if (e_z_global)
//      {
//        match_mode[modepos++] = 'g';
//      }
//
//      if (r_s_global)
//      {
//        match_mode[modepos++] = 'G';
//      }
//
//      match_mode[modepos] = '\0';
//    }
//
//#ifdef WIN32
//  cw_info =
//    spawn_cow_orker ("matchmol", match_mode, VARDATA (molfile_needle),
//                   datalen_molfile_needle, VARDATA (molfile_haystack),
//                   datalen_molfile_haystack);
//
//  process = cw_info->process;
//#else
//  sdfile =
//    (char *) palloc (datalen_molfile_needle + strlen (PGCHEM_SDF_SEPARATOR) +
//                   datalen_molfile_haystack + 1);
//
//  memcpy (sdfile, VARDATA (molfile_needle), datalen_molfile_needle);
//  sdfile[datalen_molfile_needle] = '\0';      //terminate properly
//  strcat (sdfile, PGCHEM_SDF_SEPARATOR);
//  strncat (sdfile, VARDATA (molfile_haystack), datalen_molfile_haystack);
//
//  cw_info =
//    spawn_cow_orker ("matchmol", match_mode, (const char *) sdfile,
//                   strlen (sdfile));
//
//  process = cw_info->process;
//
//  //pfree (sdfile);
//#endif
//
//  if (process != NULL)
//    {
//      rvbuffer = (char *) palloc (rvbuffer_len);
//
//      fgets (rvbuffer, rvbuffer_len, process);
//      retval = (bool) (strstr (rvbuffer, "1:T") != NULL);
//
//      //pfree (rvbuffer);
//    }
//  else
//    {
//      ereport (ERROR, (errmsg ("matchmol start failed")));
//    }
//
//  eno = cow_orker_cleanup (cw_info);
//
//  if (eno == 127)
//    {
//      ereport (ERROR, (errmsg ("matchmol start failed: 127")));
//    }
//
//  /* return as psql type by reference */
//  return retval;
//}

//inline static bool
//is_checkmol_molfile (bytea * arg_molfile)
//{
//  bool retval = FALSE;
//  char *molfile = VARDATA (arg_molfile);
//
//  retval = (strstr (molfile, "CheckMol") != NULL);
//
//  return retval;
//}

//PG_FUNCTION_INFO_V1 (pgchem_is_fp_ss_candidate);

/*Datum
pgchem_is_fp2_ss_candidate_SSE (PG_FUNCTION_ARGS)
{
  bytea *fp2_sarg = PG_GETARG_BYTEA_P (0);
  bytea *fp2_test = PG_GETARG_BYTEA_P (1);
  unsigned int i,j,k;
  DQWORD result;
  DQWORD *fpdata_sarg;
  DQWORD *fpdata_test;
  
  j = VARSIZE (fp2_sarg) - VARHDRSZ;
  k = VARSIZE (fp2_test) - VARHDRSZ;

  if (j < 1 || k < 1 || j != k)
    {
      ereport (ERROR,
	       (errmsg
		("inconsistent fingerprint sizes read. aborting query.")));
    }
  else
    {
      fpdata_sarg = (DQWORD *) VARDATA (fp2_sarg);
      fpdata_test = (DQWORD *) VARDATA (fp2_test);
      
      j=0;
      
for(i=0;i<2;i++) {
               
     __asm__ ("movups %1, %%xmm0;" // load 16 byte of search argument from [j] into xmm0
              "movups %2, %%xmm1;" // load 16 byte of test argument from [j] into xmm1
              "movups %3, %%xmm2;" // load 16 byte of search argument from [j+1] into xmm2
              "movups %4, %%xmm3;" // load 16 byte of test argument from [j+1] into xmm3
              "movups %5, %%xmm4;" // and so on...
              "movups %6, %%xmm5;"
              "movups %7, %%xmm6;"
              "movups %8, %%xmm7;"
              
              "andps %%xmm0, %%xmm1;" // parallel AND xmm0 with xmm1 and store the result in xmm1
              "andps %%xmm2, %%xmm3;" // ...
              "andps %%xmm4, %%xmm5;"
              "andps %%xmm6, %%xmm7;"
              
              "xorps %%xmm0, %%xmm1;" // parallel XOR xmm0 with xmm1 and store the result in xmm1
              "xorps %%xmm2, %%xmm3;" // ...
              "xorps %%xmm4, %%xmm5;"
              "xorps %%xmm6, %%xmm7;"
              
              "orps %%xmm7, %%xmm1;" // parallel OR xmm7, xmm5 and xmm3 with xmm1
              "orps %%xmm5, %%xmm1;" // so we get a single 16 byte result in xmm1
              "orps %%xmm3, %%xmm1;"
              
              "movups %%xmm1, %0;" // store the 16 byte from xmm1 into the 'result' DQWORD structure
              
               : "=m" (result)
               : "m" (fpdata_sarg[j]), "m" (fpdata_test[j]), "m" (fpdata_sarg[j+1]), "m" (fpdata_test[j+1]), "m" (fpdata_sarg[j+2]), "m" (fpdata_test[j+2]), "m" (fpdata_sarg[j+3]), "m" (fpdata_test[j+3])
               : "xmm0", "xmm1", "xmm2", "xmm3","xmm4", "xmm5", "xmm6", "xmm7"
               ); 
               
               if(result.a !=0 || result.b !=0) return false; // if one of the two members of 'result' is != 0 we have a mismatch.
               
               j+=4;
         }  
    }

  return true;
} */

/*#ifdef HAS_INTEL_SSE

Datum
pgchem_is_fp_ss_candidate (PG_FUNCTION_ARGS)
{
  MOLFP2 *fp2_sarg = PG_GETARG_MOLFP2_P (0);
  MOLFP2 *fp2_test = PG_GETARG_MOLFP2_P (1);
  int i;
  DQWORD result __attribute__ ((aligned (16)));
  DQWORD fpdata_sarg[8] __attribute__ ((aligned (16)));
  DQWORD fpdata_test[8] __attribute__ ((aligned (16)));


  //fpdata_sarg = (DQWORD *) VARDATA (fp2_sarg);
  //fpdata_test = (DQWORD *) VARDATA (fp2_test);

  memcpy (&fpdata_sarg, fp2_sarg, FPSIZE2);
  memcpy (&fpdata_test, fp2_test, FPSIZE2);

  for (i = 0; i < (FPSIZE2 / 16); i++)
    {

      __asm__ ("movaps %1, %%xmm0;"	// load 16 byte of search argument from [i] into xmm0
	       "movaps %2, %%xmm1;"	// load 16 byte of test argument from [i] into xmm1
	       "andps %%xmm0, %%xmm1;"	// parallel AND xmm0 with xmm1 and store the result in xmm1
	       "xorps %%xmm0, %%xmm1;"	// parallel XOR xmm0 with xmm1 and store the result in xmm1
	       "movaps %%xmm1, %0;"	// store the 16 byte from xmm1 into the 'result' DQWORD structure
    : "=m" (result): "m" (fpdata_sarg[i]), "m" (fpdata_test[i]):"xmm0", "xmm1",
	       "memory");

      if (result.a != 0 || result.b != 0)
	return false;		// if one of the two members of 'result' is != 0 we have a mismatch.
    }


  return true;
}

#else

Datum
pgchem_is_fp_ss_candidate (PG_FUNCTION_ARGS)
{
  MOLFP2 fp2_sarg = PG_GETARG_MOLFP2_P (0);
  MOLFP2 fp2_test = PG_GETARG_MOLFP2_P (1);

  int i;

  for (i = 0; i < FPSIZE2; i++)
    {
      if (((*fp2_sarg & *fp2_test++) ^ *fp2_sarg++) != 0)
	return false;
    }

}

return true;
}

#endif */

/*
* Returns the version of pgchem.
*/
PG_FUNCTION_INFO_V1 (pgchem_version);

Datum
pgchem_version (PG_FUNCTION_ARGS)
{
  PG_RETURN_CSTRING (PGCHEM_VERSION);
}

/*
* Returns the version of barsoi.
*/
PG_FUNCTION_INFO_V1 (pgchem_barsoi_version);

Datum
pgchem_barsoi_version (PG_FUNCTION_ARGS)
{
  char *barsoi_version_buffer;

  barsoi_version_buffer = (char *) palloc (256);

  xm_version (barsoi_version_buffer);

  PG_RETURN_CSTRING (barsoi_version_buffer);
}


//PG_FUNCTION_INFO_V1 (pgchem_validate_molfile);
//
//Datum
//pgchem_validate_molfile (PG_FUNCTION_ARGS)
//{
//  bool retval = false;
//  bytea *arg_molfile = PG_GETARG_BYTEA_P (0);
//
//  char *molfile = VARDATA (arg_molfile);
//
//  retval = (strstr (molfile, "V2000") != NULL
//          && strstr (molfile, "M  END") != NULL
//          && strstr (molfile, "\n") != NULL);
//
//  PG_RETURN_BOOL (retval);
//}

//PG_FUNCTION_INFO_V1 (pgchem_is_checkmol_molfile);
//
//Datum
//pgchem_is_checkmol_molfile (PG_FUNCTION_ARGS)
//{
//  bool retval = false;
//  bytea *arg_molfile = PG_GETARG_BYTEA_P (0);
//
//  retval = is_checkmol_molfile (arg_molfile);
//
//  PG_RETURN_BOOL (retval);
//}

/*
* Check if molecule is a nostruct.
*/
PG_FUNCTION_INFO_V1 (pgchem_is_nostruct);

Datum
pgchem_is_nostruct (PG_FUNCTION_ARGS)
{
  bool retval = false;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval =
    (ob_is_nostruct (MFPTR (arg_molecule)) != 0);

  PG_RETURN_BOOL (retval);
}

//PG_FUNCTION_INFO_V1 (pgchem_fg_fingerprint);
//
//Datum
//pgchem_fg_fingerprint (PG_FUNCTION_ARGS)
//{
//  FILE *process;
//  char *fgbuffer;
//  int eno;
//  cow_orker_info *cw_info;
//
//  const int fgbuffer_len = 32;
//
//  bytea *arg_molfile = PG_GETARG_BYTEA_P (0);
//
//  int32 datalen_molfile = VARSIZE (arg_molfile) - VARHDRSZ;
//
//  bytea *retval = (bytea *) palloc (fgbuffer_len + VARHDRSZ);
//
//  VARATT_SIZEP (retval) = fgbuffer_len + VARHDRSZ;
//
//  fgbuffer = (char *) palloc (fgbuffer_len);
//
//#ifdef WIN32
//  cw_info =
//    spawn_cow_orker ("checkmol", "-b", (const char *) VARDATA (arg_molfile),
//                   datalen_molfile, NULL, 0);
//#else
//  cw_info =
//    spawn_cow_orker ("checkmol", "-b", (const char *) VARDATA (arg_molfile),
//                   datalen_molfile);
//#endif
//
//  process = cw_info->process;
//
//  if (process != NULL)
//    {
//
//      fread (fgbuffer, fgbuffer_len, 1, process);
//
//    }
//  else
//    {
//      ereport (ERROR, (errmsg ("checkmol start failed")));
//    }
//
//  eno = cow_orker_cleanup (cw_info);
//
//  if (eno == 127)
//    {
//      ereport (ERROR, (errmsg ("checkmol start failed: 127")));
//    }
//
//
//
//  /* copy the return value to a variable of psql type text */
//  memcpy (VARDATA (retval), fgbuffer, fgbuffer_len);
//
//  //pfree (fgbuffer);
//
//  /* return as psql type by reference */
//  PG_RETURN_BYTEA_P (retval);
//}

/*PG_FUNCTION_INFO_V1 (pgchem_fingerprint2);

Datum
pgchem_fingerprint2 (PG_FUNCTION_ARGS)
{
  bytea *arg_molfile;
  MOLFP2 retval;

  //const int fpbuffer_len = 128;
  //char fpbuffer[FPSIZE];

  arg_molfile = PG_GETARG_BYTEA_P (0);

  /* This is a bitfield fingerprint, so just to be sure we don't get garbage 
  //memset (fpbuffer, 0, FPSIZE);

  retval = (MOLFP2) palloc (FPSIZE2);

  ob_fp2 (VARDATA (arg_molfile), (unsigned int *) retval);

  //VARATT_SIZEP (retval) = fpbuffer_len + VARHDRSZ;

  /* copy the return value to a variable of psql type text 
  //memcpy (retval, fpbuffer, FPSIZE);

  /* return as psql type by reference 
  PG_RETURN_MOLFP2_P (retval);
}*/

/*Datum
pgchem_fingerprint2 (PG_FUNCTION_ARGS)
{
  bytea *arg_molfile;
  bytea *retval;

  const int fpbuffer_len = 128;
  char fpbuffer[fpbuffer_len];

  arg_molfile = PG_GETARG_BYTEA_P (0);

  /* This is a bitfield fingerprint, so just to be sure we don't get garbage 
  memset (fpbuffer, 0, fpbuffer_len);

  ob_fp2 (VARDATA (arg_molfile), (unsigned int *) fpbuffer);

  retval = (bytea *) palloc (fpbuffer_len + VARHDRSZ);

  VARATT_SIZEP (retval) = fpbuffer_len + VARHDRSZ;

  /* copy the return value to a variable of psql type text 
  memcpy (VARDATA (retval), fpbuffer, fpbuffer_len);

  /* return as psql type by reference 
  PG_RETURN_BYTEA_P (retval);
}*/

/*PG_FUNCTION_INFO_V1 (pgchem_tanimoto);

Datum
pgchem_tanimoto (PG_FUNCTION_ARGS)
{
  float8 tanimoto = 0.0;

  MOLFP2 *arg_fp_needle = PG_GETARG_MOLFP2_P (0);
  MOLFP2 *arg_fp_haystack = PG_GETARG_MOLFP2_P (1);

  tanimoto =
    (float8) ob_tanimoto ((unsigned int *) arg_fp_needle,
			  (unsigned int *) arg_fp_haystack);

  // return as psql type 
  PG_RETURN_FLOAT8 (tanimoto);
}*/

/*Datum
pgchem_tanimoto (PG_FUNCTION_ARGS)
{
  float8 tanimoto = 0.0;
  int j, k;

  bytea *arg_fp_needle = PG_GETARG_BYTEA_P (0);
  bytea *arg_fp_haystack = PG_GETARG_BYTEA_P (1);

  j = VARSIZE (arg_fp_needle) - VARHDRSZ;
  k = VARSIZE (arg_fp_haystack) - VARHDRSZ;

  if (j < 1 || k < 1 || j != k)
    {
      ereport (ERROR,
	       (errmsg
		("inconsistent fingerprint sizes read. aborting query.")));
    }
  else
    {
      tanimoto =
	(float8) ob_tanimoto ((unsigned int *) VARDATA (arg_fp_needle),
			      (unsigned int *) VARDATA (arg_fp_haystack));
    }

  // return as psql type 
  PG_RETURN_FLOAT8 (tanimoto);
}*/


/*PG_FUNCTION_INFO_V1 (pgchem_match_substruct_nhrr);

 Datum
pgchem_match_substruct_nhr (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (match_substruct
		  (PG_GETARG_BYTEA_P (0), PG_GETARG_BYTEA_P (1), false, false));
}


PG_FUNCTION_INFO_V1 (pgchem_match_substruct_hnrr);

 Datum
pgchem_match_substruct_hnr (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (match_substruct
		  (PG_GETARG_BYTEA_P (1), PG_GETARG_BYTEA_P (0), false, false));
}

PG_FUNCTION_INFO_V1 (pgchem_match_substruct_nhrs);

 Datum
pgchem_match_substruct_nhrs (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (match_substruct
		  (PG_GETARG_BYTEA_P (0), PG_GETARG_BYTEA_P (1), false, true));
}

PG_FUNCTION_INFO_V1 (pgchem_match_substruct_hnrs);

 Datum
pgchem_match_substruct_hnrs (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (match_substruct
		  (PG_GETARG_BYTEA_P (1), PG_GETARG_BYTEA_P (0), false, true));
}

PG_FUNCTION_INFO_V1 (pgchem_match_substruct_nhsr);

  Datum
pgchem_match_substruct_nhsr (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (match_substruct
		  (PG_GETARG_BYTEA_P (0), PG_GETARG_BYTEA_P (1), true, false));
}

PG_FUNCTION_INFO_V1 (pgchem_match_substruct_hnsr);

 Datum
pgchem_match_substruct_hnsr (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (match_substruct
		  (PG_GETARG_BYTEA_P (1), PG_GETARG_BYTEA_P (0), true, false));
}

PG_FUNCTION_INFO_V1 (pgchem_match_substruct_nhss);
  Datum
pgchem_match_substruct_nhss (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (match_substruct
		  (PG_GETARG_BYTEA_P (0), PG_GETARG_BYTEA_P (1), true, true));
}


PG_FUNCTION_INFO_V1 (pgchem_match_substruct_hnss);

 Datum
pgchem_match_substruct_hnss (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (match_substruct
		  (PG_GETARG_BYTEA_P (1), PG_GETARG_BYTEA_P (0), true, true));
} */

//PG_FUNCTION_INFO_V1 (pgchem_match_substruct);
//
//Datum
//pgchem_match_substruct (PG_FUNCTION_ARGS)
//{
//  PG_RETURN_BOOL (match_substruct
//                (PG_GETARG_BYTEA_P (0), PG_GETARG_BYTEA_P (1),
//                 PG_GETARG_BOOL (2), PG_GETARG_BOOL (3),
//                 PG_GETARG_BOOL (4)));
//}
//
//PG_FUNCTION_INFO_V1 (pgchem_match_exact);
//
//Datum
//pgchem_match_exact (PG_FUNCTION_ARGS)
//{
//  PG_RETURN_BOOL (match_exact
//                (PG_GETARG_BYTEA_P (0), PG_GETARG_BYTEA_P (1),
//                 PG_GETARG_BOOL (2), PG_GETARG_BOOL (3)));
//}

/*PG_FUNCTION_INFO_V1 (pgchem_not_match_exact);

 Datum
pgchem_not_match_exact (PG_FUNCTION_ARGS)
{
 PG_RETURN_BOOL (!match_exact(PG_GETARG_BYTEA_P(0),PG_GETARG_BYTEA_P(1)));
}*/

//PG_FUNCTION_INFO_V1 (pgchem_ms_fingerprint_long);
//
//Datum
//pgchem_ms_fingerprint_long (PG_FUNCTION_ARGS)
//{
//  cow_orker_info *cw_info;
//  FILE *process;
//  char *fpbuffer;
//  int eno;
//  text *retval;
//
//  const int fpbuffer_len = 768;
//
//  bytea *arg_molfile = PG_GETARG_BYTEA_P (0);
//
//  int32 datalen_molfile = VARSIZE (arg_molfile) - VARHDRSZ;
//
//  fpbuffer = (char *) palloc (fpbuffer_len);
//
//#ifdef WIN32
//  cw_info =
//    spawn_cow_orker ("checkmol", "-x", (const char *) VARDATA (arg_molfile),
//                   datalen_molfile, NULL, 0);
//#else
//  cw_info =
//    spawn_cow_orker ("checkmol", "-x", (const char *) VARDATA (arg_molfile),
//                   datalen_molfile);
//#endif
//
//  process = cw_info->process;
//
//  if (process != NULL)
//    {
//      memset (fpbuffer, 0, sizeof (fpbuffer));
//
//      fgets (fpbuffer, fpbuffer_len, process);
//    }
//  else
//    {
//      ereport (ERROR, (errmsg ("checkmol start failed")));
//    }
//
//  eno = cow_orker_cleanup (cw_info);
//
//  retval = (text *) palloc (strlen (fpbuffer) + VARHDRSZ);
//  //memset(retval,0,VARHDRSZ+sizeof(fpbuffer));
//  VARATT_SIZEP (retval) = strlen (fpbuffer) + VARHDRSZ;
//
//  if (eno == 127)
//    {
//      ereport (ERROR, (errmsg ("checkmol start failed: 127")));
//    }
//
//
//
//  /* copy the return value to a variable of psql type text */
//  memcpy (VARDATA (retval), fpbuffer, strlen (fpbuffer));
//
//  //pfree (fpbuffer);
//  /* return as psql type by reference */
//  PG_RETURN_TEXT_P (retval);
//}
//
//PG_FUNCTION_INFO_V1 (pgchem_ms_fingerprint_short);
//
//Datum
//pgchem_ms_fingerprint_short (PG_FUNCTION_ARGS)
//{
//  cow_orker_info *cw_info;
//  FILE *process;
//  char *fpbuffer;
//  int eno;
//  text *retval;
//
//  const int fpbuffer_len = 288;
//
//  bytea *arg_molfile = PG_GETARG_BYTEA_P (0);
//
//  int32 datalen_molfile = VARSIZE (arg_molfile) - VARHDRSZ;
//
//  fpbuffer = (char *) palloc (fpbuffer_len);
//
//#ifdef WIN32
//  cw_info =
//    spawn_cow_orker ("checkmol", "-X", (const char *) VARDATA (arg_molfile),
//                   datalen_molfile, NULL, 0);
//#else
//  cw_info =
//    spawn_cow_orker ("checkmol", "-X", (const char *) VARDATA (arg_molfile),
//                   datalen_molfile);
//#endif
//
//  process = cw_info->process;
//
//  if (process != NULL)
//    {
//      memset (fpbuffer, 0, sizeof (fpbuffer));
//
//      fgets (fpbuffer, fpbuffer_len, process);
//    }
//  else
//    {
//      ereport (ERROR, (errmsg ("checkmol start failed")));
//    }
//
//  eno = cow_orker_cleanup (cw_info);
//
//  retval = (text *) palloc (strlen (fpbuffer) + VARHDRSZ);
//  //memset(retval,0,VARHDRSZ+sizeof(fpbuffer));
//  VARATT_SIZEP (retval) = strlen (fpbuffer) + VARHDRSZ;
//
//  if (eno == 127)
//    {
//      ereport (ERROR, (errmsg ("checkmol start failed: 127")));
//    }
//
//
//
//  /* copy the return value to a variable of psql type text */
//  memcpy (VARDATA (retval), fpbuffer, strlen (fpbuffer));
//
//  //pfree (fpbuffer);
//
//  /* return as psql type by reference */
//  PG_RETURN_TEXT_P (retval);
//}

//PG_FUNCTION_INFO_V1 (pgchem_fgroup_codes);
//
//Datum
//pgchem_fgroup_codes (PG_FUNCTION_ARGS)
//{
//
//  cow_orker_info *cw_info;
//  FILE *process;
//  char *fpbuffer;
//  int eno;
//  text *retval;
//
//  const int fpbuffer_len = 1024;
//
//  bytea *arg_molfile = PG_GETARG_BYTEA_P (0);
//
//  int32 datalen_molfile = VARSIZE (arg_molfile) - VARHDRSZ;
//
//  fpbuffer = (char *) palloc (fpbuffer_len);
//
//#ifdef WIN32
//  cw_info =
//    spawn_cow_orker ("checkmol", "-c", (const char *) VARDATA (arg_molfile),
//                   datalen_molfile, NULL, 0);
//#else
//  cw_info =
//    spawn_cow_orker ("checkmol", "-c", (const char *) VARDATA (arg_molfile),
//                   datalen_molfile);
//#endif
//
//  process = cw_info->process;
//
//  if (process != NULL)
//    {
//      memset (fpbuffer, 0, sizeof (fpbuffer));
//
//      fgets (fpbuffer, fpbuffer_len, process);
//    }
//  else
//    {
//      ereport (ERROR, (errmsg ("checkmol start failed")));
//    }
//
//  eno = cow_orker_cleanup (cw_info);
//
//  retval = (text *) palloc (strlen (fpbuffer) + VARHDRSZ);
//  //memset(retval,0,VARHDRSZ+sizeof(fpbuffer));
//  VARATT_SIZEP (retval) = strlen (fpbuffer) + VARHDRSZ;
//
//  if (eno == 127)
//    {
//      ereport (ERROR, (errmsg ("checkmol start failed: 127")));
//    }
//
//
//
//  /* copy the return value to a variable of psql type text */
//  memcpy (VARDATA (retval), fpbuffer, strlen (fpbuffer));
//
//  //pfree (fpbuffer);
//
//  /* return as psql type by reference */
//  PG_RETURN_TEXT_P (retval);
//}

//PG_FUNCTION_INFO_V1 (pgchem_tweak_molecule);
//
//Datum
//pgchem_tweak_molecule (PG_FUNCTION_ARGS)
//{
//  FILE *process;
//  char *molbuffer;
//  char linebuffer[256];
//  int eno;
//  int molbuffer_size;
//  cow_orker_info *cw_info;
//  bytea *retval;
//  bytea *arg_molfile = PG_GETARG_BYTEA_P (0);
//  bool is_forced = PG_GETARG_BOOL (1);
//
//
//  bool already_processed = FALSE;
//  already_processed = is_checkmol_molfile (arg_molfile);
//
//  if (!is_forced && already_processed)
//    {
//      retval = arg_molfile;
//    }
//  else
//    {
//
//      int32 datalen_molfile = VARSIZE (arg_molfile) - VARHDRSZ;
//
//      molbuffer_size = datalen_molfile * 2;
//
//      molbuffer = (char *) palloc (molbuffer_size);
//
//#ifdef WIN32
//
//      cw_info =
//      spawn_cow_orker ("checkmol", "-m",
//                       (const char *) VARDATA (arg_molfile),
//                       datalen_molfile, NULL, 0);
//#else
//
//      cw_info =
//      spawn_cow_orker ("checkmol", "-m",
//                       (const char *) VARDATA (arg_molfile),
//                       datalen_molfile);
//#endif
//
//      process = cw_info->process;
//
//      if (process != NULL)
//      {
//        memset (molbuffer, 0, molbuffer_size);
//        linebuffer[0] = '\0';
//
//        while (fgets (linebuffer, sizeof (linebuffer), process) != NULL)
//          {
//            strcat (molbuffer, linebuffer);
//            memset (linebuffer, 0, sizeof (linebuffer));
//          }
//      }
//      else
//      {
//        ereport (ERROR, (errmsg ("checkmol start failed")));
//      }
//
//      eno = cow_orker_cleanup (cw_info);
//
//
//      retval = (bytea *) palloc (strlen (molbuffer) + VARHDRSZ);
//      //memset(retval,0,VARHDRSZ+sizeof(fpbuffer));
//      VARATT_SIZEP (retval) = strlen (molbuffer) + VARHDRSZ;
//
//      if (eno == 127)
//      {
//        ereport (ERROR, (errmsg ("checkmol start failed: 127")));
//      }
//
//      /* copy the return value to a variable of psql type bytea */
//      memcpy (VARDATA (retval), molbuffer, strlen (molbuffer));
//
//      //pfree (molbuffer);
//
//    }
//
//
//  /* return as psql type by reference */
//  PG_RETURN_BYTEA_P (retval);
//}

/*
* Convert molecule to SMILES.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_smiles);

Datum
pgchem_molecule_to_smiles (PG_FUNCTION_ARGS)
{
  char *tmpSmiles=NULL;
  text *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  bool omit_iso_and_chiral_markings = PG_GETARG_BOOL (1);
  int len;

  tmpSmiles =
    ob_smiles_to_smiles (SMIPTR (arg_molecule), omit_iso_and_chiral_markings ? 1 : 0);
				
  if(tmpSmiles == NULL) {
              goto smiles_fail;
        }  else if (!strlen(tmpSmiles)) {
              free(tmpSmiles);
              smiles_fail: elog (ERROR, "SMILES generation failed! Offender was :\n %s",MFPTR (arg_molecule));
        }    

  len = strlen (tmpSmiles);

  retval = (text *) palloc (len + VARHDRSZ);
  memset(retval,0x0,len + VARHDRSZ);
  
  SET_VARSIZE (retval,(len + VARHDRSZ));

  strncpy (VARDATA (retval), tmpSmiles, len);

  free (tmpSmiles);

  PG_RETURN_TEXT_P (retval);
}

/*
* Convert molecule to canonical SMILES.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_canonical_smiles);

Datum
pgchem_molecule_to_canonical_smiles (PG_FUNCTION_ARGS)
{
  char *tmpSmiles=NULL;
  text *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  bool omit_iso_and_chiral_markings = PG_GETARG_BOOL (1);
  int len;

  tmpSmiles =
    ob_smiles_to_canonical_smiles (SMIPTR (arg_molecule), omit_iso_and_chiral_markings ? 1 : 0);
				
  if(tmpSmiles == NULL) {
              goto smiles_fail;
        }  else if (!strlen(tmpSmiles)) {
              free(tmpSmiles);
              smiles_fail: elog (ERROR, "Canonical SMILES generation failed! Offender was :\n %s",MFPTR (arg_molecule));
        }    

  len = strlen (tmpSmiles);

  retval = (text *) palloc (len + VARHDRSZ);
  memset(retval,0x0,len + VARHDRSZ);
  
  SET_VARSIZE (retval,(len + VARHDRSZ));

  strncpy (VARDATA (retval), tmpSmiles, len);

  free (tmpSmiles);

  PG_RETURN_TEXT_P (retval);
}

/*
* Convert molecule to InChI.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_inchi);

Datum
pgchem_molecule_to_inchi (PG_FUNCTION_ARGS)
{
  char *tmpInChI = NULL;
  text *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  int len;

  tmpInChI = ob_smiles_to_inchi (SMIPTR (arg_molecule));
  
  if(tmpInChI == NULL) {
              goto inchi_fail;
        }  else if (!strlen(tmpInChI) || strstr (tmpInChI, "InChI=") == NULL) {
              free(tmpInChI);
              inchi_fail: elog (ERROR, "InChI generation failed! Offender was :\n %s",MFPTR (arg_molecule));
        }    

  len = strlen (tmpInChI);

  retval = (text *) palloc (len + VARHDRSZ);
  memset(retval,0x0,len + VARHDRSZ);
  
  SET_VARSIZE (retval,(len + VARHDRSZ));

  strncpy (VARDATA (retval), tmpInChI, len);

  free (tmpInChI);

  PG_RETURN_TEXT_P (retval);
}

/*
* Convert molecule to V3000.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_V3000);

Datum
pgchem_molecule_to_V3000 (PG_FUNCTION_ARGS)
{
  char *tmpV3000=NULL;
  text *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  int len;

  tmpV3000 = ob_mol_to_V3000 (MFPTR (arg_molecule));
  
   if(tmpV3000 == NULL) {
              goto v3000_fail;
        }  else if (!strlen(tmpV3000) || strstr (tmpV3000, "V3000") == NULL) {
              free(tmpV3000);
              v3000_fail: elog (ERROR, "V3000 generation failed! Offender was :\n %s",MFPTR (arg_molecule));
        }      

  len = strlen (tmpV3000);

  retval = (text *) palloc (len + VARHDRSZ);
  memset(retval,0x0,len+VARHDRSZ);
   
  SET_VARSIZE (retval,(len + VARHDRSZ));

 strncpy (VARDATA (retval), tmpV3000, len);

  free (tmpV3000);

  PG_RETURN_TEXT_P (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_molecule_to_molfile);

Datum
pgchem_molecule_to_molfile (PG_FUNCTION_ARGS)
{
  MOLECULE *molecule = PG_GETARG_MOLECULE_P (0);
  text *retval;
  int len = (molecule->sizemf)-1;

  retval = (text *) palloc (len+VARHDRSZ);
  memset(retval,0x0,len+VARHDRSZ);

  strncpy (VARDATA (retval), MFPTR(molecule), len);
  
  SET_VARSIZE (retval,(len + VARHDRSZ));
  
  //printf(result);

  PG_RETURN_TEXT_P (retval);
}


/*****************************************************************
 * This function Copyright (c) 2005
 * by Bayer Business Services GmbH
 *****************************************************************/

/*
* Calculate molweight of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_molweight);

Datum
pgchem_molweight (PG_FUNCTION_ARGS)
{
  float8 retval = 0.0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_molweight (SMIPTR (arg_molecule));

  PG_RETURN_FLOAT8 (retval);
}

/*
* Calculate Hill formula of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_hillformula);

Datum
pgchem_hillformula (PG_FUNCTION_ARGS)
{
  char *tmpFormula;
  text *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  int len;

  tmpFormula = ob_hillformula (SMIPTR (arg_molecule));

  len = strlen (tmpFormula);

  retval = (text *) palloc (len + VARHDRSZ);
  memset(retval,0x0,len + VARHDRSZ);
  
  SET_VARSIZE (retval,(len + VARHDRSZ));

  strncpy (VARDATA (retval), tmpFormula, len);

  free (tmpFormula);

  PG_RETURN_TEXT_P (retval);
}

/*
* Calculate exact mass of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_exactmass);

Datum
pgchem_exactmass (PG_FUNCTION_ARGS)
{
  float8 retval = 0.0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_exactmass (SMIPTR (arg_molecule));

  PG_RETURN_FLOAT8 (retval);
}

/*
* Calculate total charge of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_total_charge);

Datum
pgchem_total_charge (PG_FUNCTION_ARGS)
{
  int32 retval = 0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_total_charge (SMIPTR (arg_molecule));

  PG_RETURN_INT32 (retval);
}

/*
* Calculate number of atoms of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_num_atoms);

Datum
pgchem_num_atoms (PG_FUNCTION_ARGS)
{
  uint32 retval = 0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_num_atoms (SMIPTR (arg_molecule));

  PG_RETURN_UINT32 (retval);
}

/*
* Calculate number of heavy atoms of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_num_heavy_atoms);

Datum
pgchem_num_heavy_atoms (PG_FUNCTION_ARGS)
{
  uint32 retval = 0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_num_heavy_atoms (SMIPTR (arg_molecule));

  PG_RETURN_UINT32 (retval);
}

/*
* Calculate number of bonds of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_num_bonds);

Datum
pgchem_num_bonds (PG_FUNCTION_ARGS)
{
  uint32 retval = 0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_num_bonds (SMIPTR (arg_molecule));

  PG_RETURN_UINT32 (retval);
}

/*
* Calculate number of rotatable bonds of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_num_rotatable_bonds);

Datum
pgchem_num_rotatable_bonds (PG_FUNCTION_ARGS)
{
  uint32 retval = 0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval =
    ob_num_rotatable_bonds (SMIPTR (arg_molecule));

  PG_RETURN_UINT32 (retval);
}

/*
* Determine chirality of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_is_chiral);

Datum
pgchem_is_chiral (PG_FUNCTION_ARGS)
{
  bool retval = false;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval =
    (ob_is_chiral (MFPTR (arg_molecule)) ==
     1) ? true : false;

  PG_RETURN_BOOL (retval);
}

/*
* Molecule has 2D coordinates.
*/
PG_FUNCTION_INFO_V1 (pgchem_2D);

Datum
pgchem_2D (PG_FUNCTION_ARGS)
{
  bool retval = false;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval =
    (ob_2D (MFPTR (arg_molecule)) == 1) ? true : false;

  PG_RETURN_BOOL (retval);
}

/*
* Molecule has 3D coordinates.
*/
PG_FUNCTION_INFO_V1 (pgchem_3D);

Datum
pgchem_3D (PG_FUNCTION_ARGS)
{
  bool retval = false;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval =
    (ob_3D (MFPTR (arg_molecule)) == 1) ? true : false;

  PG_RETURN_BOOL (retval);
}



//inline static bool
//match_exact_a (bytea * molfile_needle, bytea * molfile_haystack,
//             bool strict_typing, bool e_z_global, bool r_s_global,
//             bool strict_chg, bool strict_iso, bool strict_rad)
//{
//
//  bool retval = false;
//
//  //char molfile[MAX_MOL_SIZE + 1];
//
//  int32 size_needle;
//  int32 size_haystack;
//
//  size_needle = VARSIZE (molfile_needle) - VARHDRSZ;
//  size_haystack = VARSIZE (molfile_haystack) - VARHDRSZ;
//
//  if (size_needle > MAX_MOL_SIZE || size_haystack > MAX_MOL_SIZE)
//    ereport (ERROR,
//           (errmsg ("impossibly large molecule read. aborting query.")));
//
//  molfile_buf[0] = '\0';
//
//  strncat (molfile_buf, VARDATA (molfile_needle), size_needle);
//
//  xm_set_ring_perception_algorithm (RPA_SAR);
//  xm_set_strict_typing (strict_typing ? FEATURE_ON : FEATURE_OFF);
//  mm_set_r_s_check (r_s_global ? FEATURE_ON : FEATURE_OFF);
//  mm_set_e_z_check (e_z_global ? FEATURE_ON : FEATURE_OFF);
//  mm_set_chg_check (strict_chg ? FEATURE_ON : FEATURE_OFF);
//  mm_set_iso_check (strict_iso ? FEATURE_ON : FEATURE_OFF);
//  mm_set_rad_check (strict_rad ? FEATURE_ON : FEATURE_OFF);
//  mm_set_exact_match (FEATURE_ON);
//
//  if (strcmp (molfile_buf, current_molfile) != 0)
//    {
//      strncpy (current_molfile, molfile_buf, MAX_MOL_SIZE);
//      mm_set_mol (molfile_buf);
//      mm_set_current_mol_as_query ();
//    }
//
//  molfile_buf[0] = '\0';
//  strncat (molfile_buf, VARDATA (molfile_haystack), size_haystack);
//  mm_set_mol (molfile_buf);
//
//  retval = (mm_match () != 0);
//
//  return retval;
//}
//
//PG_FUNCTION_INFO_V1 (pgchem_match_exact_a);
//
//Datum
//pgchem_match_exact_a (PG_FUNCTION_ARGS)
//{
//  PG_RETURN_BOOL (match_exact_a
//                (PG_GETARG_BYTEA_P (0), PG_GETARG_BYTEA_P (1),
//                 PG_GETARG_BOOL (2), PG_GETARG_BOOL (3),
//                 PG_GETARG_BOOL (4), PG_GETARG_BOOL (5), PG_GETARG_BOOL (6),
//                 PG_GETARG_BOOL (7)));
//}
//
//inline static bool
//match_substruct_a (bytea * molfile_needle, bytea * molfile_haystack,
//                 bool strict_typing, bool e_z_global, bool r_s_global,
//                 bool strict_chg, bool strict_iso, bool strict_rad)
//{
//
//  bool retval = false;
//
//  //char molfile[MAX_MOL_SIZE + 1];
//
//  int32 size_needle;
//  int32 size_haystack;
//
//  size_needle = VARSIZE (molfile_needle) - VARHDRSZ;
//  size_haystack = VARSIZE (molfile_haystack) - VARHDRSZ;
//
//  if (size_needle > MAX_MOL_SIZE || size_haystack > MAX_MOL_SIZE)
//    ereport (ERROR,
//           (errmsg ("impossibly large molecule read. aborting query.")));
//
//  molfile_buf[0] = '\0';
//
//  strncat (molfile_buf, VARDATA (molfile_needle), size_needle);
//
//  xm_set_ring_perception_algorithm (RPA_SAR);
//  xm_set_strict_typing (strict_typing ? FEATURE_ON : FEATURE_OFF);
//  mm_set_r_s_check (r_s_global ? FEATURE_ON : FEATURE_OFF);
//  mm_set_e_z_check (e_z_global ? FEATURE_ON : FEATURE_OFF);
//  mm_set_chg_check (strict_chg ? FEATURE_ON : FEATURE_OFF);
//  mm_set_iso_check (strict_iso ? FEATURE_ON : FEATURE_OFF);
//  mm_set_rad_check (strict_rad ? FEATURE_ON : FEATURE_OFF);
//  mm_set_exact_match (FEATURE_OFF);
//
//  if (strcmp (molfile_buf, current_molfile) != 0)
//    {
//      strncpy (current_molfile, molfile_buf, MAX_MOL_SIZE);
//      mm_set_mol (molfile_buf);
//      mm_set_current_mol_as_query ();
//    }
//
//  molfile_buf[0] = '\0';
//  strncat (molfile_buf, VARDATA (molfile_haystack), size_haystack);
//  mm_set_mol (molfile_buf);
//
//  retval = (mm_match () != 0);
//
//  return retval;
//}
//
//PG_FUNCTION_INFO_V1 (pgchem_match_substruct_a);
//
//Datum
//pgchem_match_substruct_a (PG_FUNCTION_ARGS)
//{
//  PG_RETURN_BOOL (match_substruct_a
//                (PG_GETARG_BYTEA_P (0), PG_GETARG_BYTEA_P (1),
//                 PG_GETARG_BOOL (2), PG_GETARG_BOOL (3),
//                 PG_GETARG_BOOL (4), PG_GETARG_BOOL (5), PG_GETARG_BOOL (6),
//                 PG_GETARG_BOOL (7)));
//}

/*
* Calculate molstatistics fingerprint of molecule.
* Long version: n_atoms:6;n_bonds:6...
*/
PG_FUNCTION_INFO_V1 (pgchem_ms_fingerprint_long_a);

Datum
pgchem_ms_fingerprint_long_a (PG_FUNCTION_ARGS)
{
  MOLECULE *arg_molecule;
  text *retval;
  int len;

  bool strict_chg = PG_GETARG_BOOL (1);
  bool strict_iso = PG_GETARG_BOOL (2);
  bool strict_rad = PG_GETARG_BOOL (3);
  //char molfile[MAX_MOL_SIZE + 1];

  const int fpbuffer_len = 1024;
  char fpbuffer[fpbuffer_len];

  arg_molecule = PG_GETARG_MOLECULE_P (0);

  mm_set_chg_check (strict_chg ? FEATURE_ON : FEATURE_OFF);
  mm_set_iso_check (strict_iso ? FEATURE_ON : FEATURE_OFF);
  mm_set_rad_check (strict_rad ? FEATURE_ON : FEATURE_OFF);
  xm_set_ring_perception_algorithm (RPA_SAR);
  cm_set_mol (MFPTR (arg_molecule), FEATURE_ON);
  cm_molstat (fpbuffer);

  len = strlen (fpbuffer);

  retval = (text *) palloc (len + VARHDRSZ);

  SET_VARSIZE (retval,(len + VARHDRSZ));

  /* copy the return value to a variable of psql type text */
  memcpy (VARDATA (retval), fpbuffer, len);

  /* return as psql type by reference */
  PG_RETURN_TEXT_P (retval);
}

/*
* Calculate molstatistics fingerprint of molecule.
* Short version: 6,6,1,0,0,0,0,6,6,0,0...
*/
PG_FUNCTION_INFO_V1 (pgchem_ms_fingerprint_short_a);

Datum
pgchem_ms_fingerprint_short_a (PG_FUNCTION_ARGS)
{
  MOLECULE *arg_molecule;
  text *retval;
  int len;

  //char molfile[MAX_MOL_SIZE + 1];

  char fpbuffer[1024];

  arg_molecule = PG_GETARG_MOLECULE_P (0);
  
  memset(fpbuffer,0x0,1024*sizeof(char));

  xm_set_ring_perception_algorithm (RPA_SAR);
  cm_set_mol (MFPTR (arg_molecule), FEATURE_ON);
  cm_molstat_X (fpbuffer);

  len = strlen (fpbuffer);

  retval = (text *) palloc (len + VARHDRSZ);

  SET_VARSIZE (retval,(len + VARHDRSZ));

  /* copy the return value to a variable of psql type text */
  memcpy (VARDATA (retval), fpbuffer, len);

  /* return as psql type by reference */
  PG_RETURN_TEXT_P (retval);
}
/*
* Calculate functional group codes of molecule.
* 0000A000;...
*/
PG_FUNCTION_INFO_V1 (pgchem_fgroup_codes_a);

Datum
pgchem_fgroup_codes_a (PG_FUNCTION_ARGS)
{
  MOLECULE *arg_molecule;
  text *retval;
  int len;

  //char molfile[MAX_MOL_SIZE + 1];

  const int fpbuffer_len = 1024;
  char fpbuffer[fpbuffer_len];

  arg_molecule = PG_GETARG_MOLECULE_P (0);

  xm_set_ring_perception_algorithm (RPA_SAR);
  cm_set_mol (MFPTR (arg_molecule), FEATURE_ON);
  cm_fg_codes (fpbuffer);

  len = strlen (fpbuffer);

  retval = (text *) palloc (len + VARHDRSZ);

  SET_VARSIZE (retval,(len + VARHDRSZ));

  /* copy the return value to a variable of psql type text */
  memcpy (VARDATA (retval), fpbuffer, len);

  /* return as psql type by reference */
  PG_RETURN_TEXT_P (retval);
}

//PG_FUNCTION_INFO_V1 (pgchem_tweak_molecule_a);
//
//Datum
//pgchem_tweak_molecule_a (PG_FUNCTION_ARGS)
//{
//
//  char *molbuffer;
//  int molbuffer_size;
//  bytea *retval;
//  int32 size_needle;
//  //char molfile[MAX_MOL_SIZE + 1];
//
//  bytea *arg_molfile;
//  bool is_forced;
//
//  arg_molfile = PG_GETARG_BYTEA_P (0);
//  is_forced = PG_GETARG_BOOL (1);
//
//  if (!is_forced && is_checkmol_molfile (arg_molfile))
//    {
//      retval = arg_molfile;
//    }
//  else
//    {
//
//      size_needle = VARSIZE (arg_molfile) - VARHDRSZ;
//
//      if (size_needle > MAX_MOL_SIZE)
//      ereport (ERROR,
//               (errmsg
//                ("impossibly large molecule read. aborting query.")));
//
//
//      molfile_buf[0] = '\0';
//
//      strncat (molfile_buf, VARDATA (arg_molfile), size_needle);
//
//      molbuffer_size = size_needle * 2;
//
//      molbuffer = (char *) palloc (molbuffer_size);
//
//      xm_set_ring_perception_algorithm (RPA_SAR);
//      cm_set_mol (molfile_buf, FEATURE_OFF);
//      cm_tweak_molfile (molbuffer);
//
//      retval = (bytea *) palloc (strlen (molbuffer) + VARHDRSZ);
//      //memset(retval,0,VARHDRSZ+sizeof(fpbuffer));
//      VARATT_SIZEP (retval) = strlen (molbuffer) + VARHDRSZ;
//
//      /* copy the return value to a variable of psql type bytea */
//      memcpy (VARDATA (retval), molbuffer, strlen (molbuffer));
//
//      //pfree (ne);
//      //pfree (molbuffer);
//
//    }
//
//
//  /* return as psql type by reference */
//  PG_RETURN_BYTEA_P (retval);
//}

/*****************************************************************
 * This function Copyright (c) 2005
 * by Bayer Business Services GmbH
 *****************************************************************/

/*PG_FUNCTION_INFO_V1 (pgchem_check_small_fragment);

Datum
pgchem_check_small_fragment (PG_FUNCTION_ARGS)
{
  bool retval = FALSE;
  int size_needle;
  int size_small_fragment;
  bytea *molfile_needle;
  text *which_small_fragment;
  char small_fragment_name[MAX_SMALL_FRAGMENT_NAME_SIZE + 1];
  bool exact = FALSE;
  //char molfile[MAX_MOL_SIZE + 1];

  molfile_needle = PG_GETARG_BYTEA_P (0);
  which_small_fragment = PG_GETARG_TEXT_P (1);
  exact = PG_GETARG_BOOL (2);

  size_needle = VARSIZE (molfile_needle) - VARHDRSZ;
  size_small_fragment = VARSIZE (which_small_fragment) - VARHDRSZ;

  if (size_needle > MAX_MOL_SIZE)
    ereport (ERROR,
	     (errmsg ("impossibly large molecule read. aborting query.")));


  molfile_buf[0] = '\0';
  small_fragment_name[0] = '\0';

  strncat (molfile_buf, VARDATA (molfile_needle), size_needle);
  strncat (small_fragment_name, VARDATA (which_small_fragment),
	   (size_small_fragment <
	    MAX_SMALL_FRAGMENT_NAME_SIZE) ? size_small_fragment :
	   MAX_SMALL_FRAGMENT_NAME_SIZE);

  xm_set_ring_perception_algorithm (RPA_SAR);

  if (strcasecmp (small_fragment_name, "cyclopentadiene") == 0)
    mm_set_mol (sf_cyclopentadiene);
  else if (strcasecmp (small_fragment_name, "benzene") == 0)
    mm_set_mol (sf_benzene);
  else
    ereport (ERROR,
	     (errmsg
	      ("Small fragment '%s' is unknown.", small_fragment_name)));

  mm_set_current_mol_as_query ();
  mm_set_mol (molfile_buf);

  if (exact)
    retval = (mm_match_exact (0, 0) != 0);
  else
    retval = (mm_match_substruct (0, 0, 0) != 0);

  //pfree(hs);

  PG_RETURN_BOOL (retval);
}*/

/* PG_FUNCTION_INFO_V1 (pgchem_fastmatch_a);

Datum
pgchem_fastmatch_a (PG_FUNCTION_ARGS)
{

  bool retval = false;
  bytea *molfile_haystack;
  int32 size_haystack;

  molfile_haystack = PG_GETARG_BYTEA_P (0);

  size_haystack = VARSIZE (molfile_haystack) - VARHDRSZ;

  if (size_haystack > MAX_MOL_SIZE)
    ereport (ERROR,
	     (errmsg ("impossibly large molecule read. aborting query.")));

  molfile_buf[0] = '\0';
  strncat (molfile_buf, VARDATA (molfile_haystack), size_haystack);
  mm_set_mol (molfile_buf);

  retval = (mm_match () != 0);

  PG_RETURN_BOOL (retval);
}

inline static bool
register_query_4_fastmatch_a (bytea * arg_molfile, bool exact,
			      bool strict_typing, bool e_z_global,
			      bool r_s_global, bool strict_chg,
			      bool strict_iso, bool strict_rad)
{
  int32 size_needle;

  size_needle = VARSIZE (arg_molfile) - VARHDRSZ;

  if (size_needle > MAX_MOL_SIZE)
    ereport (ERROR,
	     (errmsg ("impossibly large molecule read. aborting query.")));

  molfile_buf[0] = '\0';

  strncat (molfile_buf, VARDATA (arg_molfile), size_needle);

  xm_set_ring_perception_algorithm (RPA_SAR);
  xm_set_strict_typing (strict_typing ? FEATURE_ON : FEATURE_OFF);
  mm_set_r_s_check (r_s_global ? FEATURE_ON : FEATURE_OFF);
  mm_set_e_z_check (e_z_global ? FEATURE_ON : FEATURE_OFF);
  mm_set_chg_check (strict_chg ? FEATURE_ON : FEATURE_OFF);
  mm_set_iso_check (strict_iso ? FEATURE_ON : FEATURE_OFF);
  mm_set_rad_check (strict_rad ? FEATURE_ON : FEATURE_OFF);
  mm_set_exact_match (exact ? FEATURE_ON : FEATURE_OFF);
  mm_set_mol (molfile_buf);
  mm_set_current_mol_as_query ();

  PG_RETURN_BOOL (true);
}

PG_FUNCTION_INFO_V1 (pgchem_fastmatch_set_exact_query_a);

Datum
pgchem_fastmatch_set_exact_query_a (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (register_query_4_fastmatch_a
		  (PG_GETARG_BYTEA_P (0), TRUE,
		   PG_GETARG_BOOL (1), PG_GETARG_BOOL (2),
		   PG_GETARG_BOOL (3), PG_GETARG_BOOL (4), PG_GETARG_BOOL (5),
		   PG_GETARG_BOOL (6)));
}

PG_FUNCTION_INFO_V1 (pgchem_fastmatch_set_substructure_query_a);

Datum
pgchem_fastmatch_set_substructure_query_a (PG_FUNCTION_ARGS)
{
  PG_RETURN_BOOL (register_query_4_fastmatch_a
		  (PG_GETARG_BYTEA_P (0), FALSE,
		   PG_GETARG_BOOL (1), PG_GETARG_BOOL (2),
		   PG_GETARG_BOOL (3), PG_GETARG_BOOL (4), PG_GETARG_BOOL (5),
		   PG_GETARG_BOOL (6)));
} */

//PG_FUNCTION_INFO_V1 (pgchem_match_substruct_ob_a);
//
//Datum
//pgchem_match_substruct_ob_a (PG_FUNCTION_ARGS)
//{
//  bool retval = false;
//  text *arg_smarts = PG_GETARG_TEXT_P (0);
//  bytea *arg_molfile = PG_GETARG_BYTEA_P (1);
//  int smarts_string_len = VARSIZE (arg_smarts) - VARHDRSZ;
//  char *tmpSMARTS = (char *) palloc (smarts_string_len + 1);
//  tmpSMARTS[0] = '\0';
//
//  strncat (tmpSMARTS, VARDATA (arg_smarts), smarts_string_len);
//
//  retval = (ob_SSS_SMARTS (tmpSMARTS, VARDATA (arg_molfile)) != 0);
//
//  PG_RETURN_BOOL (retval);
//}
//


/*PG_FUNCTION_INFO_V1 (pgchem_VF2);

Datum
pgchem_VF2 (PG_FUNCTION_ARGS)
{
  bool retval = false;
  MOLECULE *arg_molecule_sarg = PG_GETARG_MOLECULE_P (0);
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (1);
  bool exact = PG_GETARG_BOOL (2);
  int match = 0;
  register int i = FPSIZE2;
  register unsigned int *p = arg_molecule->fp;
  register unsigned int *q = arg_molecule_sarg->fp;
  
  if(exact) {
      while(i--) {
            if (*p != *q) {PG_RETURN_BOOL (false);}
            p++;
            q++;
        }  

         match = ob_ESS_VF2_bin (ANCPTR(arg_molecule_sarg), ANCPTR(arg_molecule));
  } else {    
   while(i--) {
      if ((*p & *q) != *q) {PG_RETURN_BOOL (false);}
      p++;
      q++;
      }

       match = ob_SSS_VF2_bin (ANCPTR(arg_molecule_sarg), ANCPTR(arg_molecule));
  }    

  retval = (match != 0);

  PG_RETURN_BOOL (retval);
}*/

/*
* Match molecule against SMARTS pattern, as requested by Mr. Ertl from Novartis. :-)
*/
PG_FUNCTION_INFO_V1 (pgchem_smartsfilter);

Datum
pgchem_smartsfilter (PG_FUNCTION_ARGS)
{
  bool retval = false;
  text *arg_smarts = PG_GETARG_TEXT_P (0);
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (1);
  int smarts_string_len = VARSIZE (arg_smarts) - VARHDRSZ;
  int match = 0;
  char *tmpSMARTS = (char *) palloc (smarts_string_len + 1);
  tmpSMARTS[0] = '\0';

  strncat (tmpSMARTS, VARDATA (arg_smarts), smarts_string_len);
  
  if(strstr(tmpSMARTS,"@")!=NULL || strstr(tmpSMARTS,"/")!=NULL || strstr(tmpSMARTS,"\\")!=NULL)
    match = ob_SSS_SMARTS_native (tmpSMARTS, SMIPTR(arg_molecule));
    else
    match = ob_SSS_SMARTS_native_bin (tmpSMARTS, ANCPTR(arg_molecule));
  
  if(match<0) elog (ERROR, "Invalid SMARTS pattern: %s",tmpSMARTS);

  retval = (match != 0);

  PG_RETURN_BOOL (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_smartsfilter_count);

Datum
pgchem_smartsfilter_count (PG_FUNCTION_ARGS)
{
  int retval = 0;
  text *arg_smarts = PG_GETARG_TEXT_P (0);
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (1);
  int smarts_string_len = VARSIZE (arg_smarts) - VARHDRSZ;
  char *tmpSMARTS = (char *) palloc (smarts_string_len + 1);
  tmpSMARTS[0] = '\0';

  strncat (tmpSMARTS, VARDATA (arg_smarts), smarts_string_len);

if(strstr(tmpSMARTS,"@")!=NULL || strstr(tmpSMARTS,"/")!=NULL || strstr(tmpSMARTS,"\\")!=NULL)
    retval = ob_SSS_SMARTS_native_count (tmpSMARTS, SMIPTR(arg_molecule));
else
  retval = ob_SSS_SMARTS_native_count_bin (tmpSMARTS, ANCPTR(arg_molecule));
  
  if(retval<0) elog (ERROR, "Invalid SMARTS pattern: %s",tmpSMARTS);

  PG_RETURN_INT32 (retval);
}

/*
* Predict polar surface area (PSA) of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_PSA);

Datum
pgchem_PSA (PG_FUNCTION_ARGS)
{
  float8 retval = 0.0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_PSA (SMIPTR (arg_molecule));

  PG_RETURN_FLOAT8 (retval);
}

/*
* Predict molar refractivity (MR) of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_MR);

Datum
pgchem_MR (PG_FUNCTION_ARGS)
{
  float8 retval = 0.0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_MR (SMIPTR (arg_molecule));

  PG_RETURN_FLOAT8 (retval);
}

/*
* Predict logP of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_logP);

Datum
pgchem_logP (PG_FUNCTION_ARGS)
{
  float8 retval = 0.0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_logP (SMIPTR (arg_molecule));

  PG_RETURN_FLOAT8 (retval);
}

//#ifndef WIN32
//PG_FUNCTION_INFO_V1 (pgchem_lower_connection_priority);
//
//Datum
//pgchem_lower_connection_priority (PG_FUNCTION_ARGS)
//{
//  int32 priority;
//  int result;
//
//  priority = PG_GETARG_INT32 (0);
//
//  if (priority < 0)
//    ereport (ERROR,
//           (errmsg ("cannot raise connection priority, only lower it")));
//
//  result = setpriority (PRIO_PROCESS, 0, priority);
//
//  if (result == 0)
//    {
//      PG_RETURN_BOOL (1);     /* true */
//    }
//  else
//    {
//      PG_RETURN_BOOL (0);     /* false */
//    }
//}
//#endif

/*
* Calculate md5 molhash of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_inchikey);

Datum
pgchem_molecule_to_inchikey (PG_FUNCTION_ARGS)
{
  text *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = (text *) palloc (INCHIKEYSZ + VARHDRSZ);
  
  memset(retval,0x0,INCHIKEYSZ + VARHDRSZ);

  memcpy (VARDATA (retval), arg_molecule->inchikey, INCHIKEYSZ);

  SET_VARSIZE (retval,(INCHIKEYSZ + VARHDRSZ));

  PG_RETURN_TEXT_P (retval);
}

/*
* Get binary fingerprint as bit varying
*/
PG_FUNCTION_INFO_V1 (pgchem_fp_out);

Datum
pgchem_fp_out (PG_FUNCTION_ARGS)
{
  VarBit *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  int len = FPSIZE*sizeof(uint32);

  retval = (VarBit *) palloc (len + VARBITHDRSZ);
  
  memcpy(VARBITS(retval),arg_molecule->fp,len);
  
  VARBITLEN(retval) = len*8; //8 bit chars

  SET_VARSIZE(retval,(len+VARBITHDRSZ));

  PG_RETURN_VARBIT_P (retval);
}

/*
* Get MACCS binary fingerprint as bit varying
*/
PG_FUNCTION_INFO_V1 (pgchem_fp_MACCS);

Datum
pgchem_fp_MACCS (PG_FUNCTION_ARGS)
{
  VarBit *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  int len = FPSIZE_MACCS*sizeof(uint32);
  uint32 *tmp_maccs;
  
  retval = (VarBit *) palloc (len + VARBITHDRSZ);
  
  tmp_maccs = (uint32*) palloc(len);
  
  ob_fp_MACCS_bin(ANCPTR(arg_molecule), tmp_maccs);
  
  memcpy(VARBITS(retval),tmp_maccs,len);
  
  VARBITLEN(retval) = len*8; //8 bit chars

  SET_VARSIZE(retval,(len+VARBITHDRSZ));

  PG_RETURN_VARBIT_P (retval);
}

/*
* Get reaction binary fingerprint as bit varying
*/
PG_FUNCTION_INFO_V1 (pgchem_r_fp_out);

Datum
pgchem_r_fp_out (PG_FUNCTION_ARGS)
{
  VarBit *retval;
  REACTION *arg_reaction = PG_GETARG_REACTION_P (0);
  int len = 2*FPSIZE2*sizeof(uint32);

  retval = (VarBit *) palloc (len + VARBITHDRSZ);
  
  memcpy(VARBITS(retval),arg_reaction->fp,len);
  
  VARBITLEN(retval) = len*8; //8 bit chars

  SET_VARSIZE(retval,(len+VARBITHDRSZ));

  PG_RETURN_VARBIT_P (retval);
}
/*
* Molecule is Benzene.
*/
/*PG_FUNCTION_INFO_V1 (pgchem_isbz);

Datum
pgchem_isbz (PG_FUNCTION_ARGS)
{
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  if (arg_molecule->isbz == true)
    PG_RETURN_BOOL (true);

  PG_RETURN_BOOL (false);
}*/

/*
* Molecule contains Benzene.
*/
/*PG_FUNCTION_INFO_V1 (pgchem_hasbz);

Datum
pgchem_hasbz (PG_FUNCTION_ARGS)
{
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  if (arg_molecule->nobz == false)
    PG_RETURN_BOOL (true);

  PG_RETURN_BOOL (false);
}*/

PG_FUNCTION_INFO_V1 (pgchem_num_H_donors);

Datum
pgchem_num_H_donors (PG_FUNCTION_ARGS)
{
  int32 retval = 0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_num_H_donors (SMIPTR (arg_molecule));

  PG_RETURN_UINT32 (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_num_H_acceptors);

Datum
pgchem_num_H_acceptors (PG_FUNCTION_ARGS)
{
  int32 retval = 0;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

  retval = ob_num_H_acceptors (SMIPTR (arg_molecule));

  PG_RETURN_UINT32 (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_mutate_fp);

Datum
pgchem_mutate_fp (PG_FUNCTION_ARGS)
{
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  uint32 *offset = arg_molecule->fp+OFFSET;
  
  ob_fp3_bin(ANCPTR(arg_molecule), offset);
  
  PG_RETURN_MOLECULE_P (arg_molecule);
}

PG_FUNCTION_INFO_V1 (pgchem_blank_fp);

Datum
pgchem_blank_fp (PG_FUNCTION_ARGS)
{
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  uint32 *offset = arg_molecule->fp+OFFSET;
  
  memset(offset,0x0,FPSIZE3*(sizeof(uint32)));
  
  PG_RETURN_MOLECULE_P (arg_molecule);
}

//PG_FUNCTION_INFO_V1 (pgchem_fp6);

/*Datum
pgchem_fp6 (PG_FUNCTION_ARGS)
{
  VarBit *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  int len = FPSIZE3*4*sizeof(unsigned int);
  char buffer[len];

  retval = (VarBit *) palloc (len + VARBITHDRSZ);
  
  ob_fp6(MFPTR(arg_molecule),(unsigned int *)&buffer);
  
  memcpy(VARBITS(retval),buffer,len);
  
  VARBITLEN(retval) = len*8;

  SET_VARSIZE(retval,(len+VARBITHDRSZ));

  PG_RETURN_VARBIT_P (retval);
}*/

PG_FUNCTION_INFO_V1(pgchem_nbits_set);
Datum
pgchem_nbits_set(PG_FUNCTION_ARGS)
{
/* how many bits are set in a bitstring? */

         VarBit     *a = PG_GETARG_VARBIT_P(0);
         unsigned char *ap = VARBITS(a);
         
       /*  int n=0;
         int i;
         unsigned char aval;
         for (i=0; i < VARBITBYTES(a); ++i) {
                 aval = *ap; ++ap;
                 if (aval == 0) continue;
                 if (aval & 1) ++n;
                 if (aval & 2) ++n;
                 if (aval & 4) ++n;
                 if (aval & 8) ++n;
                 if (aval & 16) ++n;
                 if (aval & 32) ++n;
                 if (aval & 64) ++n;
                 if (aval & 128) ++n;
         }*/
         
         PG_RETURN_INT32(ob_popcount(ap,VARBITBYTES(a)));
} 

PG_FUNCTION_INFO_V1 (pgchem_disconnected);

Datum
pgchem_disconnected (PG_FUNCTION_ARGS)
{
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  
  PG_RETURN_BOOL (arg_molecule->disconnected != 0);
}

PG_FUNCTION_INFO_V1 (pgchem_r_num_reactants);

Datum
pgchem_r_num_reactants(PG_FUNCTION_ARGS){
        REACTION *arg_reaction = PG_GETARG_REACTION_P (0);
       
        PG_RETURN_INT32(arg_reaction->num_reactants);
}

PG_FUNCTION_INFO_V1 (pgchem_r_num_products);

Datum
pgchem_r_num_products(PG_FUNCTION_ARGS){
        REACTION *arg_reaction = PG_GETARG_REACTION_P (0);
       
        PG_RETURN_INT32(arg_reaction->num_products);
}

PG_FUNCTION_INFO_V1 (pgchem_r_molecule_at);

Datum
pgchem_r_molecule_at(PG_FUNCTION_ARGS){
        REACTION *arg_reaction = PG_GETARG_REACTION_P (0);
        int4 arg_position = PG_GETARG_INT32 (1);
        MOLECULE *retval;
        char *offset = MOLARRAYPTR(arg_reaction);
        int i,len;
        
        if(arg_position < 1 || arg_position > arg_reaction->num_products+arg_reaction->num_reactants) elog (ERROR, "Molecule index out of bounds: %d",arg_position);
        
        for(i=1;i<arg_position;i++) {
            len = VARSIZE((MOLECULE*)offset)*sizeof(char);
            offset+=len;
        }
        
        len = VARSIZE((MOLECULE*)offset)*sizeof(char);
        
        retval = (MOLECULE*) palloc(len);
        
        memset(retval,0x0,len);
        
        memcpy(retval,(MOLECULE*)offset,len);                
       
        PG_RETURN_MOLECULE_P(retval);
}

PG_FUNCTION_INFO_V1 (pgchem_r_reaction_to_smiles);

Datum
pgchem_r_reaction_to_smiles (PG_FUNCTION_ARGS)
{
  REACTION *reaction = PG_GETARG_REACTION_P (0);
  MOLECULE *tmpMol;
  char* offset = MOLARRAYPTR(reaction);
  text *result;
  char* tmpBuf;
  int i, size=0;
  
  for(i=0;i<reaction->num_products+reaction->num_reactants;i++) {
      size+=((MOLECULE*) offset)->sizesmi*sizeof(char);
      offset+=VARSIZE((MOLECULE*) offset);
  }
  
  offset = MOLARRAYPTR(reaction);   
  
  tmpBuf = (char *) palloc (size+sizeof(char));
  
  memset(tmpBuf,0x0,size+sizeof(char));

   for(i=0;i<reaction->num_reactants;i++) {
    tmpMol = (MOLECULE*) offset;
    if(strstr(SMIPTR(tmpMol),"\r\n") != NULL) {
    strncat(tmpBuf,SMIPTR(tmpMol),tmpMol->sizesmi-3);
} else if(strstr(SMIPTR(tmpMol),"\n") != NULL) {
    strncat(tmpBuf,SMIPTR(tmpMol),tmpMol->sizesmi-2);
}    
    if(i<reaction->num_reactants-1) strncat(tmpBuf,".",sizeof(char));
    offset+=VARSIZE(tmpMol)*sizeof(char);
}

strncat(tmpBuf,">>",2*sizeof(char));   

for(i=0;i<(reaction->num_products);i++) {
    tmpMol = (MOLECULE*) offset;
    if(strstr(SMIPTR(tmpMol),"\r\n") != NULL) {
    strncat(tmpBuf,SMIPTR(tmpMol),tmpMol->sizesmi-3);
} else if(strstr(SMIPTR(tmpMol),"\n") != NULL) {
    strncat(tmpBuf,SMIPTR(tmpMol),tmpMol->sizesmi-2);
}    
    if(i<reaction->num_products-1) strncat(tmpBuf,".",sizeof(char));
    offset+=VARSIZE(tmpMol)*sizeof(char);
}  

  result = (text *) palloc(strlen(tmpBuf)+VARHDRSZ);
  
  memset(result,0x0,strlen(tmpBuf)+VARHDRSZ);
  
  memcpy(VARDATA(result),tmpBuf,strlen(tmpBuf));
  
  SET_VARSIZE(result,strlen(tmpBuf)+VARHDRSZ);
  
  pfree(tmpBuf);

  PG_RETURN_TEXT_P (result);
}

PG_FUNCTION_INFO_V1 (pgchem_tversky);

Datum pgchem_tversky (PG_FUNCTION_ARGS)
{
  MOLECULE *mol1_prototype = PG_GETARG_MOLECULE_P (0);
  MOLECULE *mol2_variant = PG_GETARG_MOLECULE_P (1);
  float8 alpha_prototype = PG_GETARG_FLOAT8 (2);
  float8 beta_variant = PG_GETARG_FLOAT8 (3);

  PG_RETURN_FLOAT8 (ob_tversky
		    ((uint8 *) mol1_prototype->fp,  (uint8 *) mol2_variant->fp, FPSIZE2*sizeof(uint32), alpha_prototype, beta_variant));
}




