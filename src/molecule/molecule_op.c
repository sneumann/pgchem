/************************************************************************
 * molecule_op.c molecule operator support functions
 *
 * Copyright (c) 2007,2009 by Ernst-G. Schmid
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
#include "postgres.h"
#include "executor/executor.h"
#include "fmgr.h"
#include "molecule.h"
#include "barsoi/barsoi.h"
#include "obwrapper.h"

PG_FUNCTION_INFO_V1 (molecule_alwaystrue);
PG_FUNCTION_INFO_V1 (molecule_contained_in);
PG_FUNCTION_INFO_V1 (molecule_contains);
PG_FUNCTION_INFO_V1 (molecule_equals);
PG_FUNCTION_INFO_V1 (molecule_similarity);
PG_FUNCTION_INFO_V1 (molecule_similarity_gist);

Datum molecule_alwaystrue (PG_FUNCTION_ARGS);
Datum molecule_contained_in (PG_FUNCTION_ARGS);
Datum molecule_contains (PG_FUNCTION_ARGS);
Datum molecule_equals (PG_FUNCTION_ARGS);
Datum molecule_similarity (PG_FUNCTION_ARGS);
Datum molecule_similarity_gist (PG_FUNCTION_ARGS);

/*
* Always returns TRUE for debugging purposes
*/
Datum
molecule_alwaystrue (PG_FUNCTION_ARGS)
{
 PG_RETURN_BOOL (true);
}

/*
* Check if a query molecule is contained in a predicate molecule by performing a graph isomorphism check.
* If the query molecule is disconnected (has fragments), the check is done with checkmol/barsoi, bacause the OpenBabel matcher does not support this.
* Otherwise the faster openBabel matcher is used.
*/
Datum
molecule_contained_in (PG_FUNCTION_ARGS)
{
  MOLECULE *query = PG_GETARG_MOLECULE_P (0);
  MOLECULE *predicate = PG_GETARG_MOLECULE_P (1);
  int match;//, i, j, offset = 1;
  char* querysmi = SMIPTR(query);
  //unsigned int *efa_q, *efa_p;

  //printf("Query: %s\n Pred: %s\n",query->data, predicate->data);

/*  if (query->isbz == true)
    {
      if (predicate->nobz == true)
	PG_RETURN_BOOL (false);
      else
	PG_RETURN_BOOL (true);
    } */

  if (query->disconnected == true)
    //{
      elog (ERROR, "Disconnected molecules as query input are not supported!");
      
      /*xm_set_ring_perception_algorithm (RPA_SAR);
      xm_set_strict_typing (FEATURE_OFF);
      mm_set_r_s_check (FEATURE_OFF);
      mm_set_e_z_check (FEATURE_OFF);
      mm_set_chg_check (FEATURE_ON);
      mm_set_iso_check (FEATURE_ON);
      mm_set_rad_check (FEATURE_ON);
      mm_set_exact_match (FEATURE_OFF);

      mm_set_mol (MFPTR(query));
      mm_set_current_mol_as_query ();
      
      printf("%s\n",SMIPTR(query));

      mm_set_mol (MFPTR(predicate));
      
      printf("%s\n",SMIPTR(predicate));

      if (mm_match () != 0) {
          printf("Match TRUE!\n");
	PG_RETURN_BOOL (true);
}	
    }
  else
    {*/
    
  /*  efa_q = (unsigned int*) EFAPTR(query);
    efa_p = (unsigned int*) EFAPTR(predicate);
    
   if(efa_q[0] > efa_p[0]) {printf("EFA kill type 1"); PG_RETURN_BOOL (false);}
    else {
        for(i=1;i<efa_q[0];i+=2) {
            match = 0;
            for(j=offset;j<efa_p[0];j+=2) {
                if(efa_q[i]==efa_p[j] && efa_q[++i] <= efa_p[++j]) {match = 1; offset=j; break;}
                //printf("%d %d %d %d\n",efa_q[i],efa_p[j],efa_q[i+1],efa_p[j+1]);
            }
            if(match==0) {printf("EFA kill type 2");PG_RETURN_BOOL (false);}    
        }    
    }  */
    
    if(strstr(querysmi,"@")!=NULL || strstr(querysmi,"/")!=NULL || strstr(querysmi,"\\")!=NULL)
    match = ob_SSS_SMARTS_native (querysmi, SMIPTR(predicate));
    else
    match = ob_SSS_SMARTS_native_bin (querysmi, ANCPTR(predicate));
    
     if(match<0) elog (ERROR, "Invalid SMARTS pattern: %s",SMIPTR(query));
    
      if (match != 0)
	PG_RETURN_BOOL (true);
    //}
//printf("EFA kill type 3\n");
  PG_RETURN_BOOL (false);
}

/*
* Check if a query molecule is contained in a predicate molecule by performing a graph isomorphism check.
* If the query molecule is disconnected (has fragments), the check is done with checkmol/barsoi, bacause the OpenBabel matcher does not support this.
* Otherwise the faster openBabel matcher is used.
*/
Datum
molecule_contains (PG_FUNCTION_ARGS)
{
  MOLECULE *query = PG_GETARG_MOLECULE_P (1);
  MOLECULE *predicate = PG_GETARG_MOLECULE_P (0);
  int match;
  char* querysmi = SMIPTR(query);
  //,i,j, offset=1;
  //unsigned int *efa_q, *efa_p;
  
 //printf("Query: %s\n Pred: %s\n",query->data, predicate->data);

 /* if (query->isbz == true)
    {
      if (predicate->nobz == true)
	PG_RETURN_BOOL (false);
      else
	PG_RETURN_BOOL (true);
    } */

  if (query->disconnected == true)
    //{
     elog (ERROR, "Disconnected molecules as query input are not supported!");
      
      /*xm_set_ring_perception_algorithm (RPA_SAR);
      xm_set_strict_typing (FEATURE_OFF);
      mm_set_r_s_check (FEATURE_OFF);
      mm_set_e_z_check (FEATURE_OFF);
      mm_set_chg_check (FEATURE_ON);
      mm_set_iso_check (FEATURE_ON);
      mm_set_rad_check (FEATURE_ON);
      mm_set_exact_match (FEATURE_OFF);

     mm_set_mol (MFPTR(query));
      mm_set_current_mol_as_query ();
      
      printf("%s\n",SMIPTR(query));

      mm_set_mol (MFPTR(predicate));
      
       printf("%s\n",SMIPTR(predicate));

      if (mm_match () != 0) {
      printf("Match TRUE!\n");
	PG_RETURN_BOOL (true);}
    }
  else
    {*/
    
    /*efa_q = (unsigned int*) EFAPTR(query);
    efa_p = (unsigned int*) EFAPTR(predicate);
    
    if(efa_q[0] > efa_p[0]) {printf("EFA kill type 1"); PG_RETURN_BOOL (false);}
    else {
        for(i=1;i<efa_q[0];i+=2) {
            match = 1;
            for(j=offset;j<efa_p[0];j+=2) {
                if(efa_q[i] == efa_p[j] && efa_q[i+1] <= efa_p[j+1]) {match = 1; offset=j; break;}
                
               // printf("%d %d %d %d\n",efa_q[i],efa_p[j],efa_q[i+1],efa_p[j+1]);
            }
            if(match==0) {printf("EFA kill type 2");PG_RETURN_BOOL (false);}    
        }    
    } */ 
    
    
    if(strstr(querysmi,"@")!=NULL || strstr(querysmi,"/")!=NULL || strstr(querysmi,"\\")!=NULL)
    match = ob_SSS_SMARTS_native (querysmi, SMIPTR(predicate));
    else
    match = ob_SSS_SMARTS_native_bin (querysmi, ANCPTR(predicate));
    
    if(match<0) elog (ERROR, "Invalid SMARTS pattern: %s",SMIPTR(query));
    
      if (match != 0)
	PG_RETURN_BOOL (true);
    //}

//printf("EFA kill type 3\n");
  PG_RETURN_BOOL (false);
}

/*
* Check if a query molecule is equal to a predicate molecule by comparing their molhashes.
*/
Datum
molecule_equals (PG_FUNCTION_ARGS)
{
  MOLECULE *query = PG_GETARG_MOLECULE_P (0);
  MOLECULE *predicate = PG_GETARG_MOLECULE_P (1);
  
 /* if (query->isbz == true) {
      if (predicate->isbz == true)
    PG_RETURN_BOOL (true);
      else
    PG_RETURN_BOOL (false);
    }    */
   //if(query->popcount != predicate->popcount) PG_RETURN_BOOL (false);
  
   if (memcmp (query->inchikey, predicate->inchikey, INCHIKEYSZ) != 0)
    PG_RETURN_BOOL (false);

  PG_RETURN_BOOL (true);
} 

/*
* Returns the Tanimoto similarity of two molecules.
*/
Datum
molecule_similarity (PG_FUNCTION_ARGS)
{
  MOLECULE *mol1 = PG_GETARG_MOLECULE_P (0);
  MOLECULE *mol2 = PG_GETARG_MOLECULE_P (1);
  
  /*unsigned int andbits = 0, orbits = 0;
  int andfp, orfp;
  unsigned int *fp1, *fp2;
  unsigned int i;
  
  fp1 = mol1->fp;
  fp2 = mol2->fp;
  
  for(i = 0; i < FPSIZE; i++) {
    andfp = fp1[i] & fp2[i];
    orfp = fp1[i] | fp2[i];
    for (; andfp; andfp = andfp<<1) {
      if (andfp < 0) { ++andbits; }
    }
    for (; orfp; orfp = orfp<<1) {
      if (orfp < 0) { ++orbits; }
    }
  }
  
  if (orbits == 0) {
    PG_RETURN_FLOAT8(0.0);
  }
  
  PG_RETURN_FLOAT8((double) andbits / (double) orbits);*/

  PG_RETURN_FLOAT8 (ob_tanimoto
		    ((uint8 *) mol1->fp,  (uint8 *) mol2->fp, FPSIZE2*sizeof(uint32)));
		    
  /*PG_RETURN_FLOAT8 (ob_tanimoto_n
		    ((uint8 *) mol1->fp,  (uint8 *) mol2->fp));	  */ 
		    
  /*PG_RETURN_FLOAT8 (ob_tanimoto
		    ((unsigned char *) mol1->(fp2+OFFSET),  (unsigned char *) mol2->(fp2+OFFSET), FPSIZE2*sizeof(unsigned int)));*/
		    
		    /*PG_RETURN_FLOAT8(0.0);*/
}

/*
* Returns the Tanimoto similarity of two molecules with index support.
*/
Datum
molecule_similarity_gist (PG_FUNCTION_ARGS)
{
  MOLECULE *mol1 = PG_GETARG_MOLECULE_P (0);
  HeapTupleHeader  t = PG_GETARG_HEAPTUPLEHEADER(1);
  bool isnull;
  float8 similarity;

  MOLECULE *mol2 = NULL;
  char *oper = NULL;
  float8 threshold;

  mol2 = (MOLECULE*) DatumGetPointer(GetAttributeByName(t, "_mq", &isnull));
  if(isnull) 
	  elog (ERROR, "Query molecule must not be NULL");

  oper = (char*) DatumGetPointer(GetAttributeByName(t, "_op", &isnull));
  if(isnull) 
	  elog (ERROR, "Query operator must not be NULL");

  threshold = DatumGetFloat8(GetAttributeByName(t, "_threshold", &isnull));
  if(isnull) 
	  elog (ERROR, "Query threshold must not be NULL");
  
  similarity = ob_tanimoto
		    ((uint8 *) mol1->fp,  (uint8 *) mol2->fp, FPSIZE2*sizeof(uint32));		    
		    
  if(oper[0]=='>' && oper[1]=='=') PG_RETURN_BOOL(similarity >= threshold);
  if(oper[0]=='<' && oper[1]=='=') PG_RETURN_BOOL(similarity <= threshold);
  if(oper[1]=='>') PG_RETURN_BOOL(similarity > threshold);
  if(oper[1]=='<') PG_RETURN_BOOL(similarity < threshold);
  if(oper[1]=='=') PG_RETURN_BOOL(similarity == threshold);

  PG_RETURN_BOOL(false);
}
