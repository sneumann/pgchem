/************************************************************************
 * reaction_op.c molecule operator support functions
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
#include "reaction.h"
//#include "barsoi/barsoi.h"
#include "obwrapper.h"
#include "fmgr.h"

PG_FUNCTION_INFO_V1 (reaction_alwaystrue);
PG_FUNCTION_INFO_V1 (reaction_contained_in);
PG_FUNCTION_INFO_V1 (reaction_contains);
PG_FUNCTION_INFO_V1 (reaction_equals);
PG_FUNCTION_INFO_V1 (reaction_equals_exact);
PG_FUNCTION_INFO_V1 (reaction_equals_products_exact);
PG_FUNCTION_INFO_V1 (reaction_similarity);
PG_FUNCTION_INFO_V1 (reaction_similarity_reactants);
PG_FUNCTION_INFO_V1 (reaction_similarity_products);

Datum reaction_alwaystrue (PG_FUNCTION_ARGS);
Datum reaction_contained_in (PG_FUNCTION_ARGS);
Datum reaction_contains (PG_FUNCTION_ARGS);
Datum reaction_equals_exact (PG_FUNCTION_ARGS);
Datum reaction_equals_products_exact (PG_FUNCTION_ARGS);
Datum reaction_equals (PG_FUNCTION_ARGS);
Datum reaction_similarity (PG_FUNCTION_ARGS);
Datum reaction_similarity_reactants (PG_FUNCTION_ARGS);
Datum reaction_similarity_products (PG_FUNCTION_ARGS);

static Datum rss_match(REACTION *query, REACTION *predicate) {
  MOLECULE *tmpMol_q = NULL, *tmpMol_p = NULL;
  char *offset_q = MOLARRAYPTR(query);
  char *offset_p = MOLARRAYPTR(predicate);
  char *offset_q_products = NULL; 
  int i,j;//,query_len = predicate->num_products+predicate->num_reactants;
  uint8 r_matches[predicate->num_reactants];
  uint8 p_matches[predicate->num_products];
  int match = 0;
 
  if(query->num_products > predicate->num_products || query->num_reactants > predicate->num_reactants) PG_RETURN_BOOL (false);
  
  if(query->num_products +  predicate->num_products + query->num_reactants + predicate->num_reactants == 0) PG_RETURN_BOOL (true);
  
  if(query->num_products + query->num_reactants == 0) PG_RETURN_BOOL (false);
  
  memset(&r_matches,0x0,predicate->num_reactants*sizeof(uint8));
  memset(&p_matches,0x0,predicate->num_products*sizeof(uint8));
  
  if (query->num_reactants != 0) {
  for (i=0;i<predicate->num_reactants;i++) {
      
      offset_q = MOLARRAYPTR(query);
      tmpMol_p = (MOLECULE*) offset_p;
      
      for(j=0;j<query->num_reactants;j++) {
          
            tmpMol_q = (MOLECULE*) offset_q;
            
            if (tmpMol_q->disconnected == true) elog (ERROR, "Disconnected reaction elements as query input are not supported!");
            
            match = ob_SSS_SMARTS_native (SMIPTR(tmpMol_q), SMIPTR(tmpMol_p));
            
            if(match<0) elog (ERROR, "Invalid SMARTS pattern: %s",SMIPTR(tmpMol_q));
      
            if ( match != 0) {
                r_matches[i]=1;
            }    
      
            offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }
      
      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
  }
}  else {
    for (i=0;i<predicate->num_reactants;i++) {
        tmpMol_p = (MOLECULE*) offset_p;
        offset_p+=VARSIZE(tmpMol_p)*sizeof(char);  //Must advance to products in any case
}    
} 

 j=0;
  
  for(i=0;i<predicate->num_reactants;i++) {
      if(r_matches[i] == 1) j++;
  }
  
  if(j<query->num_reactants) PG_RETURN_BOOL (false);     
  
  offset_q_products = offset_q;
  //offset_p = MOLARRAYPTR(predicate)+predicate->num_reactants*sizeof(MOLECULE*);
   
  if(query->num_products != 0) { 
  for (i=0;i<predicate->num_products;i++) {
      
      //offset_q = MOLARRAYPTR(query);
      offset_q = offset_q_products;
      tmpMol_p = (MOLECULE*) offset_p;
      
      /*for(j=0;j<query->num_reactants;j++) {
          tmpMol_q = (MOLECULE*) offset_q;
          offset_q+=tmpMol_q->len*sizeof(char);
          } */           
     
      for(j=0;j<query->num_products;j++) {
          
            tmpMol_q = (MOLECULE*) offset_q;
            
            if (tmpMol_q->disconnected == true) elog (ERROR, "Disconnected reaction elements as query input are not supported!");
            
            match = ob_SSS_SMARTS_native (SMIPTR(tmpMol_q), SMIPTR(tmpMol_p));
            
            if(match<0) elog (ERROR, "Invalid SMARTS pattern: %s",SMIPTR(tmpMol_q));
      
            if (match != 0) {
                p_matches[i]=1;
            }    
      
            offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }
      
      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
  } 
}  

  j=0;
  
  for(i=0;i<predicate->num_products;i++) {
      if(p_matches[i] == 1) j++;
  }
  
  if(j<query->num_products) PG_RETURN_BOOL (false);  
  
  PG_RETURN_BOOL (true);
}  

/*
* Always returns TRUE for debugging purposes
*/
Datum
reaction_alwaystrue (PG_FUNCTION_ARGS)
{
 PG_RETURN_BOOL (true);
}

/*
* Check if a query reaction is contained in a predicate reaction by performing a graph isomorphism check.
* If the query reaction is disconnected (has fragments), the check is done with checkmol/barsoi, bacause the OpenBabel matcher does not support this.
* Otherwise the faster openBabel matcher is used.
*/
Datum
reaction_contained_in (PG_FUNCTION_ARGS)
{
  REACTION *query = PG_GETARG_REACTION_P (0);
  REACTION *predicate = PG_GETARG_REACTION_P (1);
    
      
  PG_RETURN_BOOL (rss_match(query,predicate));
}

/*
* Check if a query reaction is contained in a predicate reaction by performing a graph isomorphism check.
* If the query reaction is disconnected (has fragments), the check is done with checkmol/barsoi, bacause the OpenBabel matcher does not support this.
* Otherwise the faster openBabel matcher is used.
*/
Datum
reaction_contains (PG_FUNCTION_ARGS)
{
  REACTION *query = PG_GETARG_REACTION_P (1);
  REACTION *predicate = PG_GETARG_REACTION_P (0);
  
  PG_RETURN_BOOL (rss_match(query,predicate));
}

/*
* Check if a query reaction is exactly equal to a predicate reaction by comparing their hashes and positions.
*/
Datum
reaction_equals_exact (PG_FUNCTION_ARGS)
{
  REACTION *query = PG_GETARG_REACTION_P (0);
  REACTION *predicate = PG_GETARG_REACTION_P (1);
  MOLECULE *tmpMol_q, *tmpMol_p;
  char *offset_q = MOLARRAYPTR(query);
  char *offset_p = MOLARRAYPTR(predicate);
  int i;
  
  if(query->num_products != predicate->num_products || query->num_reactants != predicate->num_reactants) PG_RETURN_BOOL (false);
  
  for (i=0;i<query->num_products+query->num_reactants;i++) {
      tmpMol_q = (MOLECULE*) offset_q;
      tmpMol_p = (MOLECULE*) offset_p;
      
      if (memcmp (tmpMol_q->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) != 0) PG_RETURN_BOOL (false);
      
      offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
  }    
      
  PG_RETURN_BOOL (true);
} 

/*
* Check if a query reaction is exactly equal on the product side to a predicate reaction by comparing all hashes and positions of the products.
*/
Datum
reaction_equals_products_exact (PG_FUNCTION_ARGS)
{
  REACTION *query = PG_GETARG_REACTION_P (0);
  REACTION *predicate = PG_GETARG_REACTION_P (1);
  MOLECULE *tmpMol_q = NULL, *tmpMol_p = NULL;
  char *offset_q = MOLARRAYPTR(query);
  char *offset_p = MOLARRAYPTR(predicate);
  int i,j;
  uint8 r_matches[predicate->num_reactants];
  
  if(query->num_products != predicate->num_products || query->num_reactants != predicate->num_reactants) PG_RETURN_BOOL (false);
  
  memset(&r_matches,0x0,predicate->num_reactants*sizeof(uint8));
  
  if (query->num_reactants != 0) {
  for (i=0;i<predicate->num_reactants;i++) {
      
      offset_q = MOLARRAYPTR(query);
      tmpMol_p = (MOLECULE*) offset_p;
      
      for(j=0;j<query->num_reactants;j++) {
          
            tmpMol_q = (MOLECULE*) offset_q;
      
            if (memcmp (tmpMol_q->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) == 0) {
                r_matches[i]=1;
            }    
      
            offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }
      
      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
  }
}  else {
    for (i=0;i<predicate->num_reactants;i++) {
        tmpMol_p = (MOLECULE*) offset_p;
        offset_p+=VARSIZE(tmpMol_p)*sizeof(char);  //Must advance to products in any case
}    
} 

 j=0;
  
  for(i=0;i<predicate->num_reactants;i++) {
      if(r_matches[i] == 1) j++;
  }
  
  if(j!=query->num_reactants) PG_RETURN_BOOL (false);
  
   /*for(j=0;j<query->num_reactants;j++) {
          tmpMol_q = (MOLECULE*) offset_q;
          offset_q+=tmpMol_q->len*sizeof(char);
          }*/
  
   if(query->num_products != 0) { 
  for (i=0;i<query->num_products;i++) {
      tmpMol_q = (MOLECULE*) offset_q;
      tmpMol_p = (MOLECULE*) offset_p;
      
      if (memcmp (tmpMol_q->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) != 0) PG_RETURN_BOOL (false);
      
      offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
  } 
}     

  PG_RETURN_BOOL (true);
} 

/*
* Check if a query reaction is equal on the product side to a predicate reaction by comparing all hashes.
*/
Datum
reaction_equals (PG_FUNCTION_ARGS)
{
  REACTION *query = PG_GETARG_REACTION_P (0);
  REACTION *predicate = PG_GETARG_REACTION_P (1);
  MOLECULE *tmpMol_q = NULL, *tmpMol_p = NULL;
  char *offset_q = MOLARRAYPTR(query);
  char *offset_p = MOLARRAYPTR(predicate);
  char *offset_q_products = NULL; 
  int i,j;//,query_len = predicate->num_products+predicate->num_reactants;
  uint8 r_matches[predicate->num_reactants];
  uint8 p_matches[predicate->num_products];
 
  if(query->num_products != predicate->num_products || query->num_reactants != predicate->num_reactants) PG_RETURN_BOOL (false);
  
  memset(&r_matches,0x0,predicate->num_reactants*sizeof(uint8));
  memset(&p_matches,0x0,predicate->num_products*sizeof(uint8));
  
  if (query->num_reactants != 0) {
  for (i=0;i<predicate->num_reactants;i++) {
      
      offset_q = MOLARRAYPTR(query);
      tmpMol_p = (MOLECULE*) offset_p;
      
      for(j=0;j<query->num_reactants;j++) {
          
            tmpMol_q = (MOLECULE*) offset_q;
      
            if (memcmp (tmpMol_q->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) == 0) {
                r_matches[i]=1;
            }    
      
            offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }
      
      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
  }
}  else {
    for (i=0;i<predicate->num_reactants;i++) {
        tmpMol_p = (MOLECULE*) offset_p;
        offset_p+=VARSIZE(tmpMol_p)*sizeof(char);  //Must advance to products in any case
}    
} 

 j=0;
  
  for(i=0;i<predicate->num_reactants;i++) {
      if(r_matches[i] == 1) j++;
  }
  
  if(j!=query->num_reactants) PG_RETURN_BOOL (false);     
  
  offset_q_products = offset_q;
   
  if(query->num_products != 0) { 
  for (i=0;i<predicate->num_products;i++) {
      
      offset_q = offset_q_products;
      tmpMol_p = (MOLECULE*) offset_p;   
     
      for(j=0;j<query->num_products;j++) {
          
            tmpMol_q = (MOLECULE*) offset_q;
      
            if ( memcmp (tmpMol_p->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) == 0) {
                p_matches[i]=1;
            }    
      
            offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }
      
      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
  } 
}  

  j=0;
  
  for(i=0;i<predicate->num_products;i++) {
      if(p_matches[i] == 1) j++;
  }
  
  if(j!=query->num_products) PG_RETURN_BOOL (false);  
  
  PG_RETURN_BOOL (true);
}  

/*
* Returns the Tanimoto similarity of two reactions.
*/
Datum
reaction_similarity (PG_FUNCTION_ARGS)
{
  REACTION *rxn1 = PG_GETARG_REACTION_P (0);
  REACTION *rxn2 = PG_GETARG_REACTION_P (1);
  

  PG_RETURN_FLOAT8 (ob_tanimoto
		    ((uint8 *) rxn1->fp,  (uint8 *) rxn2->fp, 2*FPSIZE2*sizeof(uint32)));
}

Datum
reaction_similarity_reactants (PG_FUNCTION_ARGS)
{
  REACTION *rxn1 = PG_GETARG_REACTION_P (0);
  REACTION *rxn2 = PG_GETARG_REACTION_P (1);
  

  PG_RETURN_FLOAT8 (ob_tanimoto
		    ((uint8 *) rxn1->fp,  (uint8 *) rxn2->fp, FPSIZE2*sizeof(uint32)));
}

Datum
reaction_similarity_products (PG_FUNCTION_ARGS)
{
  REACTION *rxn1 = PG_GETARG_REACTION_P (0);
  REACTION *rxn2 = PG_GETARG_REACTION_P (1);
  

  PG_RETURN_FLOAT8 (ob_tanimoto
		    ((uint8 *) rxn1->fp+FPSIZE2,  (uint8 *) rxn2->fp+FPSIZE2, FPSIZE2*sizeof(uint32)));
}
