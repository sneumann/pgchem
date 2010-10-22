/************************************************************************
 * molecule_io.c molecule input/output support functions
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
#include "libpq/md5.h"
#include "fmgr.h"
#include "libpq/pqformat.h"	/* needed for send/recv functions */
#include "molecule.h"
#include "obwrapper.h"
#include "barsoi/barsoi.h"

Datum molecule_in (PG_FUNCTION_ARGS);
Datum molecule_in_text (PG_FUNCTION_ARGS);
Datum molecule_in_varchar (PG_FUNCTION_ARGS);
Datum molecule_in_bytea (PG_FUNCTION_ARGS);
Datum molecule_out (PG_FUNCTION_ARGS);
Datum molecule_recv (PG_FUNCTION_ARGS);
Datum molecule_send (PG_FUNCTION_ARGS);
Datum pgchem_strip_salts (PG_FUNCTION_ARGS);
Datum pgchem_add_hydrogens (PG_FUNCTION_ARGS);
Datum pgchem_remove_hydrogens (PG_FUNCTION_ARGS);


MOLECULE *
new_molecule (char *smiles, char *molfile)
{
  unsigned int sizemf;
  unsigned int sizesmi;
  size_t totalsize;
  MOLECULE *result;
  char *inchikey = NULL;
  char *ancillarydata = NULL;
  char *aidata = NULL;
  uint32 ancsize = 0;
  
  ancillarydata = ob_lyophilize_molecule(smiles);
  
  if(ancillarydata == NULL) {
    elog (ERROR, "Molecule lyophilization failed! SMILES:\n %s", smiles);
  }    
  
  ancsize = *(unsigned int*) ancillarydata;
  sizemf = strlen (molfile)+1;
  sizesmi = strlen (smiles)+1;
  totalsize = CALCDATASZ (sizemf, sizesmi, ancsize);

  result = (MOLECULE *) palloc (totalsize);
  memset (result, 0x0, totalsize);

  if (strchr (smiles, '.') != NULL)
    result->disconnected = true;

  result->sizemf = sizemf;
  result->sizesmi = sizesmi;

  strncpy (SMIPTR(result), smiles, sizesmi);
  strncpy (MFPTR(result), molfile, sizemf);
  
  aidata = (char*) &((unsigned int*)ancillarydata)[1];
  memcpy(ANCPTR(result), aidata, ancsize);

  inchikey = ob_smiles_to_inchikey (smiles);
    
  if(inchikey == NULL || strlen(inchikey) != INCHIKEYSZ) {
    elog (ERROR, "InChI key generation failled ! SMILES:\n %s", smiles);
  } else {
    memcpy(result->inchikey, inchikey, INCHIKEYSZ);
  }   

  if (inchikey != NULL) {
    free (inchikey);
    inchikey=NULL;
  }
  
  ob_fp_bin(aidata, result->fp);
  
  if(ancillarydata != NULL) {
    free(ancillarydata);
    ancillarydata=NULL;
  }  

  SET_VARSIZE (result,totalsize);
  return result;
}

static MOLECULE *make_molecule(char *raw_input, int size) {
  MOLECULE *result;
  char *input = NULL;
  char *molfile = NULL;
  char *smiles = NULL;
  char *endptr;
  unsigned int new_len;
  
  if(strstr (raw_input, "M  END") != NULL) {
    input = palloc (size+sizeof(char));
    memcpy (input, raw_input, size);
    endptr = strstr (input, "M  END") + strlen("M  END")*sizeof(char);
    *endptr = 0x0;
    new_len = strlen(input);
    pfree (input);
    input=NULL;

    input = palloc (new_len + 1);
    strncpy(input,raw_input,new_len);
    input[new_len] = 0x0;
  } else {
    input = palloc (size+1);
    memcpy (input, raw_input, size);
    input[size]=0x0;
  }    
  
  if (strstr (input, "V2000") == NULL || strstr (input, "M  END") == NULL)
    {				//TODO: Bad Hack. If more fmts are needed, we need a dedicated converter function
      if (strstr (input, "V3000") != NULL && strstr (input, "M  END") != NULL) //V3000?
	{
	  molfile = ob_V3000_to_mol (input);
	  
	  if(molfile == NULL || !strlen(molfile) || strstr(molfile,"V3000")==NULL) {
	    if(molfile!=NULL) {
	      free (molfile);
	      molfile=NULL;
	    }
	    elog (ERROR, "Molfile generation failed! V3000:\n %s",input);
	  }  
  
	  smiles = ob_mol_to_smiles (input,0);
	  
	  if(smiles == NULL || !strlen(smiles)) {
	    elog (ERROR, "SMILES generation failed! Offender was :\n %s",input);
	  }  else if (!strlen(smiles)) {
	    elog (WARNING, "SMILES generation failed! Trying fallback...");
	    free (smiles);
	    smiles=NULL;

	    smiles = ob_mol_to_canonical_smiles (input,0);
            if(smiles == NULL) {
	      elog (ERROR, "Canonical SMILES generation finally failed! Offender was :\n %s",input);
            } else if (!strlen(smiles)) {
	      free(smiles);
	      smiles=NULL;
	      elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",input);
	    }  
	    elog (WARNING, "Fallback OK"); 
	  }  
	} else if (strstr (input, "InChI=") != NULL) { //InChI?	
	  molfile = ob_inchi_to_mol (input);
	  
	  if(molfile == NULL || !strlen(molfile) || strstr(molfile,"V2000")==NULL) {
	    if(molfile!=NULL) free (molfile);
	    elog (ERROR, "Molfile generation failed! InChI:\n %s",input);
	  }  
      
	  smiles = ob_mol_to_smiles (molfile,0);
	  
	  if(smiles == NULL || !strlen(smiles)) {
	    elog (ERROR, "SMILES generation failed! Offender was :\n %s",input);
	  }  else if (!strlen(smiles)) {
	    elog (WARNING, "SMILES generation failed! Trying fallback...");
	    free (smiles);
	    smiles=NULL;
	    smiles = ob_mol_to_canonical_smiles (input,0);
            if(smiles == NULL) {
	      elog (ERROR, "Canonical SMILES generation finally failed! InChI:\n %s",input);
            } else if (!strlen(smiles)) {
	      free(smiles);
	      smiles=NULL;
	      elog (ERROR, "SMILES generation finally failed! InChi:\n %s",input);
	    } 
            elog (NOTICE, "Fallback OK");   
	  }  
      } else { //SMILES?       
	  molfile = ob_smiles_to_mol (input);
	  
	  if(molfile == NULL || !strlen(molfile)) {
	    if(molfile!=NULL) free (molfile);
	    elog (ERROR, "Molfile generation failed! SMILES:\n %s",input);
	  }  

	  // smiles will be free()d below, but input was palloc()ed.
	  smiles = strdup(input);
      }
    }
  else //V2000
    {
      smiles = ob_mol_to_smiles (input,0);
      
      if(smiles == NULL || !strlen(smiles)) {
	elog (ERROR, "SMILES generation failed! Offender was :\n %s",input);
      }  else if (!strlen(smiles)) {
	elog (WARNING, "SMILES generation failed! Trying fallback...");
	free (smiles);
	smiles=NULL;

	smiles = ob_mol_to_canonical_smiles (input,0);
	if(smiles == NULL) {
	  elog (ERROR, "Cacnonical SMILES generation finally failed! Offender was :\n %s",input);
	} else if (!strlen(smiles)) {
	  free(smiles);
	  smiles=NULL;
	  elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",input);
	} 
	elog (NOTICE, "Fallback OK");   
      }  
      molfile = strdup(input);
    }

  if (smiles==NULL || molfile==NULL) {
    elog (ERROR,
	  "Input is not a V2000/V3000 molfile or InChI or SMILES: %s",
	  input);
  }

  result = new_molecule (smiles, molfile);

  if (smiles != NULL) {
	free (smiles);
	smiles=NULL;
  }
  if (molfile != NULL) {
	free (molfile);
	molfile=NULL;
  }    
      
  if (input != NULL) {
    pfree(input);
    input=NULL;
  }
    
  return result;
}    

/*
 * Convert a molecule in text form (V2000, SMILES, InChI) into a molecule
 */
PG_FUNCTION_INFO_V1 (molecule_in);

Datum
molecule_in (PG_FUNCTION_ARGS)
{
  char *input = PG_GETARG_CSTRING (0);
  int size = strlen(input);
  
  //printf("molecule_in\n");
 
  PG_RETURN_MOLECULE_P (make_molecule(input,size));
}

PG_FUNCTION_INFO_V1 (molecule_in_text);

Datum
molecule_in_text (PG_FUNCTION_ARGS)
{
  char *input = VARDATA(PG_GETARG_TEXT_P (0));
  int size = VARSIZE(PG_GETARG_TEXT_P (0))-VARHDRSZ;
 
  PG_RETURN_MOLECULE_P (make_molecule(input,size));
}

PG_FUNCTION_INFO_V1 (molecule_in_varchar);

Datum
molecule_in_varchar (PG_FUNCTION_ARGS)
{
  VarChar *x = PG_GETARG_VARCHAR_P (0);  
  char *input = VARDATA(x);
  int size = VARSIZE(x)-VARHDRSZ;
  
  PG_RETURN_MOLECULE_P (make_molecule(input,size));
}

PG_FUNCTION_INFO_V1 (molecule_in_bytea);

Datum
molecule_in_bytea (PG_FUNCTION_ARGS)
{
  char *input = VARDATA(PG_GETARG_BYTEA_P (0));
  int size = VARSIZE(PG_GETARG_BYTEA_P (0))-VARHDRSZ;
  
  PG_RETURN_MOLECULE_P (make_molecule(input,size));
}

/*
 * Output a molecule in cstring form as molfile
 */
PG_FUNCTION_INFO_V1 (molecule_out);

Datum
molecule_out (PG_FUNCTION_ARGS)
{
  MOLECULE *molecule = PG_GETARG_MOLECULE_P (0);
  //int sizesmi = molecule->sizesmi;
  //int i;
  //unsigned int *arr = (unsigned int*) EFAPTR(molecule);

  char *result = (char *) palloc (molecule->sizesmi);
  
  memset(result,0x0,molecule->sizesmi);

  strncpy (result, SMIPTR(molecule), molecule->sizesmi);
  
  //for(i=0;i<arr[0];i++) printf("%d ",arr[i]);
  //printf("\n");

  PG_RETURN_CSTRING (result);
}

/*****************************************************************************
 * Binary Input/Output functions
 *
 * These are optional.
 *****************************************************************************/

PG_FUNCTION_INFO_V1 (molecule_recv);

Datum
molecule_recv (PG_FUNCTION_ARGS)
{
  StringInfo buf = (StringInfo) PG_GETARG_POINTER (0);
  int len = buf->len;
  const char *str = pq_getmsgbytes (buf, len);
  MOLECULE *result = (MOLECULE *) palloc (len);

  //SET_VARSIZE (result,(buf->len + VARHDRSZ));
  
  memset(result,0x0,len);

  memcpy (result, str, len);

  PG_RETURN_POINTER (result);
}

PG_FUNCTION_INFO_V1 (molecule_send);

Datum
molecule_send (PG_FUNCTION_ARGS)
{
  MOLECULE *molecule = PG_GETARG_MOLECULE_P (0);
    
  /*StringInfoData buf;

    pq_begintypsend(&buf);
	
    pq_sendbytes(&buf,(const char*) molecule,sizeof(molecule));
	
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));*/
	
  PG_RETURN_BYTEA_P(molecule);
}

/*PG_FUNCTION_INFO_V1 (pgchem_V3000_to_molecule);

  Datum
  pgchem_V3000_to_molecule (PG_FUNCTION_ARGS)
  {
  MOLECULE *retval;
  text *arg_V3000;
  char *tmpV3000;
  char *molfile;
  char *smiles;
  int V3000_string_len;

  arg_V3000 = PG_GETARG_TEXT_P (0);

  V3000_string_len = VARSIZE (arg_V3000) - VARHDRSZ;

  tmpV3000 = (char *) palloc (V3000_string_len + 1);
  tmpV3000[0] = '\0';

  strncat (tmpV3000, VARDATA (arg_V3000), V3000_string_len);
  molfile = ob_V3000_to_mol (tmpV3000);
  smiles = ob_mol_to_canonical_smiles (molfile,1);

  retval = new_molecule(smiles,molfile);
  
  free(molfile);
  free(smiles);

  PG_RETURN_MOLECULE_P (retval);
  }

  PG_FUNCTION_INFO_V1 (pgchem_inchi_to_molecule);

  Datum
  pgchem_inchi_to_molecule (PG_FUNCTION_ARGS)
  {
  MOLECULE *retval;
  text *arg_inchi;
  char *tmpinchi;
  char *molfile;
  char *smiles;
  int inchi_string_len;

  arg_inchi = PG_GETARG_TEXT_P (0);

  inchi_string_len = VARSIZE (arg_inchi) - VARHDRSZ;

  tmpinchi = (char *) palloc (inchi_string_len + 1);
  tmpinchi[0] = '\0';

  strncat (tmpinchi, VARDATA (arg_inchi), inchi_string_len);
  molfile = ob_inchi_to_mol (tmpinchi);
  smiles = ob_mol_to_canonical_smiles (molfile,1);

  retval = new_molecule(smiles,molfile);
  
  free(molfile);
  free(smiles);

  PG_RETURN_MOLECULE_P (retval);
  } */

/*
 * Strip smaller fragments from a molecule, leaving only the larggest one
 */
PG_FUNCTION_INFO_V1 (pgchem_strip_salts);

Datum
pgchem_strip_salts (PG_FUNCTION_ARGS)
{
  char *molfile = NULL;
  char *smiles = NULL;
  //unsigned int *efa_array = NULL;
  MOLECULE *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  bool neutralize_residue = PG_GETARG_BOOL (1);

  smiles = ob_strip_salts (SMIPTR (arg_molecule), neutralize_residue ? 1 : 0);
  
  molfile = ob_smiles_to_mol (smiles);
  
  if(molfile == NULL || !strlen(molfile)) {
    elog (ERROR, "Molfile generation failed! Offender was :\n %s",smiles);
  } 
   
  //efa_array = ob_efa_array(smiles);           

  retval = new_molecule (smiles, molfile);

  free (molfile);
  free (smiles);
  //free (efa_array);

  PG_RETURN_MOLECULE_P (retval);
}

/*
 * Add hydrogens to a molecule
 */
PG_FUNCTION_INFO_V1 (pgchem_add_hydrogens);

Datum
pgchem_add_hydrogens (PG_FUNCTION_ARGS)
{
  char *molfile = NULL;
  char *smiles = NULL;
  //unsigned int *efa_array = NULL;
  MOLECULE *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  bool polaronly = PG_GETARG_BOOL (1);
  bool correct4PH = PG_GETARG_BOOL (2);

  smiles =
    ob_add_hydrogens (SMIPTR (arg_molecule),
		      polaronly ? 1 : 0, correct4PH ? 1 : 0);

  molfile = ob_smiles_to_mol (smiles);
  
  if(molfile == NULL || !strlen(molfile)) {
    elog (ERROR, "Molfile generation failed! Offender was :\n %s",smiles);
  } 

  //efa_array = ob_efa_array(smiles);           

  retval = new_molecule (smiles, molfile);

  free (molfile);
  free (smiles);
  //free (efa_array);

  PG_RETURN_MOLECULE_P (retval);
}

/*
 * Remove hydrogens from a molecule
 */
PG_FUNCTION_INFO_V1 (pgchem_remove_hydrogens);

Datum
pgchem_remove_hydrogens (PG_FUNCTION_ARGS)
{
  char *molfile = NULL;
  char *smiles = NULL;
  //unsigned int *efa_array = NULL;
  MOLECULE *retval;
  MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
  bool nonpolaronly = PG_GETARG_BOOL (1);

  smiles =
    ob_delete_hydrogens (SMIPTR (arg_molecule),
			 nonpolaronly ? 1 : 0);

  molfile = ob_smiles_to_mol (smiles);
  
  if(molfile == NULL || !strlen(molfile)) {
    elog (ERROR, "Molfile generation failed! Offender was :\n %s",smiles);
  } 
  //efa_array = ob_efa_array(smiles);           

  retval = new_molecule (smiles, molfile);

  free (molfile);
  free (smiles);
  //free (efa_array);

  PG_RETURN_MOLECULE_P (retval);
}
