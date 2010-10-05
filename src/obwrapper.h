/**********************************************************************
 * obwrapper.h OpenBabel wrapper functions
 *
 * Copyright (c) 2004,2005 by Ernst-G. Schmid
 * Copyright (c) 2004,2005 by Bayer Business Services GmbH
 * for explicitly marked functions
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

//#include "molecule/molecule.h"
#include "molecule/molecule_limits.h"

typedef struct 
{
    unsigned int idx;
    int hybridization;
    unsigned char atomicnum;
    unsigned short isotope;
    char formalcharge;
    unsigned char stereo;
    unsigned char spinmultiplicity;
    unsigned char aromatic;
} _ATOM;  

typedef struct 
{
    unsigned int beginidx;
    unsigned int endidx;
    unsigned char order;
    unsigned int stereo; 
    unsigned char aromatic;
    //int flags; 
} _BOND;

#ifdef __cplusplus
extern "C"
{
#endif
  char *ob_mol_to_smiles (char *molfile, int omit_iso_and_chiral_markings);
  char *ob_mol_to_canonical_smiles (char *molfile, int omit_iso_and_chiral_markings);
  char *ob_smiles_to_smiles (char *smiles, int omit_iso_and_chiral_markings);
  char *ob_smiles_to_canonical_smiles (char *smiles, int omit_iso_and_chiral_markings);
  char *ob_smiles_to_mol (char *smiles);
  char *ob_inchi_to_mol (char *inchi);
  char *ob_smiles_to_inchi (char *smiles);
  char *ob_smiles_to_inchikey (char *smiles);
  char *ob_molfile_to_inchikey (char *molfile);
  double ob_MR (char *smiles);
  double ob_PSA (char *smiles);
  double ob_logP (char *smiles);
  unsigned int ob_num_H_donors (char *smiles);
  unsigned int ob_num_H_acceptors (char *smiles);
  double ob_tanimoto (unsigned char *fp1, unsigned char *fp2, unsigned short size);
  double ob_tversky (unsigned char *fp1_prototype, unsigned char *fp2_variant, unsigned short size, double alpha_prototype, double beta_variant);
  //double ob_tanimoto_n (unsigned char *fp1, unsigned char *fp2);
  unsigned int ob_popcount (unsigned char *fp, unsigned short size);
  void ob_fp2 (char *molfile, unsigned int *fp);
  void ob_fp3_bin (char *serializedInput, unsigned int *fp);
  void ob_fp3 (char *molfile, unsigned int *fp);
  void ob_fp_MACCS_bin (char *serializedInput, unsigned int *fp);
  void ob_fp_MACCS (char *molfile, unsigned int *fp);
  void ob_fp_bin (char *serializedInput, unsigned int *fp);
  int ob_total_charge (char *smiles);
  unsigned int ob_num_rotatable_bonds (char *smiles);
  unsigned int ob_num_bonds (char *smiles);
  unsigned int ob_num_atoms (char *smiles);
  unsigned int ob_num_heavy_atoms (char *smiles);
  unsigned int ob_is_chiral (char *molfile);
  unsigned int ob_is_nostruct (char *molfile);
  char *ob_add_hydrogens (char *smiles, int polaronly, int correct4PH);
  char *ob_delete_hydrogens (char *smiles, int nonpolaronly);
  char *ob_strip_salts (char *smiles, int neutralize_residue);
  int ob_2D (char *molfile);
  int ob_3D (char *molfile);
  char *ob_mol_to_V2000 (char *molfile);
  char *ob_mol_to_V3000 (char *molfile);
  char *ob_V3000_to_mol (char *molfile);
  //unsigned int *ob_efa_array (char *smiles);
  //unsigned int ob_SSS_SMARTS (const char *smarts_pattern, char *molfile);
  /*int ob_SSS_VF2 (char *molfile_sarg,
		 		 		 		      char *molfile);
  int ob_SSS_VF2_bin (char *serializedInput_sarg,
		 		 		 		      char *serializedInput);
		 		 		 		      
  int ob_ESS_VF2_bin (char *serializedInput_sarg,
		 		 		 		      char *serializedInput);*/
  
  
  int ob_SSS_SMARTS_native_bin (const char *smarts_pattern,
		 		 		 		      char *serializedInput);
		 		 		 		      
  int ob_SSS_SMARTS_native (const char *smarts_pattern,
		 		 		 		      char *smiles);
		 		 		 		      
  int ob_SSS_SMARTS_native_count_bin (const char *smarts_pattern, char *serializedInput);		 		 		 		      
		 		 		 		      
  int ob_SSS_SMARTS_native_count (const char *smarts_pattern,
		 		 		 		      char *smiles);
/*****************************************************************
 * This function Copyright (c) 2004
 * by Bayer Business Services GmbH
 *****************************************************************/
  double ob_molweight (char *smiles);
/*****************************************************************
 * This function Copyright (c) 2004
 * by Bayer Business Services GmbH
 *****************************************************************/
  char *ob_hillformula (char *smiles);
/*****************************************************************
 * This function Copyright (c) 2005
 * by Bayer Business Services GmbH
 *****************************************************************/
  double ob_exactmass (char *smiles);
  
  char *ob_lyophilize_molecule(char* smiles);
#ifdef __cplusplus
}
#endif

