/************************************************************************
 * obwrapper.cpp OpenBabel wrapper functions
 *
 * Copyright (c) 2004,2009 by Ernst-G. Schmid
 * Copyright (c) 2004,2009 by Bayer Business Services GmbH
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
//#include <vector>

#if HAVE_SSTREAM
#include <sstream>
#elif HAVE_SSTREAM_H
#include <sstream.h>
#endif

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#endif

#include <mol.h>
#include <fingerprint.h>
#include <obconversion.h>
#include <openbabel/groupcontrib.h>
#include <obiter.h>
//#include <vf2/vf2.h>
#include "obwrapper.h"

using namespace OpenBabel;
using namespace std;
//using namespace XchemTigress;

static const int popcount_counts[] = {
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
};

static bool ob_rehydrate_molecule(OBBase* pOb, char *serializedInput);

extern "C" char *
ob_mol_to_smiles (char *molfile, int omit_iso_and_chiral_markings)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream smilesstream;
  char *tmpSmiles;

  conv.SetInAndOutFormats ("MDL", "SMI");
  conv.AddOption ("n", OBConversion::OUTOPTIONS);
  
  if (omit_iso_and_chiral_markings != 0)
    {
      conv.AddOption ("i", OBConversion::OUTOPTIONS);
    }

  conv.Read (&mol, &molstream);

  if (mol.Empty ())
    return NULL;

  conv.Write (&mol, &smilesstream);

  outstring = smilesstream.str ();

  //outstring = outstring.substr (0, outstring.length ()-1);

  tmpSmiles = strdup (outstring.c_str ());
  
  //printf("Orig: %s\n",tmpSmiles);

  return (tmpSmiles);
}

extern "C" char *
ob_mol_to_canonical_smiles (char *molfile, int omit_iso_and_chiral_markings)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream smilesstream;
  char *tmpSmiles;

  conv.SetInAndOutFormats ("MDL", "CAN");
  conv.AddOption ("n", OBConversion::OUTOPTIONS);
  //conv.AddOption ("c", OBConversion::OUTOPTIONS);
  if (omit_iso_and_chiral_markings != 0)
    {
      conv.AddOption ("i", OBConversion::OUTOPTIONS);
    }

  conv.Read (&mol, &molstream);

  if (mol.Empty ())
    return NULL;

  conv.Write (&mol, &smilesstream);

  outstring = smilesstream.str ();

  outstring = outstring.substr (0, outstring.length () - 1);

  tmpSmiles = strdup (outstring.c_str ());

  return (tmpSmiles);
}

extern "C" char *
ob_smiles_to_smiles (char *smiles, int omit_iso_and_chiral_markings)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream smilesstream;
  char *tmpSmiles;

  conv.SetInAndOutFormats ("SMI", "SMI");
  conv.AddOption ("n", OBConversion::OUTOPTIONS);
  
  if (omit_iso_and_chiral_markings != 0)
    {
      conv.AddOption ("i", OBConversion::OUTOPTIONS);
    }

  conv.Read (&mol, &molstream);

  if (mol.Empty ())
    return NULL;

  conv.Write (&mol, &smilesstream);

  outstring = smilesstream.str ();

  //outstring = outstring.substr (0, outstring.length ()-1);

  tmpSmiles = strdup (outstring.c_str ());
  
  //printf("Orig: %s\n",tmpSmiles);

  return (tmpSmiles);
}

extern "C" char *
ob_smiles_to_canonical_smiles (char *smiles, int omit_iso_and_chiral_markings)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream smilesstream;
  char *tmpSmiles;

  conv.SetInAndOutFormats ("SMI", "CAN");
  conv.AddOption ("n", OBConversion::OUTOPTIONS);
  //conv.AddOption ("c", OBConversion::OUTOPTIONS);
  if (omit_iso_and_chiral_markings != 0)
    {
      conv.AddOption ("i", OBConversion::OUTOPTIONS);
    }

  conv.Read (&mol, &molstream);

  if (mol.Empty ())
    return NULL;

  conv.Write (&mol, &smilesstream);

  outstring = smilesstream.str ();

  outstring = outstring.substr (0, outstring.length () - 1);

  tmpSmiles = strdup (outstring.c_str ());

  return (tmpSmiles);
}

extern "C" char *
ob_smiles_to_mol (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  string outstring;
  istringstream smilesstream (tmpStr);
  ostringstream molstream;
  char *tmpMolfile;

  conv.SetInAndOutFormats ("SMI", "MDL");

  conv.Read (&mol, &smilesstream);

  if (mol.Empty ())
    return NULL;

  conv.Write (&mol, &molstream);

  outstring = molstream.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpMolfile = strdup (outstring.c_str ());

  return (tmpMolfile);
}

/*extern "C" unsigned int *
ob_efa_array (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream smilesstream (tmpStr);
  std::vector<unsigned int> atom_numbers;
  std::vector<unsigned int> efa;
  unsigned int atomic_num;
  unsigned int *efa_array = NULL;
  unsigned int i;
  unsigned int size = 1;
  bool endloop = false;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &smilesstream);
  
  printf("%d\n",sizeof(mol));

if (!mol.Empty()) {
    
  FOR_ATOMS_OF_MOL(atom, mol) {
      atomic_num=atom->GetAtomicNum();
         if(atomic_num>1) atom_numbers.push_back(atomic_num);
         }
         
  sort( atom_numbers.begin(), atom_numbers.end());

      vector<unsigned int>::const_iterator it=atom_numbers.begin();
      
      atomic_num = *it;
      
      while(true) { 
      i = 0;   

      while(*it == atomic_num) {
           i++;
           it++;
       if(it==atom_numbers.end()) {
           endloop = true;
           break;
       }    
    } 
    
     efa.push_back(atomic_num);
     efa.push_back(i);
     
     if(endloop) break;
     
     atomic_num = *it;  
  }
  
  size = efa.size()+1; 
}  
  
  efa_array = new unsigned int[size];
  
  efa_array[0] = size;
  
  if (size > 1) {
      memcpy((unsigned int*)&efa_array[1],(unsigned int*)&efa.front(),(size-1)*sizeof(unsigned int));
      }
  else
    cout << "Warning: EFA array contains no elements." << endl;
         
  return efa_array;
}*/

extern "C" char *
ob_inchi_to_mol (char *inchi)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (inchi);
  string outstring;
  istringstream inchistream (tmpStr);
  ostringstream molstream;
  char *tmpMolfile;

  conv.SetInAndOutFormats ("INCHI", "MDL");

  conv.Read (&mol, &inchistream);

  if (mol.Empty ())
    return NULL;

  conv.Write (&mol, &molstream);

  outstring = molstream.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpMolfile = strdup (outstring.c_str ());

  return (tmpMolfile);
}

extern "C" char *
ob_mol_to_V2000 (char *molfile)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream V2000stream;
  char *tmpV2000;

  conv.SetInAndOutFormats ("MDL", "MDL");
  conv.AddOption ("2", OBConversion::INOPTIONS);
  conv.AddOption ("2", OBConversion::OUTOPTIONS);

  conv.Read (&mol, &molstream);

  if (mol.Empty ())
    return NULL;

  conv.Write (&mol, &V2000stream);

  outstring = V2000stream.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpV2000 = strdup (outstring.c_str ());

  return (tmpV2000);
}

extern "C" char *
ob_mol_to_V3000 (char *molfile)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream V3000stream;
  char *tmpV3000;

  conv.SetInAndOutFormats ("MDL", "MDL");
  conv.AddOption ("2", OBConversion::INOPTIONS);
  conv.AddOption ("3", OBConversion::OUTOPTIONS);

  conv.Read (&mol, &molstream);

  if (mol.Empty ())
    return NULL;

  conv.Write (&mol, &V3000stream);

  outstring = V3000stream.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpV3000 = strdup (outstring.c_str ());

  return (tmpV3000);
}

extern "C" char *
ob_V3000_to_mol (char *V3000)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (V3000);
  string outstring;
  istringstream V3000stream (tmpStr);
  ostringstream molstream;
  char *tmpMolfile;

  conv.SetInAndOutFormats ("MDL", "MDL");
  conv.AddOption ("3", OBConversion::INOPTIONS);
  conv.AddOption ("2", OBConversion::OUTOPTIONS);

  conv.Read (&mol, &V3000stream);
  conv.Write (&mol, &molstream);

  outstring = molstream.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpMolfile = strdup (outstring.c_str ());

  return (tmpMolfile);
}

/*****************************************************************
 * This function Copyright (c) 2004
 * by Bayer Business Services GmbH
 *****************************************************************/
extern "C" double
ob_molweight (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  double molweight = 0.0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);

  molweight = mol.GetMolWt ();

  return (molweight);
}

extern "C" char *
ob_hillformula (char *smiles)
{

  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  string molfmla;
  OBMol mol;
  OBConversion conv;
  char *tmpFormula;

  conv.SetInAndOutFormats ("SMI", "SMI");
  conv.Read (&mol, &molstream);

  molfmla = mol.GetFormula ();

  tmpFormula = strdup (molfmla.c_str ());

  return (tmpFormula);
}

extern "C" unsigned int
ob_num_atoms (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  unsigned int numatoms = 0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);
  //mol.AddHydrogens (false, false);

  numatoms = mol.NumAtoms ();

  return (numatoms);
}

extern "C" unsigned int
ob_num_heavy_atoms (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  unsigned int numheavyatoms = 0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);

  numheavyatoms = mol.NumHvyAtoms ();

  return (numheavyatoms);
}

extern "C" unsigned int
ob_num_bonds (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  unsigned int numbonds = 0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);
  //mol.AddHydrogens (false, false);

  numbonds = mol.NumBonds ();

  return (numbonds);
}

extern "C" unsigned int
ob_num_rotatable_bonds (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  unsigned int numrotors = 0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);

  numrotors = mol.NumRotors ();

  return (numrotors);
}

extern "C" int
ob_total_charge (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  int totalcharge = 0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);

  totalcharge = mol.GetTotalCharge ();

  return (totalcharge);
}

extern "C" unsigned int
ob_is_nostruct (char *molfile)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  istringstream molstream (tmpStr);

  conv.SetInAndOutFormats ("MDL", "MDL");

  conv.Read (&mol, &molstream);

  return (mol.Empty ())? 1 : 0;
}

extern "C" unsigned int
ob_is_chiral (char *molfile)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  istringstream molstream (tmpStr);
  int chiral = 0;

  conv.SetInAndOutFormats ("MDL", "MDL");

  conv.Read (&mol, &molstream);
  mol.FindChiralCenters ();

  chiral = mol.IsChiral ()? 1 : 0;

  return (chiral);
}

extern "C" int
ob_2D (char *molfile)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  istringstream molstream (tmpStr);
  int twod = 0;

  conv.SetInAndOutFormats ("MDL", "MDL");

  conv.Read (&mol, &molstream);

  if (mol.Has2D ())
    twod++;

  return (twod);
}

extern "C" int
ob_3D (char *molfile)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  istringstream molstream (tmpStr);
  int threed = 0;

  conv.SetInAndOutFormats ("MDL", "MDL");

  conv.Read (&mol, &molstream);

  if (mol.Has3D ())
    threed++;

  return (threed);
}

extern "C" char *
ob_add_hydrogens (char *smiles, int polaronly, int correct4PH)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  string outstring;
  istringstream molstream1 (tmpStr);
  ostringstream molstream2;
  char *tmpMolfile;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream1);

  mol.AddHydrogens (polaronly != 0, correct4PH != 0);

  conv.Write (&mol, &molstream2);

  outstring = molstream2.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpMolfile = strdup (outstring.c_str ());

  return (tmpMolfile);
}

extern "C" char *
ob_delete_hydrogens (char *smiles, int nonpolaronly)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  string outstring;
  istringstream molstream1 (tmpStr);
  ostringstream molstream2;
  char *tmpMolfile;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream1);

if(mol.NumHvyAtoms () > 0) {
  if (nonpolaronly != 0)
    {
      mol.DeleteNonPolarHydrogens ();
    }
  else
    {
      mol.DeleteHydrogens ();
    }
} else {
    cout << "Warning: Cannot remove hydrogens. Resulting molecule would be empty!" << endl;
}        

  conv.Write (&mol, &molstream2);

  outstring = molstream2.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpMolfile = strdup (outstring.c_str ());

  return (tmpMolfile);
}

/* extern "C" char *
ob_strip_salts (char *molfile, int neutralize_residue)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  string outstring;
  istringstream molstream1 (tmpStr);
  ostringstream molstream2;
  char *tmpMolfile;

  conv.SetInAndOutFormats ("MDL", "MDL");

  conv.Read (&mol, &molstream1);

  mol.StripSalts (0);
  
  if(neutralize_residue != 0) {
      FOR_ATOMS_OF_MOL(a, mol) {
         a->SetFormalCharge(0);
         }
  }    

  conv.Write (&mol, &molstream2);

  outstring = molstream2.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpMolfile = strdup (outstring.c_str ());

  return (tmpMolfile);
} */

extern "C" char *
ob_strip_salts (char *smiles, int neutralize_residue)
{
  OBAtom atom;
  OBMol mol, largestFragment;
  OBConversion conv;
  string tmpStr (smiles);
  string outstring;
  istringstream molstream1 (tmpStr);
  ostringstream molstream2;
  vector<OBMol> fragments;
  vector <OBMol>::const_iterator i;
  char *tmpMolfile;
  int max = 0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream1);

  fragments = mol.Separate();
  
  for( i = fragments.begin(); i != fragments.end(); i++ ) {
      if (i->NumAtoms() > max) {
          max=i->NumAtoms();
          largestFragment = *i;
      }
  }       
 
  if(neutralize_residue != 0) {
      largestFragment.ConvertDativeBonds();
      FOR_ATOMS_OF_MOL(atom, largestFragment) {
         atom->SetFormalCharge(0);
         }
  }    

  conv.Write (&largestFragment, &molstream2);

  outstring = molstream2.str ();

  // remove the trailling $$$$\n from the SDFile
  if (outstring.find ("$$$$\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 5);
    }
  else if (outstring.find ("$$$$\r\n", 0) != string::npos)
    {
      outstring = outstring.substr (0, outstring.length () - 6);
    }

  tmpMolfile = strdup (outstring.c_str ());

  return (tmpMolfile);
}

/*****************************************************************
 * This function Copyright (c) 2005
 * by Bayer Business Services GmbH
 *****************************************************************/
extern "C" double
ob_exactmass (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  double exactmass = 0.0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);

  exactmass = mol.GetExactMass ();

  return (exactmass);
}

extern "C" void
ob_fp2 (char *molfile, unsigned int *fp)
{
  OBMol mol;
  OBConversion conv;
  OBFingerprint *fprint;
  string tmpStr (molfile);
  //string fp_mode ("FP2");
  istringstream molstream (tmpStr);
  vector < unsigned int >vfp;		 
  fprint = OBFingerprint::FindFingerprint ("FP2");
  vector < unsigned int >::iterator i;

  conv.SetInFormat ("MDL");

  conv.Read (&mol, &molstream);

  fprint->GetFingerprint (&mol, vfp);

  /*  for (i = vfp.begin();i != vfp.end();i++) {
     fp[j++]=*i;
     } */
     
  memset(fp,0x0, FPSIZE2*sizeof(unsigned int));

  memcpy (fp, &vfp[0], FPSIZE2*sizeof(unsigned int));

}

extern "C" void
ob_fp3_bin (char *serializedInput, unsigned int *fp)
{
  OBMol mol;
  vector < unsigned int >vfp;
  OBFingerprint *fprint = OBFingerprint::FindFingerprint ("FPPC8");
  //vector < unsigned int >::iterator i;
  int len;
  
   memset(fp,0x0, FPSIZE3*sizeof(unsigned int));    
  
  if (fprint != NULL) {
      
       ob_rehydrate_molecule(&mol, serializedInput);
       
        fprint->GetFingerprint (&mol, vfp);

   /* for (i = vfp.begin();i != vfp.end();i++) {
     fp[j++]=*i;
     printf("%i",*i);
     }  */
     
  len = vfp.size();
  
  //cout << len << endl;
  
  if (len > FPSIZE3) {
          len = FPSIZE3;
          cout << "Warning: Index dictionary size exceeded. Only the first " << FPSIZE3*sizeof(unsigned int)*8 << " bits will be used!" << endl; 
          }
          
     

  memcpy (fp, &vfp[0], len*sizeof(unsigned int));
  } else {
    cout << "FPPC8 fingerprint not found!" << endl;
}   

}

extern "C" void
ob_fp3 (char *molfile, unsigned int *fp)
{
  OBMol mol;
  OBConversion conv;
  OBFingerprint *fprint;
  string tmpStr (molfile);
  //string fp_mode ("evoFP");
  //string fp_mode ("FPPC8");
  istringstream molstream (tmpStr);
  vector < unsigned int >vfp;
  fprint = OBFingerprint::FindFingerprint ("FPPC8");
  //vector < unsigned int >::iterator i;
  int len;
  
  memset(fp,0x0, FPSIZE3*sizeof(unsigned int));    

if (fprint != NULL) {
    
  conv.SetInFormat ("MDL");

  conv.Read (&mol, &molstream);

  fprint->GetFingerprint (&mol, vfp);

   /* for (i = vfp.begin();i != vfp.end();i++) {
     fp[j++]=*i;
     printf("%i",*i);
     }  */
     
  len = vfp.size();
  
  //cout << len << endl;
  
  if (len > FPSIZE3) {
          len = FPSIZE3;
          cout << "Warning: Index dictionary size exceeded. Only the first " << FPSIZE3*sizeof(unsigned int)*8 << " bits will be used!" << endl; 
          }    

  memcpy (fp, &vfp[0], len*sizeof(unsigned int));
  } else {
    cout << "FPPC8 fingerprint not found!" << endl;
}   

}

extern "C" void
ob_fp_MACCS_bin (char *serializedInput, unsigned int *fp)
{
  OBMol mol;
  vector < unsigned int >vfp;
  OBFingerprint *fprint = NULL;
  
  fprint = OBFingerprint::FindFingerprint ("MACCS");
  //vector < unsigned int >::iterator i;
  int len;
  
  memset(fp,0x0, FPSIZE_MACCS*sizeof(unsigned int));    
  
  if (fprint != NULL) {
      
       ob_rehydrate_molecule(&mol, serializedInput);
       
            //cout << mol.NumHvyAtoms() << endl;
  
  //cout << fprint->GetID() << endl;
      
       fprint->GetFingerprint (&mol, vfp);

   /* for (i = vfp.begin();i != vfp.end();i++) {
     fp[j++]=*i;
     printf("%i",*i);
     }  */
     
  len = vfp.size();
  
  //cout << len << endl;
  
  if (len > FPSIZE_MACCS) {
          len = FPSIZE_MACCS;
          cout << "Warning: Index dictionary size exceeded. Only the first " << FPSIZE_MACCS*sizeof(unsigned int)*8 << " bits will be used!" << endl; 
          }
          
      

  memcpy (fp, &vfp[0], len*sizeof(unsigned int));
  } else {
    cout << "MACCS fingerprint not found!" << endl;
}   

}

extern "C" void
ob_fp_MACCS (char *molfile, unsigned int *fp)
{
  OBMol mol;
  OBConversion conv;
  OBFingerprint *fprint=NULL;
  string tmpStr (molfile);
  //string fp_mode ("evoFP");
  //string fp_mode ("FPPC8");
  istringstream molstream (tmpStr);
  vector < unsigned int >vfp;
  fprint = OBFingerprint::FindFingerprint ("MACCS");
  //vector < unsigned int >::iterator i;
  int len;
  
  memset(fp,0x0, FPSIZE_MACCS*sizeof(unsigned int));    

if (fprint != NULL) {
    
  conv.SetInFormat ("MDL");

  conv.Read (&mol, &molstream);
  
  //cout << mol.NumHvyAtoms() << endl;
  
  //cout << fprint->GetID() << endl;

  fprint->GetFingerprint (&mol, vfp);

   /* for (i = vfp.begin();i != vfp.end();i++) {
     fp[j++]=*i;
     printf("%i",*i);
     }  */
     
  len = vfp.size();
  
  //cout << len << endl;
  
  if (len > FPSIZE_MACCS) {
          len = FPSIZE_MACCS;
          cout << "Warning: Index dictionary size exceeded. Only the first " << FPSIZE_MACCS*sizeof(unsigned int)*8 << " bits will be used!" << endl; 
          }    

  memcpy (fp, &vfp[0], len*sizeof(unsigned int));
  } else {
    cout << "MACCS fingerprint not found!" << endl;
}   

}

extern "C" void
ob_fp_bin (char *serializedInput, unsigned int *fp) {
  OBMol mol;
  OBFingerprint *fprint2, *fprint3;
  vector < unsigned int >vfp; 
  fprint2 = OBFingerprint::FindFingerprint ("FP2");
  fprint3 = OBFingerprint::FindFingerprint ("FPPC8");
  vector < unsigned int >::iterator i;
  unsigned int *offset;
  int len;
  
  ob_rehydrate_molecule(&mol, serializedInput);

  fprint2->GetFingerprint (&mol, vfp);
  
  memset(fp,0x0, FPSIZE*sizeof(unsigned int));   

  memcpy (fp, &vfp[0], FPSIZE2*sizeof(unsigned int));
  
  if (fprint3 != NULL) {
  
  vfp.clear();
  
  fprint3->GetFingerprint (&mol, vfp);    
  
  len = vfp.size();  
  
  if (len > FPSIZE3) {
          len = FPSIZE3;
          cout << "Warning: Index dictionary size exceeded. Only the first " << FPSIZE3*sizeof(unsigned int)*8 << " bits will be used!" << endl; 
          }
          
  offset = fp+OFFSET;          
          
  memcpy (offset, &vfp[0], len*sizeof(unsigned int)); 
  } else {
    cout << "FPPC8 fingerprint not found!" << endl;
}   
}    

//extern "C" void
//ob_fp6 (char *molfile, unsigned int *fp)
//{
//  OBMol mol;
//  OBConversion conv;
//  OBFingerprint *fprint;
//  string tmpStr (molfile);
//  string fp_mode ("FPPC8");
//  istringstream molstream (tmpStr);
//  vector < unsigned int >vfp;
//  fprint = OBFingerprint::FindFingerprint (fp_mode);
//  //vector < unsigned int >::iterator i;
//  int len;
//
//  conv.SetInAndOutFormats ("MDL", "MDL");
//
//  conv.Read (&mol, &molstream);
//
//  fprint->GetFingerprint (&mol, vfp);
//
//   /* for (i = vfp.begin();i != vfp.end();i++) {
//     fp[j++]=*i;
//     printf("%i",*i);
//     }  */
//     
//  len = vfp.size();
//
//  memcpy (fp, &vfp[0], len*sizeof(unsigned int));
//
//}

//extern "C" void
//ob_fp7 (char *molfile, unsigned int *fp)
//{
//  OBMol mol;
//  OBConversion conv;
//  OBFingerprint *fprint;
//  string tmpStr (molfile);
//  string fp_mode ("FP7");
//  istringstream molstream (tmpStr);
//  vector < unsigned int >vfp;		 
//  fprint = OBFingerprint::FindFingerprint (fp_mode);
//  vector < unsigned int >::iterator i;
//  int j = 0;
//
//  conv.SetInAndOutFormats ("MDL", "MDL");
//
//  conv.Read (&mol, &molstream);
//
//  fprint->GetFingerprint (&mol, vfp);
//
//  /*  for (i = vfp.begin();i != vfp.end();i++) {
//     fp[j++]=*i;
//     } */
//
//  memcpy (fp, &vfp[0], FPSIZE2*sizeof(unsigned int));
//
//}

//extern "C" double
//ob_tanimoto_n (unsigned char *fp1, unsigned char *fp2)
//{
///* Fast tanimoto code with 8 bit LUT by Andrew Dalke as published in http://www.dalkescientific.com/writings/diary/archive/2008/06/27/computing_tanimoto_scores.html.
//   Used with permission */
//  
//  int and_count=0, or_count=0;
//  register int i=FPSIZE*sizeof(unsigned int);
//  
//  while (i--) {
//    //cout << i << endl;
//    or_count  += popcount_counts[*fp2 | *fp1];
//    and_count += popcount_counts[*fp2 & *fp1];
//    fp1++;
//    fp2++;
//  }
      
  /*i=FPSIZE3*sizeof(unsigned int);
  
  while (i--) {
    //cout << (*fp2 & 0x1) << endl;
    //if ((*fp2 != 0) || (*fp1 != 0)) or_count++;
    //if ((*fp2 != 0) && (*fp2 == *fp1)) and_count++;
    //if ((*fp2 | *fp1) & 0x1) or_count++;
    //if ((*fp2 & *fp1) & 0x1) and_count++;
    or_count  += popcount_counts[(*fp2 & 0x1) | (*fp1 & 0x1)];
    and_count += popcount_counts[(*fp2 & 0x1) & (*fp1 & 0x1)];
    //or_count  += popcount_counts[(*fp2 & 0x2) | (*fp1 & 0x2)];
    //and_count += popcount_counts[(*fp2 & 0x2) & (*fp1 & 0x2)];
    fp1++;
    fp2++;
  }*/

//  return static_cast<double>(and_count)/or_count;
//}

extern "C" double
ob_tanimoto (unsigned char *fp1, unsigned char *fp2, unsigned short size)
{
/* Fast tanimoto code with 8 bit LUT by Andrew Dalke as published in http://www.dalkescientific.com/writings/diary/archive/2008/06/27/computing_tanimoto_scores.html.
   Used with permission */
  
  unsigned int and_count=0, or_count=0;
  register int i=size;
  
  while (i--) {
      //cout << (int) (fp2[i] | fp1[i]) << endl;
    or_count  += popcount_counts[*fp2 | *fp1];
    and_count += popcount_counts[*fp2 & *fp1];
    fp1++;
    fp2++;
  }

  return static_cast<double>(and_count)/or_count;
}

extern "C" double
ob_tversky (unsigned char *fp1_prototype, unsigned char *fp2_variant, unsigned short size, double alpha_prototype, double beta_variant)
{  
  unsigned int a=0,b=0,c=0;
  register int i=size;
  
  if(alpha_prototype<0.0f) alpha_prototype=0.0f;
  if(alpha_prototype>1.0f) alpha_prototype=1.0f;
  if(beta_variant<0.0f) beta_variant=0.0f;
  if(beta_variant>1.0f) beta_variant=1.0f;
  
  while (i--) {
    a += popcount_counts[*fp1_prototype];
    b += popcount_counts[*fp2_variant];
    c += popcount_counts[*fp2_variant & *fp1_prototype];
    fp1_prototype++;
    fp2_variant++;
  }
  
  return static_cast<double>(c)/((a - c)*alpha_prototype+(b - c)*beta_variant+c);
}


extern "C" unsigned int
ob_popcount (unsigned char *fp, unsigned short size)
{  
  unsigned int popcount=0, i;
  
  for (i=0; i<size; i++) {
    popcount  += popcount_counts[fp[i]];
  }

  return popcount;
}

/* Copyright © The International Union of Pure and Applied Chemistry 2005: IUPAC
  International Chemical Identifier (InChI) (contact: secretariat@iupac.org) */

extern "C" char *
ob_smiles_to_inchi (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream inchistream;
  char *tmpInChI;

  conv.SetInAndOutFormats ("SMI", "INCHI");
  conv.AddOption ("w", OBConversion::OUTOPTIONS);

  conv.Read (&mol, &molstream);
  conv.Write (&mol, &inchistream);

  //cout << inchistream.str();

  outstring = inchistream.str ();

  outstring = outstring.substr (0, outstring.length () - 1);

  tmpInChI = strdup (outstring.c_str ());

  return (tmpInChI);
}

extern "C" char *
ob_smiles_to_inchikey (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream inchikeystream;
  char *tmpInChIkey;

  conv.SetInAndOutFormats ("SMI", "INCHI");
  conv.AddOption ("w", OBConversion::OUTOPTIONS);
  conv.AddOption ("K", OBConversion::OUTOPTIONS);

  conv.Read (&mol, &molstream);
  conv.Write (&mol, &inchikeystream);

  //cout << inchistream.str();

  outstring = inchikeystream.str ();

  outstring = outstring.substr (0, outstring.length () - 1);

  tmpInChIkey = strdup (outstring.c_str ());

  return (tmpInChIkey);
}

extern "C" char *
ob_molfile_to_inchikey (char *molfile)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (molfile);
  string outstring;
  istringstream molstream (tmpStr);
  ostringstream inchikeystream;
  char *tmpInChIkey;

  conv.SetInAndOutFormats ("MDL", "INCHI");
  conv.AddOption ("w", OBConversion::OUTOPTIONS);
  conv.AddOption ("K", OBConversion::OUTOPTIONS);

  conv.Read (&mol, &molstream);
  conv.Write (&mol, &inchikeystream);

  //cout << inchistream.str();

  outstring = inchikeystream.str ();

  outstring = outstring.substr (0, outstring.length () - 1);

  tmpInChIkey = strdup (outstring.c_str ());

  return (tmpInChIkey);
}

/*extern "C" unsigned int
ob_SSS_SMARTS (const char *smarts_pattern, char *molfile)
{
  OBMol mol;
  OBConversion conv;
  OBSmartsPattern sp;
  string tmpStr (molfile);
  istringstream molstream (tmpStr);
  vector<vector<int> > maplist;

  conv.SetInAndOutFormats ("MDL", "MDL");

  conv.Read (&mol, &molstream);
  
  if (mol.Empty()) return 0;
  
  sp.Init(smarts_pattern);
  
  //if(mol.NumHvyAtoms() < sp.NumAtoms()) return 0;
  
  sp.Match(mol);

  maplist = sp.GetUMapList();

  return maplist.size();
}*/

/*extern "C" int
ob_SSS_VF2 (char *molfile_sarg, char *molfile)
{
  OBMol mol, mol_sarg;
  bool match;
  VF2 vf2;
  OBConversion conv;
  string tmpStr (molfile), tmpStr_sarg(molfile_sarg);
  istringstream molstream (tmpStr), molstream_sarg(tmpStr_sarg);

  conv.SetInFormat ("MDL");

  conv.Read (&mol_sarg, &molstream_sarg);
  
  if (mol_sarg.Empty ())
    return 0;

  conv.Read (&mol, &molstream);
  
  if (mol.Empty ())
    return 0;

  match = vf2.Match (mol,mol_sarg);

  //maplist = sp.GetUMapList ();

  //return maplist.size ();
  
  return match ? 1 : 0;
}

extern "C" int
ob_SSS_VF2_bin (char *serializedInput_sarg, char *serializedInput)
{
  OBMol mol, mol_sarg;
  bool match, mol_ok;
  VF2 vf2;

   mol_ok = ob_rehydrate_molecule(&mol_sarg, serializedInput_sarg);
  
  if (!mol_ok || mol_sarg.Empty ())
    return 0;

  mol_ok = ob_rehydrate_molecule(&mol, serializedInput);
  
  if (!mol_ok || mol.Empty ())
    return 0;

  match = vf2.Match (mol,mol_sarg);

  //maplist = sp.GetUMapList ();

  //return maplist.size ();
  
  return match ? 1 : 0;
}

extern "C" int
ob_ESS_VF2_bin (char *serializedInput_sarg, char *serializedInput)
{
  OBMol mol, mol_sarg;
  bool match, mol_ok;
  VF2 vf2;

   mol_ok = ob_rehydrate_molecule(&mol_sarg, serializedInput_sarg);
  
  if (!mol_ok || mol_sarg.Empty ())
    return 0;

  mol_ok = ob_rehydrate_molecule(&mol, serializedInput);
  
  if (!mol_ok || mol.Empty ())
    return 0;

  match = vf2.ExactMatch (mol,mol_sarg);

  //maplist = sp.GetUMapList ();

  //return maplist.size ();
  
  return match ? 1 : 0;
}*/

extern "C" int
ob_SSS_SMARTS_native_bin (const char *smarts_pattern, char *serializedInput)
{
  OBMol mol;
  OBSmartsPattern sp;
  bool match, mol_ok;
 
  mol_ok = ob_rehydrate_molecule(&mol, serializedInput);

  if (!mol_ok || mol.Empty ())
    return 0;

  if (!sp.Init (smarts_pattern)) return -1;

  if(mol.NumHvyAtoms() < sp.NumAtoms()) return 0;

  match = sp.Match (mol,true);

  //maplist = sp.GetUMapList ();

  //return maplist.size ();
  
  return match ? 1 : 0;
}

extern "C" int
ob_SSS_SMARTS_native (const char *smarts_pattern, char *smiles)
{
  OBMol mol;
  OBSmartsPattern sp;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  bool match, mol_ok;

  conv.SetInFormat ("SMI");

  conv.Read (&mol, &molstream);

  if (mol.Empty ())
    return 0;

  if (!sp.Init (smarts_pattern)) return -1;

  if(mol.NumHvyAtoms() < sp.NumAtoms()) return 0;

  match = sp.Match (mol,true);

  //maplist = sp.GetUMapList ();

  //return maplist.size ();
  
  return match ? 1 : 0;
}

extern "C" int
ob_SSS_SMARTS_native_count_bin (const char *smarts_pattern, char *serializedInput)
{
  OBMol mol;
  OBSmartsPattern sp;
  unsigned int matchcount = 0;
  bool mol_ok;
  
  mol_ok = ob_rehydrate_molecule(&mol, serializedInput);

  if (!mol_ok || mol.Empty ())
    return matchcount;

  if(!sp.Init (smarts_pattern)) return -1;

  if(mol.NumHvyAtoms() < sp.NumAtoms()) return 0;

  if (sp.Match(mol,false)) {
      matchcount = sp.GetUMapList().size();
      //matchcount=sp.NumMatches();
  }    

  //maplist = sp.GetUMapList ();

  //return maplist.size ();
  
  return matchcount;
}

extern "C" int
ob_SSS_SMARTS_native_count (const char *smarts_pattern, char *smiles)
{
  OBMol mol;
  OBConversion conv;
  OBSmartsPattern sp;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  //vector <vector<int> > maplist;
  unsigned int matchcount = 0;
  

  conv.SetInFormat ("SMI");

  conv.Read (&mol, &molstream);

  if (mol.Empty ())
    return matchcount;

  if(!sp.Init (smarts_pattern)) return -1;

  if(mol.NumHvyAtoms() < sp.NumAtoms()) return 0;

  if (sp.Match(mol,false)) {
      matchcount = sp.GetUMapList().size();
      //matchcount=sp.NumMatches();
  }    

  //maplist = sp.GetUMapList ();

  //return maplist.size ();
  
  return matchcount;
}

extern "C" double
ob_MR (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  //OBMR mr; //2.1.1
  OBGroupContrib theMR("MR", "mr.txt", "molar refractivity"); 
  OBDescriptor* pDescr = OBDescriptor::FindType("MR");
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  double MR = 0.0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);
  mol.AddHydrogens (false, false);

  //MR = mr.Predict (mol);
   if(pDescr)
    MR = pDescr->Predict(&mol);
  

  return (MR);
}

extern "C" double
ob_PSA (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  //OBPSA psa; //2.1.1
  OBGroupContrib theTPSA("TPSA", "psa.txt", "topological polar surface area"); 
  OBDescriptor* pDescr = OBDescriptor::FindType("TPSA");
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  double PSA = 0.0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);
  mol.AddHydrogens (false, false);

  //PSA = psa.Predict (mol);
   if(pDescr)
    PSA = pDescr->Predict(&mol);

  return (PSA);
}

extern "C" double
ob_logP (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  //OBLogP logP; //2.1.1
  OBGroupContrib theLOGP("LOGP", "logp.txt", "log P"); 
  OBDescriptor* pDescr = OBDescriptor::FindType("LOGP");
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  double LOGP = 0.0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);
  mol.AddHydrogens (false, false);

  //LOGP = logP.Predict (mol);
  if(pDescr)
    LOGP = pDescr->Predict(&mol);

  return (LOGP);
}

extern "C" unsigned int
ob_num_H_donors (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  unsigned int numHdonors = 0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);
  
      FOR_ATOMS_OF_MOL(a, mol)
      {
            if(a->IsHbondDonor()) numHdonors++;
      }

  return (numHdonors);
}

extern "C" unsigned int
ob_num_H_acceptors (char *smiles)
{
  OBMol mol;
  OBConversion conv;
  string tmpStr (smiles);
  istringstream molstream (tmpStr);
  unsigned int numHacceptors = 0;

  conv.SetInAndOutFormats ("SMI", "SMI");

  conv.Read (&mol, &molstream);
  
      FOR_ATOMS_OF_MOL(a, mol)
      {
      if(a->IsHbondAcceptor()) numHacceptors++;
      }

  return (numHacceptors);
}

static bool ob_rehydrate_molecule(OBBase* pOb, char *serializedInput)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    map<OBAtom*,OBChiralData*> _mapcd;

    //Define some references so we can use the old parameter names
    OBMol &mol = *pmol;
    _mapcd.clear();
    bool chiralWatch=false;
  
    //stringstream errorMsg;
    //string clearError;

    // Allows addition of further disconnected atoms to an existing molecule
    //int offset = mol.NumAtoms(); 

    int i,natoms,nbonds;
    //string r1,r2;
    
    unsigned int *intptr = (unsigned int*) serializedInput;
    
    //cout << *intptr << endl;
     
    //mol.SetDimension(*intptr);
          
    //intptr++;

    natoms = *intptr;
    
    intptr++;
    
    nbonds = *intptr;
    
    intptr++;
    
    //cout << mol.GetDimension() << natoms << nbonds << endl;
    
    _ATOM *atomptr = (_ATOM*) intptr;

    //mol.BeginModify();
    
        mol.ReserveAtoms(natoms);

        OBAtom atom;
        int stereo;

        for (i = 1;i <= natoms;i++) {
            
        //cout << atomptr->atomicnum << " " << atomptr->isotope << " " << atomptr->idx << " " << atomptr->radical << " " << atomptr->stereo << " " << atomptr->formalcharge << endl;
         
         atom.SetIdx(atomptr->idx);
         
         atom.SetHyb(atomptr->hybridization);
         
          //int iso=0;
          atom.SetAtomicNum((int) atomptr->atomicnum);
          //iso=atomptr->isotope;
          
          //if(iso)
            atom.SetIsotope((unsigned int) atomptr->isotope);

          atom.SetFormalCharge((int) atomptr->formalcharge);

         stereo = atomptr->stereo;
         
              if (stereo == 2)
                {
                  chiralWatch=true;
                  atom.SetAntiClockwiseStereo();
                }
              else if (stereo == 1)
                {
                  chiralWatch=true;
                  atom.SetClockwiseStereo();
                }
              else if(stereo == 3)
                {
                  chiralWatch=true;
                  atom.SetChiral();
                }
                
          atom.SetSpinMultiplicity((short) atomptr->spinmultiplicity);
          
          if(atomptr->aromatic != 0) atom.SetAromatic();            
                
          if (!mol.AddAtom(atom)) return (false);
          
          if(chiralWatch)  // fill the map with data for each chiral atom
            _mapcd[mol.GetAtom(i)] = new OBChiralData;
          atom.Clear();
          atomptr++;
        }
        
        _BOND *bondptr = (_BOND*) atomptr;

        unsigned int start,end,order,flags;
        
        for (i = 0;i < nbonds;i++) {
          flags = 0;
         
          start = bondptr->beginidx;
          end = bondptr->endidx;
          order = (int) bondptr->order;
          
          if (start == 0 || end == 0 || order == 0 ||
              start > natoms || end > natoms)
            return false;

          order = (unsigned int) (order == 4) ? 5 : order;
         
          stereo = bondptr->stereo;
          
          //cout << atomptr->atomicnum << " " << atomptr->isotope << " " << atomptr->idx << " " << atomptr->radical << " " << atomptr->stereo << " " << atomptr->formalcharge << endl;
            
            if (stereo) {
              if (stereo == 1) flags |= OB_WEDGE_BOND;
              if (stereo == 6) flags |= OB_HASH_BOND;
            } 
            
          if(bondptr->aromatic != 0) flags |= OB_AROMATIC_BOND;
          
          //flags |= bondptr->flags;        

          if (!mol.AddBond(start,end,order,flags)) return (false);

          // after adding a bond to atom # "start+offset"
          // search to see if atom is bonded to a chiral atom
          // HERE
          map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
          ChiralSearch = _mapcd.find(mol.GetAtom(start));
          if (ChiralSearch!=_mapcd.end())
            {
              (ChiralSearch->second)->AddAtomRef(end, input);
            }
          // after adding a bond to atom # "end + offset"
          // search to see if atom is bonded to a chiral atom
          ChiralSearch = _mapcd.find(mol.GetAtom(end));
          if (ChiralSearch!=_mapcd.end())
            {
              (ChiralSearch->second)->AddAtomRef(start, input);
            }
         bondptr++;
        }
        
    //mol.AssignSpinMultiplicity();

    //mol.EndModify();
    
    intptr = (unsigned int*) bondptr;
        
    //NE add the OBChiralData stored inside the _mapcd to the atoms now after end
    // modify so they don't get lost.
    if(_mapcd.size()>0)
      {
        OBAtom* atomptr;
        OBChiralData* cd;
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
        for(ChiralSearch=_mapcd.begin();ChiralSearch!=_mapcd.end();ChiralSearch++)
          {
            atomptr=ChiralSearch->first;
            cd=ChiralSearch->second;
            atomptr->SetData(cd);
          }    
      }
      
    //mol.SetFlags(*(int*)intptr);  
      
    mol.SetAromaticPerceived();
    mol.SetKekulePerceived();
    mol.SetChiralityPerceived();  

    return(true);
  }
  
extern "C" char *ob_lyophilize_molecule(char* smiles) {
      OBMol mol;
      OBConversion conv;
      string tmpStr (smiles);
      istringstream molstream (tmpStr);
      conv.SetInFormat("SMI");
      conv.Read (&mol, &molstream);
      if(mol.Empty()) return NULL;
      unsigned int numatoms = mol.NumAtoms();
      unsigned int numbonds = mol.NumBonds();
      unsigned int totalsize=(numatoms*sizeof(_ATOM))+(numbonds*sizeof(_BOND))+(3*sizeof(unsigned int));
      char *retval = new char[totalsize];
      int  stereo;
      _ATOM *atomptr;
      _BOND *bondptr;
      
      //mol.Kekulize();
      
      memset(retval,0x0,totalsize);
      
      unsigned int *uintptr = (unsigned int*) retval;
      
      *uintptr = totalsize-sizeof(unsigned int);
      
      uintptr++;
      
      //*uintptr = (unsigned int) mol.GetDimension();
      
      //uintptr++;
      
      *uintptr = numatoms;
      
      uintptr++;
      
      *uintptr = numbonds;
      
      uintptr++;
      
      atomptr = (_ATOM*) uintptr;
      
      //cout << mol.GetDimension() << numatoms << numbonds << endl;
      
      
        FOR_ATOMS_OF_MOL(atom, mol) {
            stereo = 0;
            atomptr->idx = atom->GetIdx();
            atomptr->hybridization = atom->GetHyb();
            atomptr->atomicnum = (unsigned char) atom->GetAtomicNum();
            atomptr->formalcharge = (char) atom->GetFormalCharge();
            atomptr->isotope = (unsigned short) atom->GetIsotope();
            
            //cout << atom->GetFormalCharge() << " " << atom->GetIsotope() << " " << atom->GetSpinMultiplicity() << endl;
            
            if(atom->IsClockwise()) stereo = 1;
            else if (atom->IsAntiClockwise()) stereo = 2;
            else if (atom->IsChiral()) stereo = 3;
            atomptr->stereo = stereo;
            atomptr->spinmultiplicity = (unsigned char) atom->GetSpinMultiplicity();
            atomptr->aromatic = atom->IsAromatic() ? 1 : 0;
            //atomptr->flags = atom->GetFlag();
            atomptr++;
            }
        
        
       bondptr = (_BOND*) atomptr;    
       
       if(numbonds>0) { 
         FOR_BONDS_OF_MOL(bond, mol) {
             stereo = 0;
             bondptr->beginidx = bond->GetBeginAtomIdx();
             bondptr->endidx = bond->GetEndAtomIdx();
             bondptr->order = (unsigned char) bond->GetBondOrder();
             if(bond->IsWedge()) stereo = 1;
             else if(bond->IsHash()) stereo = 6;
             bondptr->stereo = stereo;
             bondptr->aromatic = bond->IsAromatic() ? 1 : 0;
             //bondptr->flags = bond->GetFlags();
             bondptr++;
             }
         }
         
      //int *intptr = (int*) bondptr;
      
      //*intptr = mol.GetFlags();           
      
      return retval;
  }    
