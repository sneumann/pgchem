#include <mol.h>
#include <obconversion.h>
#include <iostream>
#include "argraph.h"
#include "argedit.h"
#include "match.h"
#include "vf2_sub_state.h"

#if HAVE_SSTREAM
#include <sstream>
#elif HAVE_SSTREAM_H
#include <sstream.h>
#endif

using namespace OpenBabel;
using namespace std;

class OBAtomDestroyer: public AttrDestroyer
  { public:
       virtual void destroy(void *a)
         { delete ((OBAtom*)a);
         }
  };
  
class OBAtomComparator: public AttrComparator
{ 
public:
  OBAtomComparator()
    { 
    }
  virtual bool compatible(void* aa, void* ab)
    {
      return (((OBAtom*)aa)->GetAtomicNum() == ((OBAtom*)ab)->GetAtomicNum());
    }
};
  
class OBBondDestroyer: public AttrDestroyer
  { public:
       virtual void destroy(void *b)
         { delete ((OBBond*) b);
         }
  };
  
  class OBBondComparator: public AttrComparator
{ 
public:
  OBBondComparator()
    {
    }
  virtual bool compatible(void* ba, void* bb)
    {
      return (((OBBond*)ba)->GetBO() == ((OBBond*)bb)->GetBO());
    }
};

int main()
{
  OBMol mol1, mol2;
  OBConversion conv;
  string tmpStr1 ("c1ccccc1");
  string tmpStr2 ("c1ccccc1F");
  istringstream smistream1(tmpStr1);
  istringstream smistream2(tmpStr2);
  int i,j,n;
  
   conv.SetInAndOutFormats ("SMI", "MDL");
   
   conv.Read (&mol1, &smistream1);
   
   conv.Read (&mol2, &smistream2);
   
   node_id ni1[mol1.NumAtoms()], ni2[mol2.NumAtoms()]; 
    
   return 0;
}    
