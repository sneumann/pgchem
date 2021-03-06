#include <mol.h>
#include <obconversion.h>


using namespace std;
using namespace OpenBabel;

int main() {  
   
     OBAtom a, b, c;
      a.SetAtomicNum(6);
      b.SetAtomicNum(8);
      c.SetAtomicNum(8);

     OBMol mol;
     mol.AddAtom(a);
     mol.AddAtom(b);
     mol.AddAtom(c);
     
     mol.AddBond(1,2,1);
     mol.AddBond(2,3,1);

      OBConversion conv;
      conv.SetOutFormat("SMI");
      cout << conv.WriteString(&mol,1) << endl;
      
     OBSmartsPattern sp;
     
     sp.Init ("O~*");
     
     sp.Match (mol);
     
     cout << sp.NumMatches() << endl;
       
     cout << sp.GetUMapList().size() << endl;
      
      return EXIT_SUCCESS;
  }    
