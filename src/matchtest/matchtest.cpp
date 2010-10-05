#include <mol.h>
#include <obconversion.h>


using namespace std;
using namespace OpenBabel;

int main() {  
   
     OBAtom a, b, c;
      a.SetAtomicNum(8);
      b.SetAtomicNum(6);
      c.SetAtomicNum(8);

     OBMol mol;
     mol.AddAtom(a);
     mol.AddAtom(b);
     mol.AddAtom(c);
     
     mol.AddBond(1,2,2);
     mol.AddBond(2,3,2);

      OBConversion conv;
      conv.SetOutFormat("SMI");
      cout << conv.WriteString(&mol,1) << endl;
      
     OBSmartsPattern sp;
     
     sp.Init ("C~*");
     
     sp.Match (mol,false);
     
       cout << sp.NumMatches() << endl;
       
        cout << sp.GetUMapList().size() << endl;
      
      return EXIT_SUCCESS;
  }    
