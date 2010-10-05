#include <map>
#include <openbabel/obiter.h>
#include <vf2/argraph.h>
#include <vf2/argedit.h>
#include <vf2/match.h>
#include <vf2/vf2_sub_state.h>
#include <vf2/vf_sub_state.h>
#include <vf2/ull_sub_state.h>
#include <vf2/vf2.h>

using namespace std;

namespace OpenBabel
{

//class OBAtomDestroyer: public AttrDestroyer
//  { public:
//       virtual void destroy(void *a)
//         { delete ((OBAtom*)a);
//         }
//  };
  
class _OBAtomComparator: public AttrComparator
{ 
public:
  _OBAtomComparator()
    { 
    }
  virtual bool compatible(void* aa, void* ab)
    {
      OBAtom* _aa = (OBAtom*)aa;
      OBAtom* _ab = (OBAtom*)ab;  
        
      if (_aa->IsInRing() && !_ab->IsInRing()) return false;
      if (_aa->IsAromatic() != _ab->IsAromatic()) return false;  
      if (_aa->GetAtomicNum() != _ab->GetAtomicNum()) return false;
      //if (_aa->GetFormalCharge() != _ab->GetFormalCharge()) return false;
      if (_aa->GetIsotope() != _ab->GetIsotope()) return false;
      //if (_aa->GetHyb() != _ab->GetHyb()) return false;
      //if (_aa->IsChiral() != _ab->IsChiral()) return false;
      
      return true;
    }
};
  
//class OBBondDestroyer: public AttrDestroyer
//  { public:
//       virtual void destroy(void *b)
//         { delete ((OBBond*) b);
//         }
//  };
  
class _OBBondComparator: public AttrComparator
{ 
public:
  _OBBondComparator()
    {
    }
  virtual bool compatible(void* ba, void* bb)
    {
      OBBond* _ba = (OBBond*)ba;
      OBBond* _bb = (OBBond*)bb;
      
      if (_ba->IsInRing() && !_bb->IsInRing()) return false;    
      if (_ba->IsAromatic() == _bb->IsAromatic()) return true;  
      if (_ba->GetBondOrder() != _bb->GetBondOrder()) return false;
      if (_ba->IsUp() != _bb->IsUp()) return false;
      if (_ba->IsDown() != _bb->IsDown()) return false;
      if (_ba->IsHash() != _bb->IsHash()) return false;
      if (_ba->IsWedge() != _bb->IsWedge()) return false;
      
      return true;
    }
};

    VF2::VF2() {};

    bool VF2::Match(OBMol mol,OBMol mol_sarg)
  {
      OBAtom *start_atom, *end_atom;
      ARGEdit ed_target, ed_sarg;//, ed_target2;
      int n,s,e;
      
      //if(mol.NumHvyAtoms() < mol_sarg.NumHvyAtoms()) return false;
      
      node_id ni_target[mol.NumAtoms()], ni_sarg[mol_sarg.NumAtoms()], start_idx, end_idx;
      
      map<int, node_id> node_id_map_target, node_id_map_sarg;
      
      map<int, node_id>::iterator i;
      
      /*if(mol_sarg.NumBonds() == 0) {*/
          
      /*FOR_ATOMS_OF_MOL(a, mol_sarg) {
          node_id_map_sarg.insert(make_pair(a->GetIdx(),ed_sarg.InsertNode((OBAtom*)&*a)));
          //node_id_map_sarg2.insert(make_pair(a->GetIdx(),ed_sarg2.InsertNode((OBAtom*)&*a)));
      } */   
    /*  }
  } else {  */ 
  
   /* FOR_ATOMS_OF_MOL(a, mol) {
        if(node_id_map_target.find(a->GetAtomicNum()) == node_id_map_target.end()) node_id_map_target.insert(make_pair(a->GetAtomicNum(),0));
    }
    
        FOR_ATOMS_OF_MOL(a, mol_sarg) {
        if (node_id_map_target.find(a->GetAtomicNum()) == node_id_map_target.end()) return false;
    }
    
    node_id_map_target.clear(); */           
      
      FOR_BONDS_OF_MOL(b, mol_sarg) {
          
         //cout << b->GetBondOrder() << endl;
           //cout << "A:" << b->IsAromatic() << endl;
        start_atom = b->GetBeginAtom();
        end_atom =  b->GetEndAtom();
        
        s = start_atom->GetIdx();
        e = end_atom->GetIdx();
        
       i = node_id_map_sarg.find(s);
        
        if(i != node_id_map_sarg.end()) {
            start_idx = i->second;
         } else {
            start_idx=ed_sarg.InsertNode(start_atom);
            node_id_map_sarg.insert(make_pair(s,start_idx));
         }
         
         i = node_id_map_sarg.find(e);  
         
         if(i != node_id_map_sarg.end()) {
              end_idx = i->second;
         } else {
             end_idx=ed_sarg.InsertNode(end_atom);
             node_id_map_sarg.insert(make_pair(e,end_idx));
         }            
                 
        //start_idx = node_id_map_sarg[b->GetBeginAtom()->GetIdx()];
        //end_idx = node_id_map_sarg[b->GetEndAtom()->GetIdx()];
        ed_sarg.InsertEdge(start_idx,end_idx,(OBBond*)&*b);
        //start_idx = node_id_map_sarg2[b->GetBeginAtom()->GetIdx()];
        //end_idx = node_id_map_sarg2[b->GetEndAtom()->GetIdx()];
        ed_sarg.InsertEdge(end_idx,start_idx,(OBBond*)&*b);
      }
  //}    
         
  /*if(mol.NumBonds() == 0) {*/
                   
      /*FOR_ATOMS_OF_MOL(a, mol) {
          node_id_map_target.insert(make_pair(a->GetIdx(),ed_target.InsertNode((OBAtom*)&*a)));
          //node_id_map_target2.insert(make_pair(a->GetIdx(),ed_target2.InsertNode((OBAtom*)&*a)));
      }*/
      
  /*} else { */   
      
      FOR_BONDS_OF_MOL(b, mol) {
        start_atom = b->GetBeginAtom();
        end_atom =  b->GetEndAtom();
        
        s = start_atom->GetIdx();
        e = end_atom->GetIdx();
        
        //already_set = (node_id_map_target.find(start_atom->GetIdx()) != node_id_map_target.end());
        
        i = node_id_map_target.find(s);
        
        if(i != node_id_map_target.end()) {
            start_idx = i->second;
         } else {
             start_idx=ed_target.InsertNode(start_atom);
             node_id_map_target.insert(make_pair(s,start_idx));  
         }
         
         i = node_id_map_target.find(e);
         
         if(i != node_id_map_target.end()) {
             end_idx = i->second;
         } else {
             end_idx=ed_target.InsertNode(end_atom);
             node_id_map_target.insert(make_pair(e,end_idx));
         }           
                 
        //start_idx = node_id_map_sarg[b->GetBeginAtom()->GetIdx()];
        //end_idx = node_id_map_sarg[b->GetEndAtom()->GetIdx()];
        ed_target.InsertEdge(start_idx,end_idx,(OBBond*)&*b);
        //start_idx = node_id_map_sarg2[b->GetBeginAtom()->GetIdx()];
        //end_idx = node_id_map_sarg2[b->GetEndAtom()->GetIdx()];
        ed_target.InsertEdge(end_idx,start_idx,(OBBond*)&*b);
      }
  //}    
      
      ARGraph<OBAtom*, OBBond*> g_target(&ed_target), g_sarg(&ed_sarg);
      
      g_sarg.SetNodeComparator(new _OBAtomComparator());  
      g_sarg.SetEdgeComparator(new _OBBondComparator());
      
      //g_sarg2.SetNodeComparator(new _OBAtomComparator());  
      //g_sarg2.SetEdgeComparator(new _OBBondComparator());   

      VF2SubState s0(&g_sarg, &g_target);
      if (match(&s0, &n, ni_sarg, ni_target)) return true;
      //VF2SubState s1(&g_sarg2, &g_target);
      //if (match(&s1, &n, ni_sarg, ni_target)) return true;   
//      VF2SubState s2(&g_sarg, &g_target2);
//      if (match(&s2, &n, ni_sarg, ni_target)) return true; 
//      VF2SubState s3(&g_sarg2, &g_target2);
//      if (match(&s3, &n, ni_sarg, ni_target)) return true; 
      
      return false;
  }     
}  
  
