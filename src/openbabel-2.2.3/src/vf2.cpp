#include <openbabel/obiter.h>
#include <vf2/argraph.h>
#include <vf2/argedit.h>
#include <vf2/match.h>
#include <vf2/vf2_sub_state.h>
#include <vf2/vf2_state.h>
#include <vf2/vf2.h>

#define VERTEX_UNUSED -1

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
        
        
     if (_aa->GetAtomicNum() != _ab->GetAtomicNum()) return false;
      
      if (_aa->IsInRing()) {
          if (!_ab->IsInRing()) return false;
      }
      
      if (_aa->IsAromatic()) {
          if (!_ab->IsAromatic()) return false;
      }
      
      //if (_aa->GetIsotope() != _ab->GetIsotope()) return false;
      
      return true;
    }
};
  
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
      
      if (_ba->IsUp() != _bb->IsUp()) return false;
       
      if (_ba->IsDown() !=_bb->IsDown()) return false;
          
      if (_ba->IsHash() !=_bb->IsHash()) return false;
          
      if (_ba->IsWedge() != _bb->IsWedge()) return false;
      
      if (_ba->IsInRing()) {
          if (!_bb->IsInRing()) return false;
      }  
      
      if (_ba->IsAromatic()) {
          if (_bb->IsAromatic()) return true;
          else
          return false;
     }    
      
      if (_ba->GetBondOrder() != _bb->GetBondOrder()) return false;   
        
      return true;
    }
};

class _OBAtomComparatorEx: public AttrComparator
{ 
public:
  _OBAtomComparatorEx()
    { 
    }
  virtual bool compatible(void* aa, void* ab)
    {
      OBAtom* _aa = (OBAtom*)aa;
      OBAtom* _ab = (OBAtom*)ab;  
        
        
     if (_aa->GetAtomicNum() != _ab->GetAtomicNum()) return false;
     if (_aa->IsInRing() != _ab->IsInRing()) return false;
     if (_aa->GetHyb() != _ab->GetHyb()) return false;
     if (_aa->IsChiral() !=_ab->IsChiral()) return false;
     if (_aa->GetIsotope() != _ab->GetIsotope()) return false;
     if (_aa->GetFormalCharge() != _ab->GetFormalCharge()) return false;
     if (_aa->IsAromatic() !=_ab->IsAromatic()) return false;   
      
     return true;
    }
};
  
class _OBBondComparatorEx: public AttrComparator
{ 
public:
  _OBBondComparatorEx()
    {
    }
  virtual bool compatible(void* ba, void* bb)
    {
      OBBond* _ba = (OBBond*)ba;
      OBBond* _bb = (OBBond*)bb;
      bool _baa, _bba;
      
     if (_ba->IsUp() != _bb->IsUp()) return false;
       
      if (_ba->IsDown() !=_bb->IsDown()) return false;
          
      if (_ba->IsHash() !=_bb->IsHash()) return false;
          
      if (_ba->IsWedge() !=_bb->IsWedge()) return false;
      
     if (_ba->IsInRing() != _bb->IsInRing()) return false;
     
     _baa = _ba->IsAromatic();
     _bba = _bb->IsAromatic();
      
      if (_baa || _bba) {
          if (_baa == _bba) return true;
      else
          return false; 
      }    
      
      if (_ba->GetBondOrder() != _bb->GetBondOrder()) return false;    
         
      return true;
    }
};

    VF2::VF2() {};

    bool VF2::Match(OBMol &mol,OBMol &mol_sarg)
  {
      OBAtom *start_atom, *end_atom;
      ARGEdit ed_target, ed_sarg;//, ed_target2;
      int n,s,e;
      unsigned int i;
      
      if((mol.NumHvyAtoms() < mol_sarg.NumHvyAtoms()) || (mol.NumBonds() < mol_sarg.NumBonds())) return false;
      
      node_id ni_target[mol.NumAtoms()], ni_sarg[mol_sarg.NumAtoms()], start_idx, end_idx;
      
      vector<unsigned int> lt_target(mol.NumAtoms()+1,VERTEX_UNUSED);
      vector<unsigned int> lt_sarg(mol_sarg.NumAtoms()+1,VERTEX_UNUSED);    
      
      //map<int, node_id> node_id_map_target, node_id_map_sarg;
      
      //map<int, node_id>::iterator i;
      
      FOR_BONDS_OF_MOL(b, mol_sarg) {
          
         //cout << b->GetBondOrder() << endl;
           //cout << "A:" << b->IsAromatic() << endl;
        start_atom = b->GetBeginAtom();
        end_atom =  b->GetEndAtom();
        
        s = start_atom->GetIdx();
        e = end_atom->GetIdx();
        
       i = lt_sarg[s];
        
        if(i!=VERTEX_UNUSED) {
            start_idx = i;
         } else {
            start_idx=ed_sarg.InsertNode(start_atom);
            //node_id_map_sarg.insert(make_pair(s,start_idx));
            lt_sarg[s]=start_idx;
         }
         
         i = lt_sarg[e];  
         
         if(i!=VERTEX_UNUSED) {
              end_idx = i;
         } else {
             end_idx=ed_sarg.InsertNode(end_atom);
             //node_id_map_sarg.insert(make_pair(e,end_idx));
             lt_sarg[e]=end_idx;
         }            
                 
        ed_sarg.InsertEdge(start_idx,end_idx,(OBBond*)&*b);
        
        ed_sarg.InsertEdge(end_idx,start_idx,(OBBond*)&*b);
      }
  //}    
      
      FOR_BONDS_OF_MOL(b, mol) {
        start_atom = b->GetBeginAtom();
        end_atom =  b->GetEndAtom();
        
        s = start_atom->GetIdx();
        e = end_atom->GetIdx();
        
        i = lt_target[s];
        
        if(i!=VERTEX_UNUSED) {
            start_idx = i;
         } else {
             start_idx=ed_target.InsertNode(start_atom);
             //node_id_map_target.insert(make_pair(s,start_idx));  
             lt_target[s]=start_idx;
         }
         
         i = lt_target[e];
         
         if(i!=VERTEX_UNUSED) {
             end_idx = i;
         } else {
             end_idx=ed_target.InsertNode(end_atom);
             //node_id_map_target.insert(make_pair(e,end_idx));
             lt_target[e]=end_idx;
         }           
                 
        ed_target.InsertEdge(start_idx,end_idx,(OBBond*)&*b);

        ed_target.InsertEdge(end_idx,start_idx,(OBBond*)&*b);
      }
  //}    
      
      ARGraph<OBAtom*, OBBond*> g_target(&ed_target), g_sarg(&ed_sarg);
      
      g_sarg.SetNodeComparator(new _OBAtomComparator());  
      g_sarg.SetEdgeComparator(new _OBBondComparator());

      VF2SubState s0(&g_sarg, &g_target);
      if (match(&s0, &n, ni_sarg, ni_target)) return true;
      
      return false;
  } 
  
  bool VF2::ExactMatch(OBMol &mol,OBMol &mol_sarg)
  {
      OBAtom *start_atom, *end_atom;
      ARGEdit ed_target, ed_sarg;//, ed_target2;
      int n,s,e;
      unsigned int i;
      
      if((mol.NumHvyAtoms() != mol_sarg.NumHvyAtoms()) || (mol.NumBonds() != mol_sarg.NumBonds())) return false;
      
      node_id ni_target[mol.NumAtoms()], ni_sarg[mol_sarg.NumAtoms()], start_idx, end_idx;
      
      //map<int, node_id> node_id_map_target, node_id_map_sarg;
      
      //map<int, node_id>::iterator i;
      
      vector<unsigned int> lt_target(mol.NumAtoms()+1,VERTEX_UNUSED);
      vector<unsigned int> lt_sarg(mol_sarg.NumAtoms()+1,VERTEX_UNUSED);    
      
      FOR_BONDS_OF_MOL(b, mol_sarg) {
          
         //cout << b->GetBondOrder() << endl;
           //cout << "A:" << b->IsAromatic() << endl;
        start_atom = b->GetBeginAtom();
        end_atom =  b->GetEndAtom();
        
        s = start_atom->GetIdx();
        e = end_atom->GetIdx();
        
       i = lt_sarg[s];
        
        if(i!=VERTEX_UNUSED) {
            start_idx = i;
         } else {
            start_idx=ed_sarg.InsertNode(start_atom);
            //node_id_map_sarg.insert(make_pair(s,start_idx));
            lt_sarg[s]=start_idx;
         }
         
         i = lt_sarg[e];  
         
         if(i != VERTEX_UNUSED) {
              end_idx = i;
         } else {
             end_idx=ed_sarg.InsertNode(end_atom);
             //node_id_map_sarg.insert(make_pair(e,end_idx));
             lt_sarg[e]=end_idx;
         }              
                 
        ed_sarg.InsertEdge(start_idx,end_idx,(OBBond*)&*b);
       
        ed_sarg.InsertEdge(end_idx,start_idx,(OBBond*)&*b);
      }
  //}    
      
      FOR_BONDS_OF_MOL(b, mol) {
        start_atom = b->GetBeginAtom();
        end_atom =  b->GetEndAtom();
        
        s = start_atom->GetIdx();
        e = end_atom->GetIdx();
        
       i = lt_target[s];
        
        if(i != VERTEX_UNUSED) {
            start_idx = i;
         } else {
             start_idx=ed_target.InsertNode(start_atom);
             //node_id_map_target.insert(make_pair(s,start_idx));  
             lt_target[s]=start_idx;
         }
         
         i = lt_target[e];
         
         if(i != VERTEX_UNUSED) {
             end_idx = i;
         } else {
             end_idx=ed_target.InsertNode(end_atom);
             //node_id_map_target.insert(make_pair(e,end_idx));
             lt_target[e]=end_idx;
         }        
                 
        ed_target.InsertEdge(start_idx,end_idx,(OBBond*)&*b);
      
        ed_target.InsertEdge(end_idx,start_idx,(OBBond*)&*b);
      }
  //}    
      
      ARGraph<OBAtom*, OBBond*> g_target(&ed_target), g_sarg(&ed_sarg);
      
      g_sarg.SetNodeComparator(new _OBAtomComparatorEx());  
      g_sarg.SetEdgeComparator(new _OBBondComparatorEx());

      VF2State s0(&g_sarg, &g_target);
      if (match(&s0, &n, ni_sarg, ni_target)) return true;
      
      return false;
  }     
}  
  
