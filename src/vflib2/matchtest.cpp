#include "argraph.h"
#include "argedit.h"
#include "match.h"
#include "vf2_sub_state.h"
#include <iostream>

using namespace std;

#define MAXNODES 8

class _Atom
  { public:
      unsigned short atomic_num;
      _Atom(unsigned short atomic_num)
        { this->atomic_num=atomic_num;
        }
  };
  
class _AtomDestroyer: public AttrDestroyer
  { public:
       virtual void destroy(void *a)
         { delete ((_Atom*)a);
         }
  };
  
class _AtomComparator: public AttrComparator
{ 
public:
  _AtomComparator()
    { 
    }
  virtual bool compatible(void* aa, void* ab)
    {
      return (((_Atom*)aa)->atomic_num == ((_Atom*)ab)->atomic_num);
    }
};

class _Bond
  { public:
      unsigned char order;
      _Bond(unsigned char order)
        { this->order=order;
        }
  };
  
class _BondDestroyer: public AttrDestroyer
  { public:
       virtual void destroy(void *b)
         { delete ((_Bond*) b);
         }
  };
  
class _BondComparator: public AttrComparator
{ 
public:
  _BondComparator()
    {
    }
  virtual bool compatible(void* ba, void* bb)
    {
      return (((_Bond*)ba)->order == ((_Bond*)bb)->order);
    }
};
  
int main()
  { ARGEdit small_ed, large_ed;
  node_id ni1[MAXNODES], ni2[MAXNODES];  // The object used to create the graph
    int i,j,n;

        // Insert the four nodes
        for(i=0; i<4; i++)
          small_ed.InsertNode(new _Atom(6)); // The inserted node will have index i.
                               // NULL stands for no semantic attribute.
                               
        for(i=0; i<8; i++)
          large_ed.InsertNode(new _Atom(6));                       

        // Insert the edges
        for(i=0; i<4; i++)
          for(j=0; j<4; j++)
            if (i!=j)
                  small_ed.InsertEdge(i, j, new _Bond(2)); // NULL stands for no sem. attribute.
                  
        for(i=0; i<8; i++)
          for(j=0; j<8; j++)
            if (i!=j)
                  large_ed.InsertEdge(i, j, new _Bond(2)); // NULL stands for no sem. attribute.          

        
        // Now the Graph can be constructed...
        ARGraph<_Atom, void> small_graph(&small_ed), large_graph(&large_ed);
        
        small_graph.SetNodeDestroyer(new _AtomDestroyer());
        large_graph.SetNodeDestroyer(new _AtomDestroyer());
        
        small_graph.SetEdgeDestroyer(new _BondDestroyer());
        large_graph.SetEdgeDestroyer(new _BondDestroyer());
        
        small_graph.SetNodeComparator(new _AtomComparator());
        
        small_graph.SetEdgeComparator(new _BondComparator());

        VF2SubState s0(&small_graph, &large_graph);

  if (!match(&s0, &n, ni1, ni2))
    { cout << "No match found!" << endl;
      return -1;
    }

  cout << "Found a match on " << n << " Atoms:" << endl;
  
  for(i=0; i<n; i++)
    cout << "\tAtom " << ni1[i] << " of molecule 1 is paired with Atom " << ni2[i] << " of molecule 2" << endl;

  return 0;
}
