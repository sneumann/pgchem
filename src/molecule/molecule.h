#include "postgres.h"
#include "molecule_limits.h"

#ifndef SET_VARSIZE // < 8.3 compatibility
#define SET_VARSIZE(v,l) (VARATT_SIZEP(v) = (l))
#endif

typedef struct
{
  uint32 fp[FPSIZE];
} MOLFP;

typedef struct
{
  int4 len;
  int4 sizemf;
  int4 sizesmi;
  int4 disconnected;
  uint32 fp[FPSIZE];
  char inchikey[INCHIKEYSZ];
  char data[1];
} MOLECULE;

//#define MOLHDRSZ		             (4*sizeof(int4))
#define SMIPTR(x)		 	         ((char*)(x->data))
#define MFPTR(x)		 		     ((char*)(x->data+x->sizesmi))
#define ANCPTR(x)		 		     ((char*)(x->data+x->sizesmi+x->sizemf))
//#define CALCDATASZ(sizemf, sizesmi)  (((FPSIZE)*sizeof(uint32)) + MOLHDRSZ + INCHIKEYSZ + sizemf + sizesmi)
#define CALCDATASZ(sizemf, sizesmi, sizeancnfo)  (sizeof(MOLECULE) + sizemf + sizesmi + sizeancnfo - sizeof(char))
//#define CALCDATASZ(sizemf, sizesmi, sizeefa)  ((sizeof(MOLECULE)) + sizemf + sizesmi + (sizeefa*sizeof(unsigned int)) - sizeof(char))

#define PG_GETARG_MOLFP_P(n)         (MOLFP *) DatumGetPointer(PG_GETARG_DATUM(n))
#define PG_RETURN_MOLFP_P(n)         PG_RETURN_POINTER(n)

#define DatumGetMoleculeP(n)         ((MOLECULE *) PG_DETOAST_DATUM(n))
#define PG_GETARG_MOLECULE_P(n)      DatumGetMoleculeP(PG_GETARG_DATUM(n))
#define PG_RETURN_MOLECULE_P(n)      PG_RETURN_POINTER(n)

//MOLECULE *new_molecule (char *smiles, char *molfile, unsigned int *efa_array);
MOLECULE *new_molecule (char *smiles, char *molfile);
