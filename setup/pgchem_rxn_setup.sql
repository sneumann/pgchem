DROP TYPE reaction CASCADE;
DROP TYPE rxnfp CASCADE;

CREATE TYPE reaction;
CREATE TYPE rxnfp;

CREATE FUNCTION rxnfp_in(cstring)
    RETURNS rxnfp
    AS 'libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION rxnfp_out(rxnfp)
    RETURNS cstring
   AS 'libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE rxnfp (
   input = rxnfp_in,
   output = rxnfp_out,
   internallength = 256,
   storage = PLAIN
   );

CREATE FUNCTION reaction_in(cstring)
    RETURNS reaction
    AS 'libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION reaction_out(reaction)
    RETURNS cstring
   AS 'libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION reaction_recv(internal)
   RETURNS reaction
  AS 'libpgchem'
   LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION reaction_send(reaction)
   RETURNS bytea
  AS 'libpgchem'
   LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE reaction (
   input = reaction_in,
   output = reaction_out,
   --receive = reaction_recv,
   --send = reaction_send,
   internallength = VARIABLE,
   storage = EXTENDED
   );

CREATE FUNCTION rxnfp_compress(internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION rxnfp_decompress(internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION rxnfp_penalty(internal,internal,internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C' WITH (isstrict);

CREATE FUNCTION rxnfp_picksplit(internal,internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C' WITH (isstrict);

CREATE FUNCTION rxnfp_union(internal, internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION rxnfp_same(internal, internal, internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION rxnfp_consistent(internal,internal,int4)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION reaction_contained_in(reaction,reaction)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION reaction_contains(reaction,reaction)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION reaction_equals(reaction,reaction)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION reaction_equals_products_exact(reaction,reaction)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION reaction_equals_exact(reaction,reaction)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION reaction_similarity(reaction,reaction)
RETURNS double precision
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION reaction_similarity_reactants(reaction,reaction)
RETURNS double precision
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION reaction_similarity_products(reaction,reaction)
RETURNS double precision
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE OPERATOR <= (
		 LEFTARG = reaction,
		 RIGHTARG = reaction,
		 PROCEDURE = reaction_contained_in,
                 COMMUTATOR = >=,
		 RESTRICT = contsel,
		 JOIN = contjoinsel
);

CREATE OPERATOR >= (
		 LEFTARG = reaction,
		 RIGHTARG = reaction,
		 PROCEDURE = reaction_contains,
                 COMMUTATOR = <=,
		 RESTRICT = contsel,
		 JOIN = contjoinsel
);

CREATE OPERATOR = (
		 LEFTARG = reaction,
		 RIGHTARG = reaction,
		 PROCEDURE = reaction_equals,
		 COMMUTATOR = =,
		 RESTRICT = eqsel,
		 JOIN = eqjoinsel
);

CREATE OPERATOR == (
		 LEFTARG = reaction,
		 RIGHTARG = reaction,
		 PROCEDURE = reaction_equals_products_exact,
		 COMMUTATOR = ==,
		 RESTRICT = eqsel,
		 JOIN = eqjoinsel
);

CREATE OPERATOR === (
		 LEFTARG = reaction,
		 RIGHTARG = reaction,
		 PROCEDURE = reaction_equals_exact,
		 COMMUTATOR = ===,
		 RESTRICT = eqsel,
		 JOIN = eqjoinsel
);

CREATE OPERATOR @ (
    leftarg = reaction,
    rightarg = reaction,
    procedure = reaction_similarity
);

CREATE OPERATOR <@ (
    leftarg = reaction,
    rightarg = reaction,
    procedure = reaction_similarity_reactants
);

CREATE OPERATOR @> (
    leftarg = reaction,
    rightarg = reaction,
    procedure = reaction_similarity_products
);

CREATE OPERATOR CLASS gist_reaction_ops
DEFAULT FOR TYPE reaction USING gist
AS
        OPERATOR        3       =        ,
        OPERATOR        5       ==       ,
        OPERATOR        6       ===      ,
        OPERATOR        7       >=       ,
        OPERATOR        8       <=       ,
        FUNCTION        1       rxnfp_consistent (internal, internal, int4),
        FUNCTION        2       rxnfp_union (internal, internal),
        FUNCTION        3       rxnfp_compress (internal),
        FUNCTION        4       rxnfp_decompress (internal),
        FUNCTION        5       rxnfp_penalty (internal, internal, internal),
        FUNCTION        6       rxnfp_picksplit (internal, internal),
        FUNCTION        7       rxnfp_same (internal, internal, internal),
        STORAGE		rxnfp;


CREATE OR REPLACE FUNCTION reaction_in(text)
  RETURNS reaction AS
'libpgchem', 'reaction_in_text'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION reaction_in(character varying)
  RETURNS reaction AS
'libpgchem', 'reaction_in_varchar'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION reaction_in(bytea)
  RETURNS reaction AS
'libpgchem', 'reaction_in_bytea'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE CAST (text AS reaction)
  WITH FUNCTION reaction_in(text);

CREATE CAST (character varying AS reaction)
  WITH FUNCTION reaction_in(character varying);

CREATE CAST (bytea AS reaction)
  WITH FUNCTION reaction_in(bytea);

CREATE OR REPLACE FUNCTION number_of_products(reaction)
  RETURNS integer AS
'libpgchem', 'pgchem_r_num_products'
  LANGUAGE 'c' IMMUTABLE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION number_of_reactants(reaction)
  RETURNS integer AS
'libpgchem', 'pgchem_r_num_reactants'
  LANGUAGE 'c' IMMUTABLE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION reaction_molecule(reaction, integer)
  RETURNS molecule AS
'libpgchem', 'pgchem_r_molecule_at'
  LANGUAGE 'c' IMMUTABLE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION fpstring(reaction)
  RETURNS bit varying AS
'libpgchem', 'pgchem_r_fp_out'
  LANGUAGE 'c' IMMUTABLE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION reactant(rxn reaction, pos integer)
  RETURNS molecule AS
$BODY$
BEGIN
IF (pos <1 OR pos>number_of_reactants(rxn)) THEN
	RAISE EXCEPTION 'Reactant index out of bounds: %', pos;
END IF;
RETURN reaction_molecule(rxn,pos);
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT
  COST 100;

CREATE OR REPLACE FUNCTION product(rxn reaction, pos integer)
  RETURNS molecule AS
$BODY$
BEGIN
IF (pos<1 OR pos>number_of_products(rxn)) THEN
	RAISE EXCEPTION 'Product index out of bounds: %', pos;
END IF;
RETURN reaction_molecule(rxn,pos+number_of_reactants(rxn));
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT
  COST 100;

CREATE OR REPLACE FUNCTION smiles(reaction)
  RETURNS text AS
'libpgchem', 'pgchem_r_reaction_to_smiles'
  LANGUAGE 'c' IMMUTABLE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION strip_rxninfo(molecule)
  RETURNS molecule AS
'libpgchem', 'pgchem_reaction_mol_strip_rxninfo'
  LANGUAGE 'c' IMMUTABLE STRICT
  COST 1;  

CREATE OR REPLACE FUNCTION reaction.t_decompose_reaction()
  RETURNS trigger AS
$BODY$
DECLARE current_mol molecule;
DECLARE num_molecules integer;
DECLARE current_position integer;
DECLARE current_id integer;
BEGIN

IF TG_OP='INSERT' OR TG_OP='UPDATE' THEN

  -- Reactants

  num_molecules := number_of_reactants(NEW.rxn);

  FOR current_position IN 1..num_molecules LOOP
	current_mol := strip_rxninfo(reactant(NEW.rxn, current_position));
	
	SELECT id INTO current_id FROM rm WHERE mol = current_mol;

	IF current_id is NULL THEN
		INSERT INTO rm (mol) VALUES (current_mol);
		current_id:=currval('rm_id_seq');
	END IF;

	INSERT INTO rmx (rid,mid,"position","role") VALUES (NEW.id, current_id, current_position, 'R');
	
  END LOOP;

  -- Products

  num_molecules := number_of_products(NEW.rxn);

  FOR current_position IN 1..num_molecules LOOP
	current_mol := strip_rxninfo(product(NEW.rxn, current_position));
	
	SELECT id INTO current_id FROM rm WHERE mol = current_mol;

	IF current_id is NULL THEN
		INSERT INTO rm (mol) VALUES (current_mol);
		current_id:=currval('rm_id_seq');
	END IF;

	INSERT INTO rmx (rid,mid,"position","role") VALUES (NEW.id, current_id, current_position, 'P');
  END LOOP;

ELSE

  RAISE EXCEPTION 'PGCHEM DECOMPOSE REACTION TRIGGER CALLED OUTSIDE INSERT OR UPDATE!';

END IF;
RETURN NULL;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE
  COST 500;


