/**********************************************************************
 * functions.sql native chemistry handling functions SQL stubs
 *
 * Copyright (c) 2004,2010 by Ernst-G. Schmid
 *
 * This file is part of the pgchem::tigress project.
 * For more information, see
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 ************************************************************************/

CREATE OR REPLACE FUNCTION add_hydrogens(molecule, boolean, boolean)
  RETURNS molecule AS
'libpgchem', 'pgchem_add_hydrogens'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION exactmass(molecule)
  RETURNS double precision AS
'libpgchem', 'pgchem_exactmass'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION fgroup_codes(molecule)
  RETURNS text AS
'libpgchem', 'pgchem_fgroup_codes_a'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION is_2d(molecule)
  RETURNS boolean AS
'libpgchem', 'pgchem_2D'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION is_3d(molecule)
  RETURNS boolean AS
'libpgchem', 'pgchem_3D'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION is_chiral(molecule)
  RETURNS boolean AS
'libpgchem', 'pgchem_is_chiral'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION is_nostruct(molecule)
  RETURNS boolean AS
'libpgchem', 'pgchem_is_nostruct'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION inchi(molecule)
  RETURNS text AS
'libpgchem', 'pgchem_molecule_to_inchi'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION canonical_smiles(molecule, boolean)
  RETURNS text AS
'libpgchem', 'pgchem_molecule_to_canonical_smiles'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION smiles(molecule, boolean)
  RETURNS text AS
'libpgchem', 'pgchem_molecule_to_smiles'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION v3000(molecule)
  RETURNS text AS
'libpgchem', 'pgchem_molecule_to_V3000'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molfile(molecule)
  RETURNS text AS
'libpgchem', 'pgchem_molecule_to_molfile'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molformula(molecule)
  RETURNS text AS
'libpgchem', 'pgchem_hillformula'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molweight(molecule)
  RETURNS double precision AS
'libpgchem', 'pgchem_molweight'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molkeys_long(molecule, boolean, boolean, boolean)
  RETURNS text AS
'libpgchem', 'pgchem_ms_fingerprint_long_a'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molkeys_long(molecule)
  RETURNS text AS
$BODY$
BEGIN
RETURN molkeys_long($1,false,false,false);
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molkeys_short(molecule)
  RETURNS text AS
'libpgchem', 'pgchem_ms_fingerprint_short_a'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_atoms(molecule)
  RETURNS integer AS
'libpgchem', 'pgchem_num_atoms'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_bonds(molecule)
  RETURNS integer AS
'libpgchem', 'pgchem_num_bonds'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_heavy_atoms(molecule)
  RETURNS integer AS
'libpgchem', 'pgchem_num_heavy_atoms'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_rotatable_bonds(molecule)
  RETURNS integer AS
'libpgchem', 'pgchem_num_rotatable_bonds'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION pgchem_barsoi_version()
  RETURNS cstring AS
'libpgchem', 'pgchem_barsoi_version'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION pgchem_version()
  RETURNS cstring AS
'libpgchem', 'pgchem_version'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION remove_hydrogens(molecule, boolean)
  RETURNS molecule AS
'libpgchem', 'pgchem_remove_hydrogens'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION strip_salts(molecule, boolean)
  RETURNS molecule AS
'libpgchem', 'pgchem_strip_salts'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION total_charge(molecule)
  RETURNS integer AS
'libpgchem', 'pgchem_total_charge'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION validate_cas_no(character varying)
  RETURNS boolean AS
$BODY$
DECLARE checksum_from_cas_no varchar;
DECLARE cas_no_left varchar;
DECLARE cas_no_right varchar;
DECLARE cas_no_full varchar;
DECLARE tmpsum int;
DECLARE position_multiplier int;
DECLARE caslen int;
BEGIN
caslen:=length($1);

IF caslen<5 OR caslen>12 THEN RETURN FALSE;
END IF;

checksum_from_cas_no:=split_part($1,'-',3)::int;
cas_no_left:=split_part($1,'-',1);
cas_no_right:=split_part($1,'-',2);
cas_no_full:=cas_no_left || cas_no_right;

if(length(cas_no_left)>7 OR length(cas_no_right)>2 OR length(checksum_from_cas_no)!=1) THEN 
return false; 
END IF;

caslen:=length(cas_no_full);
tmpsum:=0;
position_multiplier:=1;

 FOR i IN REVERSE caslen..1 LOOP
  tmpsum:=tmpsum+substr(cas_no_full,i,1)::int*position_multiplier;
  position_multiplier:=position_multiplier+1;
 END LOOP;
 RETURN tmpsum % 10 = checksum_from_cas_no::int;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION logP(molecule)
  RETURNS double precision AS
'libpgchem', 'pgchem_logP'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION MR(molecule)
  RETURNS double precision AS
'libpgchem', 'pgchem_MR'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION PSA(molecule)
  RETURNS double precision AS
'libpgchem', 'pgchem_PSA'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION SMARTSmatch(text, molecule)
  RETURNS boolean AS
'libpgchem', 'pgchem_smartsfilter'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION SMARTSmatch_count(text, molecule)
  RETURNS integer AS
'libpgchem', 'pgchem_smartsfilter_count'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION inchikey(molecule)
  RETURNS text AS
'libpgchem', 'pgchem_molecule_to_inchikey'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION fpstring(molecule)
  RETURNS bit varying AS
'libpgchem', 'pgchem_fp_out'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_h_acceptors(molecule)
  RETURNS integer AS
'libpgchem', 'pgchem_num_H_acceptors'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_h_donors(molecule)
  RETURNS integer AS
'libpgchem', 'pgchem_num_H_donors'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION lipinsky(mol molecule)
  RETURNS text AS
$BODY$
DECLARE parameters text;
BEGIN
parameters := '';

IF number_of_H_donors(mol) > 5 THEN
parameters := parameters || 'A';
END IF;
IF molweight(mol) > 500 THEN
parameters := parameters || 'B';
END IF;
IF logP(mol) > 5.0 THEN
parameters := parameters || 'C';
END IF;
IF number_of_H_acceptors(mol) > 10 THEN
parameters := parameters || 'D';
END IF;

RETURN parameters;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION mutabor(molecule)
  RETURNS molecule AS
'libpgchem', 'pgchem_mutate_fp'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION delebor(molecule)
  RETURNS molecule AS
'libpgchem', 'pgchem_blank_fp'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION fp2string(struct molecule)
  RETURNS bit varying AS
$BODY$
DECLARE fp bit varying;
BEGIN
fp := substring(fpstring(struct) for 1024);
RETURN fp;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION fp3string(struct molecule)
  RETURNS bit varying AS
$BODY$
DECLARE fp bit varying;
BEGIN
fp := substring(fpstring(struct) from 1025);
RETURN fp;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;

/*CREATE OR REPLACE FUNCTION nse(tablename text, columnname text)
  RETURNS SETOF int AS
$BODY$
DECLARE fpAllOn bit varying;
DECLARE fpAllOff bit varying;
DECLARE row record;
DECLARE i int;
BEGIN
fpAllOn := B'11111111111111111111111111111111';
fpAllOff := B'00000000000000000000000000000000';

FOR row IN EXECUTE 'SELECT fp3string(' || columnname || ') as fp from ' || tablename LOOP
fpAllOn = fpAllOn & row.fp;
fpAllOff = fpAllOff | row.fp;
END LOOP;

FOR i IN 1..32 LOOP

IF (substring(fpAllOn from i for 1)='1') THEN RETURN NEXT i; END IF;
IF (substring(fpAllOff from i for 1)='0') THEN RETURN NEXT i; END IF;

END LOOP;

RETURN;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;*/

CREATE OR REPLACE FUNCTION nbits_set(bit varying)
  RETURNS integer AS
'libpgchem', 'pgchem_nbits_set'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION disconnected(molecule)
  RETURNS boolean AS
'libpgchem', 'pgchem_disconnected'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION check_fingerprint_optimization(t_schema text, t_name text, s_column text)
  RETURNS double precision AS
$BODY$
DECLARE i float8;
DECLARE tablesize integer;
DECLARE t_q_name text;
DECLARE t_b_name text;
BEGIN

t_q_name := t_schema || '.q_c_' || t_name;
t_b_name := t_schema || '.' || t_name;

EXECUTE 'DROP TABLE IF EXISTS '||t_q_name;

EXECUTE 'CREATE TABLE ' || t_q_name || ' AS SELECT ' || s_column || ' FROM ' || t_b_name || ' WHERE fp2string('||s_column||') IN (SELECT fp2string('||s_column||') FROM ' || t_b_name || ' GROUP BY fp2string('||s_column||') HAVING (COUNT(fp2string('||s_column||'))>1))';

EXECUTE 'UPDATE '||t_q_name||' SET structure=mutabor('||s_column||')';

EXECUTE 'SELECT COUNT(1) FROM (SELECT count(fpstring('||s_column||')) as count FROM '||t_q_name||' GROUP BY fpstring('||s_column||')) as t WHERE t.count = 1' into i;

EXECUTE 'SELECT COUNT(1) FROM '||t_q_name into tablesize;

EXECUTE 'DROP TABLE IF EXISTS '||t_q_name;

RETURN i/tablesize;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION create_random_sample_q_table(t_schema text, t_name text, s_column text)
  RETURNS void AS
$BODY$
DECLARE n numeric;
DECLARE i integer;
DECLARE largeN integer;
DECLARE ts_size integer;
DECLARE t numeric;
DECLARE p numeric;
DECLARE q numeric;
DECLARE d numeric;
DECLARE numerator numeric;
DECLARE denominator numeric;
DECLARE t_q_name text;
DECLARE t_s_q_name text;
BEGIN

t:=1.96;
p:=0.5;
q:=0.5;
d:=0.05;

t_q_name := t_schema || '.q_' || t_name; 
t_s_q_name := t_schema || '.s_q_' || t_name; 

EXECUTE 'DROP TABLE IF EXISTS '||t_s_q_name;

EXECUTE 'CREATE TABLE '||t_s_q_name||' ('||s_column||' molecule) WITH (OIDS=FALSE)';

EXECUTE 'SELECT count(1) from '||t_q_name into largeN;

numerator:=t^2*(p*q);

denominator:=d^2;

n:=round((numerator/denominator),0);

--raise info 'n=%',n;

IF (n/largeN >= 0.05) THEN
n := round(((numerator/denominator)/((((numerator/denominator)-1.0)/largeN)+1.0))::numeric,0);
END IF;

--raise info 'n=%',n;

FOR i IN 1..n LOOP

EXECUTE 'INSERT INTO '||t_s_q_name||' (SELECT structure from '||t_q_name||' WHERE fp2string('||s_column||')=(SELECT fp2string('||s_column||') FROM '||t_q_name||' ORDER BY RANDOM() LIMIT 1))';

END LOOP;

EXECUTE 'SELECT count(1) from '||t_s_q_name into ts_size;

--raise info 'ts_size=%',ts_size;

IF (ts_size>=largeN) THEN
BEGIN
RAISE INFO 'Sampling has no effect. Retaining original table.';
EXECUTE 'DROP TABLE IF EXISTS '||t_s_q_name;
END;
ELSE
BEGIN
EXECUTE 'DROP TABLE IF EXISTS '||t_q_name;
EXECUTE 'ALTER TABLE '||t_s_q_name||' RENAME TO q_' || t_name;
END;
END IF;

RETURN;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION optimize_fingerprint(t_schema text, t_name text, s_column text, algorithm text, use_sampling boolean, basedict text, p_limit integer)
  RETURNS double precision AS
$BODY$
BEGIN

IF (algorithm='LP') THEN 
BEGIN
PERFORM optimize_fingerprint_step_one(t_schema, t_name, s_column, use_sampling , basedict ,p_limit);
RETURN optimize_fingerprint_step_two(t_schema, t_name, s_column);
END;
ELSIF (algorithm='GA') THEN
RAISE EXCEPTION 'Genetic optimization only externally supported'; 
ELSE RAISE EXCEPTION 'Algorithm % not supported',algorithm; 
END IF;

RAISE INFO 'Optimization on % complete. Achieved rate %',t_schema||'.'||t_name, rate;

RETURN 0.0;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION optimize_fingerprint_step_one(t_schema text, t_name text, s_column text, use_sampling boolean, basedict text, p_limit integer)
  RETURNS void AS
$BODY$
DECLARE expansionarray smallint[9];
DECLARE curr_row record;
DECLARE curr_row_i record;
DECLARE tablesize integer;
DECLARE hitcount integer;
DECLARE matchcount integer;
DECLARE expcount integer;
DECLARE j integer;
DECLARE e float8;
DECLARE x float8;
DECLARE t_dict_tmp_name text;
DECLARE t_dict_name text;
DECLARE t_sel_name text;
DECLARE t_q_name text;
DECLARE t_b_name text;
DECLARE t_h_name text;
DECLARE arraytext text;
BEGIN

t_dict_tmp_name := t_schema || '.' || t_name || '_dictionary_tmp';
t_dict_name := t_schema || '.' || t_name || '_dictionary';
t_sel_name := t_schema || '.' || t_name || '_selectivity';
t_q_name := t_schema || '.q_' || t_name; 
t_b_name := t_schema || '.' || t_name; 
t_h_name := t_schema || '.' || t_name || '_helper'; 

EXECUTE 'DROP TABLE IF EXISTS ' || t_dict_tmp_name;

EXECUTE 'DROP TABLE IF EXISTS ' || t_dict_name;

EXECUTE 'DROP TABLE IF EXISTS ' || t_sel_name;

EXECUTE 'CREATE TABLE ' || t_dict_tmp_name || ' (id serial NOT NULL, "SMARTS" text NOT NULL) WITH (OIDS=FALSE)';

EXECUTE 'COPY ' || t_dict_tmp_name || ' ("SMARTS") FROM ' ||  quote_literal(basedict);

EXECUTE 'CREATE TABLE ' || t_dict_name || ' (id serial NOT NULL, "SMARTS" text NOT NULL) WITH (OIDS=FALSE)';

EXECUTE 'INSERT INTO ' || t_dict_name ||  ' ("SMARTS") (SELECT distinct("SMARTS") FROM ' || t_dict_tmp_name || ')';

EXECUTE 'DROP TABLE IF EXISTS ' || t_dict_tmp_name;

EXECUTE 'CREATE TABLE ' || t_sel_name || '(pattern_id integer NOT NULL,expansion smallint NOT NULL DEFAULT 0,coverage double precision NOT NULL DEFAULT 0.0,weight double precision,exparray integer[],_e double precision,_x double precision,CONSTRAINT pattern_id_unique UNIQUE (pattern_id)) WITH (OIDS=FALSE)';

EXECUTE 'DROP TABLE IF EXISTS ' || t_q_name;

EXECUTE 'CREATE TABLE ' || t_q_name || ' AS SELECT ' || s_column || ' FROM ' || t_b_name || ' WHERE fp2string('||s_column||') IN (SELECT fp2string('||s_column||') FROM ' || t_b_name || ' GROUP BY fp2string('||s_column||') HAVING (COUNT(fp2string('||s_column||'))>1))';

IF (use_sampling) THEN
	PERFORM create_random_sample_q_table(t_schema, t_name, s_column);
END IF;

EXECUTE 'DELETE FROM '|| t_q_name ||' WHERE inchikey(' || s_column || ') IN ((SELECT inchikey(' || s_column || ') FROM '|| t_q_name ||' GROUP BY inchikey(' || s_column || ') HAVING (COUNT(inchikey(' || s_column || '))>1)))';

EXECUTE 'SELECT count(1) FROM '|| t_q_name into tablesize;

e := tablesize / 9;

FOR curr_row IN EXECUTE 'SELECT id, "SMARTS" FROM ' || t_dict_name LOOP

x := 0.0;

expansionarray := '{0,0,0,0,0,0,0,0,0}';

hitcount := 0;

expcount := 0;

FOR curr_row_i IN EXECUTE 'SELECT smartsmatch_count('||quote_literal(curr_row."SMARTS")||', ' || s_column || ') as mc FROM '|| t_q_name ||' WHERE smartsmatch_count('||quote_literal(curr_row."SMARTS")||', ' || s_column || ')>0' LOOP

	hitcount:=hitcount+1;

	matchcount := curr_row_i.mc;

	IF (matchcount > 8) THEN matchcount := 8; END IF;

	expansionarray[matchcount+1] := expansionarray[matchcount+1] + 1;

END LOOP;

IF (hitcount > 0) THEN

expansionarray[1] = tablesize - hitcount;

FOR i IN 1..9 LOOP

	IF (expansionarray[i] > 0) THEN 
		BEGIN
		expcount := expcount + 1;
		x := x+(power(expansionarray[i] - e,2)/e);
		END;
	END IF;

END LOOP;

arraytext :='{';

FOR i IN 1..9 LOOP
	arraytext := arraytext || expansionarray[i] || ',';
END LOOP;

arraytext := substring(arraytext FROM 1 for length(arraytext)-1) || '}';

EXECUTE 'INSERT INTO '||t_sel_name||' (pattern_id, expansion, coverage,exparray,_e,_x) VALUES ('||curr_row.id||','||expcount||','||hitcount/tablesize::float8||','||quote_literal(arraytext)||','||e||','||x||')';

END IF;

END LOOP;

EXECUTE 'UPDATE '||t_sel_name||' SET weight = _x';

EXECUTE 'CREATE TABLE '||t_h_name||' AS SELECT "SMARTS", expansion from '||t_sel_name||','|| t_dict_name||'  where id=pattern_id and pattern_id IN (SELECT pattern_id FROM '||t_sel_name||' ORDER BY weight ASC LIMIT '||p_limit||')' ;

EXECUTE 'COPY (SELECT ''#Comments after SMARTS'' UNION SELECT * from (SELECT "SMARTS" FROM '||t_h_name||') t) TO '||quote_literal('c:/tigress/obdata/dictionary.txt');

EXECUTE 'DROP TABLE IF EXISTS '||t_h_name;

EXECUTE 'DROP TABLE IF EXISTS '||t_dict_name;

EXECUTE 'DROP TABLE IF EXISTS '||t_sel_name;

END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION optimize_fingerprint_step_two(t_schema text, t_name text, s_column text)
  RETURNS double precision AS
$BODY$
DECLARE i float8;
DECLARE tablesize integer;
DECLARE t_q_name text;
BEGIN

t_q_name := t_schema || '.q_' || t_name; 

EXECUTE 'UPDATE '||t_q_name||' SET structure=mutabor('||s_column||')';

EXECUTE 'SELECT COUNT(1) FROM (SELECT count(fpstring('||s_column||')) as count FROM '||t_q_name||' GROUP BY fpstring('||s_column||')) as t WHERE t.count = 1' into i;

EXECUTE 'SELECT COUNT(1) FROM '||t_q_name into tablesize;

EXECUTE 'DROP TABLE IF EXISTS '||t_q_name;

RETURN i/tablesize;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;


