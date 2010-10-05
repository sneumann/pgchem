-- Function: public.decode_rxnfile(bytea)

-- DROP FUNCTION public.decode_rxnfile(bytea);

CREATE OR REPLACE FUNCTION public.decode_rxnfile(bytea)
  RETURNS SETOF public.reactionmap AS
$BODY$
DECLARE workfile text;
DECLARE current_line text;
DECLARE current_molfile text;
DECLARE current_fragment int4;
DECLARE rm_row reactionmap;
DECLARE num_educts int4;
DECLARE num_products int4;
DECLARE max_fragments int4;
DECLARE e_pos int4;
DECLARE p_pos int4;
DECLARE abs_pos int4;
DECLARE current_molecule bytea;
BEGIN
num_educts:=0;
num_products:=0;
e_pos:=0;
p_pos:=0;
abs_pos:=0;
workfile:=molecule_to_molfile($1);

FOR i IN 1..5 LOOP
 current_line:=split_part(workfile,'\n',i);
END LOOP;

num_educts:=substring(current_line from 1 for 3)::int2;
num_products:=substring(current_line from 4 for 3)::int2;

max_fragments:=(num_educts+num_products)+1;

	FOR current_fragment IN 2..max_fragments LOOP

	current_molfile:=split_part(workfile,'$MOL\n',current_fragment);

	current_molecule:=tweak_molecule_a(molfile_to_molecule(current_molfile),true);

	abs_pos:=abs_pos+1;

	IF abs_pos>num_educts THEN
		p_pos:=p_pos+1;
		rm_row:=ROW(2,abs_pos,p_pos,current_molecule,num_educts,num_products,num_educts+num_products,ms_fingerprint_long_a(current_molecule));
	ELSE
		e_pos:=e_pos+1;
		rm_row:=ROW(1,abs_pos,e_pos,current_molecule,num_educts,num_products,num_educts+num_products,ms_fingerprint_long_a(current_molecule));
	END IF;
		RETURN NEXT rm_row;
	END LOOP;

	RETURN;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;
ALTER FUNCTION public.decode_rxnfile(bytea) OWNER TO postgres;



-- Function: public.reaction_hash(bytea)

-- DROP FUNCTION public.reaction_hash(bytea);

CREATE OR REPLACE FUNCTION public.reaction_hash(bytea)
  RETURNS text AS
$BODY$
DECLARE molstatistics_combined_hash_lhs text;
DECLARE molstatistics_combined_hash_rhs text;
DECLARE current_row record;
BEGIN

molstatistics_combined_hash_lhs:='';
molstatistics_combined_hash_rhs:='';

FOR current_row IN select type ,md5(ms_fingerprint_long) as hash from decode_rxnfile($1) ORDER BY type,hash
LOOP

IF current_row.type=1 THEN
	molstatistics_combined_hash_lhs:=molstatistics_combined_hash_lhs || '+' || current_row.hash;
ELSIF current_row.type=2 THEN
	molstatistics_combined_hash_rhs:=molstatistics_combined_hash_rhs || '+' || current_row.hash;
END IF;

END LOOP;

molstatistics_combined_hash_lhs:=btrim(molstatistics_combined_hash_lhs,'+');
molstatistics_combined_hash_rhs:=btrim(molstatistics_combined_hash_rhs,'+');

molstatistics_combined_hash_lhs:=molstatistics_combined_hash_lhs || '->' || molstatistics_combined_hash_rhs;


RETURN molstatistics_combined_hash_lhs;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;
ALTER FUNCTION public.reaction_hash(bytea) OWNER TO postgres;


-- Function: reaction_match_exact(bytea, int4, bool, bool, bool)

-- DROP FUNCTION reaction_match_exact(bytea, int4, bool, bool, bool);

CREATE OR REPLACE FUNCTION reaction_match_exact(bytea, int4, bool, bool, bool)
  RETURNS bool AS
$BODY$
DECLARE query_row record;
DECLARE result_row record;
DECLARE educt_pos integer;
DECLARE product_pos integer;
BEGIN

educt_pos:=0;
product_pos:=0;

FOR query_row IN SELECT type,fragment,md5(ms_fingerprint_long) as hash,n_e,n_p FROM decode_rxnfile($1) ORDER BY type,hash
LOOP

IF query_row.type=1 THEN

educt_pos:=educt_pos+1;
--is_match:=match_exact_a(query_row.fragment,(select molecule from r_fragments where rxn_id=$2 and fragment_type=1 and ordered_fragment_pos=educt_pos),false,false,false);
select molecule,n_e from r_fragments,r_rf_crossref where r_rf_crossref.fragment_id = r_fragments.iiid and rxn_id=$2 and fragment_type=1 and ordered_fragment_pos=educt_pos and n_e=query_row.n_e and n_p=query_row.n_p into result_row;

--IF (query_row.n_e!=result_row.n_e OR query_row.n_p!=result_row.n_p) THEN
--RETURN FALSE;
IF (result_row.n_e IS NULL OR match_exact_a(query_row.fragment,result_row.molecule,$3,$4,$5)=false) THEN
RETURN FALSE;
END IF;

ELSIF query_row.type=2 THEN

product_pos:=product_pos+1;
select molecule,n_p from r_fragments,r_rf_crossref where r_rf_crossref.fragment_id = r_fragments.iiid and rxn_id=$2 and fragment_type=2 and ordered_fragment_pos=product_pos and n_e=query_row.n_e and n_p=query_row.n_p into result_row;

--IF (query_row.n_e!=result_row.n_e OR query_row.n_p!=result_row.n_p) THEN
--RETURN FALSE;
IF (result_row.n_p IS NULL OR match_exact_a(query_row.fragment,result_row.molecule,$3,$4,$5)=false) THEN
RETURN FALSE;
END IF;

END IF;

END LOOP;

RETURN true;

END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;
ALTER FUNCTION reaction_match_exact(bytea, int4, bool, bool, bool) OWNER TO chembank;


-- Function: public.reaction_match_exact(bytea, bytea)

-- DROP FUNCTION public.reaction_match_exact(bytea, bytea);

CREATE OR REPLACE FUNCTION public.reaction_match_exact(bytea, bytea,bool,bool,bool)
  RETURNS bool AS
$BODY$
DECLARE query_row record;
DECLARE result_row record;
DECLARE rcursor CURSOR FOR SELECT type,fragment,md5(ms_fingerprint_long) as hash,n_e,n_p FROM decode_rxnfile($1) ORDER BY type,hash FOR READ ONLY;
BEGIN

OPEN rcursor;

FOR query_row IN SELECT type,fragment,md5(ms_fingerprint_long) as hash,n_e,n_p FROM decode_rxnfile($2) ORDER BY type,hash
LOOP

FETCH rcursor INTO result_row;

IF (query_row.n_e!=result_row.n_e OR query_row.n_p!=result_row.n_p) THEN
CLOSE rcursor;
RETURN FALSE;
ELSIF (match_exact_a(query_row.fragment,result_row.fragment,$3,$4,$5)=false) THEN
CLOSE rcursor;
RETURN FALSE;
END IF;

END LOOP;

CLOSE rcursor;

RETURN true;

END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;


-- Function: public.validate_reaction(bytea)

-- DROP FUNCTION public.validate_reaction(bytea);

CREATE OR REPLACE FUNCTION public.validate_reaction(bytea)
  RETURNS bool AS
$BODY$
DECLARE query_row record;
BEGIN

FOR query_row IN SELECT n_e,n_p FROM decode_rxnfile($1)
LOOP

IF (query_row.n_e=0 OR query_row.n_p=0) THEN return false; END IF;
EXIT;

END LOOP;

RETURN true;

END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;


-- Function: public.reaction_match_substruct(bytea, bytea)

-- DROP FUNCTION public.reaction_match_substruct(bytea, bytea);

CREATE OR REPLACE FUNCTION public.reaction_match_substruct(bytea, bytea,bool,bool,bool)
  RETURNS bool AS
$BODY$
DECLARE query_row record;
DECLARE result_row record;
DECLARE update_row record;
DECLARE include text;
DECLARE max_cols text;
DECLARE csum int4;
DECLARE i int4;
DECLARE tid int4;
DECLARE total int4;
--DECLARE must_normalize bool;
BEGIN

tid:=nextval('r_matchmatrix_tid_seq');
--must_normalize:=false;
max_cols:='0';

FOR query_row IN SELECT type,fragment,n_e,n_p,abs_pos FROM decode_rxnfile($1)
LOOP

csum:=0;
PERFORM fastmatch_set_substructure_query_a(query_row.fragment,$3,$4,$5);

FOR result_row IN SELECT fragment,abs_pos FROM decode_rxnfile($2) WHERE query_row.n_e<=n_e and query_row.n_p<=n_p and type=query_row.type
LOOP
IF (fastmatch_a(result_row.fragment)=true) THEN
	IF (SELECT true FROM r_matchmatrix WHERE target_pos=result_row.abs_pos) THEN
		EXECUTE 'UPDATE r_matchmatrix set query_'||query_row.abs_pos||'=1 where target_pos='||result_row.abs_pos;
	ELSE
		EXECUTE 'INSERT INTO r_matchmatrix (query_'||query_row.abs_pos||',target_pos,tid) VALUES (1,'||result_row.abs_pos||','||tid||')';
	END IF;
	csum:=csum+1;
END IF;
END LOOP;

--RAISE NOTICE '%',csum;
IF (csum=0) THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid; return false; 
--ELSIF (csum>1) THEN must_normalize:=true;
END IF;

END LOOP;

--IF (must_normalize=false) THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid; return true; END IF;

total:=query_row.n_e+query_row.n_p;

WHILE (true)
LOOP

--retval:=false;
include:='';

FOR i IN 1..total LOOP

max_cols:=max_cols||'+query_'||i;

FOR result_row IN EXECUTE 'SELECT sum(query_' || i ||') as col_sum,count(target_pos) as row_count FROM r_matchmatrix where tid='||tid
LOOP
	IF (result_row.row_count<total OR result_row.col_sum=0) THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid; return false;
	ELSIF (result_row.col_sum>1) THEN include:=include || 'query_' || i || ',';
	END IF;
END LOOP;

--RAISE NOTICE '%',include;

--IF (retval=false) THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid; 
--return false; END IF;

END LOOP;

IF (include='') THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid;
return true; END IF;

--RAISE NOTICE '%',max_cols;

FOR result_row IN EXECUTE 'SELECT target_pos FROM r_matchmatrix where '|| max_cols ||'>1 and tid='||tid
LOOP
		EXECUTE 'UPDATE r_matchmatrix set '|| split_part(include,',',1) ||'=0 WHERE target_pos='||result_row.target_pos||' and tid='||tid;
END LOOP;


END LOOP;

DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid;
RETURN NULL;

END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE STRICT;


-- Function: reaction_match_substruct(bytea, int4, bool, bool, bool)

-- DROP FUNCTION reaction_match_substruct(bytea, int4, bool, bool, bool);

CREATE OR REPLACE FUNCTION reaction_match_substruct(bytea, int4, bool, bool, bool)
  RETURNS bool AS
$BODY$
DECLARE query_row record;
DECLARE result_row record;
DECLARE update_row record;
DECLARE include text;
DECLARE max_cols text;
DECLARE csum int4;
DECLARE i int4;
DECLARE tid int4;
DECLARE total int4;
--DECLARE must_normalize bool;
BEGIN

tid:=nextval('r_matchmatrix_tid_seq');
--must_normalize:=false;
max_cols:='0';

FOR query_row IN SELECT type,fragment,n_e,n_p,abs_pos FROM decode_rxnfile($1)
LOOP

csum:=0;
PERFORM fastmatch_set_substructure_query_a(query_row.fragment,$3,$4,$5);

FOR result_row IN SELECT molecule,abs_fragment_pos FROM r_fragments,r_rf_crossref WHERE r_rf_crossref.fragment_id = r_fragments.iiid and rxn_id=$2 AND fragment_type=query_row.type AND query_row.n_e<=n_e AND query_row.n_p<=r_rf_crossref.n_p
LOOP
IF (fastmatch_a(result_row.molecule)=true) THEN
	IF (SELECT true FROM r_matchmatrix WHERE target_pos=result_row.abs_fragment_pos) THEN
		EXECUTE 'UPDATE r_matchmatrix set query_'||query_row.abs_pos||'=1 where target_pos='||result_row.abs_fragment_pos;
	ELSE
		EXECUTE 'INSERT INTO r_matchmatrix (query_'||query_row.abs_pos||',target_pos,tid) VALUES (1,'||result_row.abs_fragment_pos||','||tid||')';
	END IF;
	csum:=csum+1;
END IF;
END LOOP;

--RAISE NOTICE '%',csum;
IF (csum=0) THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid; return false; 
--ELSIF (csum>1) THEN must_normalize:=true;
END IF;

END LOOP;

--IF (must_normalize=false) THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid; return true; END IF;

total:=query_row.n_e+query_row.n_p;

WHILE (true)
LOOP

--retval:=false;
include:='';

FOR i IN 1..total LOOP

max_cols:=max_cols||'+query_'||i;

FOR result_row IN EXECUTE 'SELECT sum(query_' || i ||') as col_sum,count(target_pos) as row_count FROM r_matchmatrix where tid='||tid
LOOP
	IF (result_row.row_count IS NULL OR result_row.row_count<total OR result_row.col_sum=0) THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid; return false;
	ELSIF (result_row.col_sum>1) THEN include:=include || 'query_' || i || ',';
	END IF;
END LOOP;

--RAISE NOTICE '%',include;

--IF (retval=false) THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid; 
--return false; END IF;

END LOOP;

IF (include='') THEN DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid;
return true; END IF;

--RAISE NOTICE '%',max_cols;

FOR result_row IN EXECUTE 'SELECT target_pos FROM r_matchmatrix where '|| max_cols ||'>1 and tid='||tid
LOOP
		EXECUTE 'UPDATE r_matchmatrix set '|| split_part(include,',',1) ||'=0 WHERE target_pos='||result_row.target_pos||' and tid='||tid;
END LOOP;


END LOOP;

DELETE FROM r_matchmatrix WHERE r_matchmatrix.tid=tid;
RETURN NULL;

END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE STRICT;
ALTER FUNCTION reaction_match_substruct(bytea, int4, bool, bool, bool) OWNER TO chembank;




