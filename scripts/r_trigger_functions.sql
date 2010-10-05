-- Function: t_decompose_reaction()

-- DROP FUNCTION t_decompose_reaction();

CREATE OR REPLACE FUNCTION t_decompose_reaction()
  RETURNS "trigger" AS
$BODY$
DECLARE molstatistics_combined_hash_lhs text;
DECLARE molstatistics_combined_hash_rhs text;
DECLARE molstatistics_values text;
DECLARE current_row record;
DECLARE product_fragment_pos integer;
DECLARE educt_fragment_pos integer;
DECLARE current_iiid integer;
BEGIN

molstatistics_combined_hash_lhs:='';
molstatistics_combined_hash_rhs:='';
product_fragment_pos=0;
educt_fragment_pos=0;

IF TG_OP='INSERT' THEN

FOR current_row IN select type ,abs_pos,rel_pos, fragment, n_e,n_p,n_total,md5(ms_fingerprint_long) as hash from decode_rxnfile(NEW.rxn) ORDER BY type,hash
LOOP

current_iiid=NULL;

molstatistics_values:=ms_fingerprint_short_a(current_row.fragment);

SELECT iiid FROM r_fragments WHERE molstatistics_values = ms_fingerprint_short_a(molecule) AND match_exact_a(current_row.fragment,molecule,false,false,false)=true INTO current_iiid;

IF (current_iiid IS NULL) THEN
INSERT INTO r_fragments (molecule) VALUES (current_row.fragment);
current_iiid:=currval('r_fragments_iiid_seq');
END IF;

IF current_row.type=1 THEN
	molstatistics_combined_hash_lhs:=molstatistics_combined_hash_lhs || '+' || current_row.hash;
	educt_fragment_pos:=educt_fragment_pos+1;
--INSERT INTO r_fragments (rxn_id, fragment_type,rel_fragment_pos,molecule,ordered_fragment_pos,n_e,n_p,abs_fragment_pos) VALUES (NEW.rxn_id, current_row.type, current_row.rel_pos, current_row.fragment, educt_fragment_pos,current_row.n_e,current_row.n_p,current_row.abs_pos);
INSERT INTO r_rf_crossref (rxn_id, fragment_type,rel_fragment_pos,fragment_id,ordered_fragment_pos,n_e,n_p,abs_fragment_pos) VALUES (NEW.rxn_id, current_row.type, current_row.rel_pos, current_iiid, educt_fragment_pos,current_row.n_e,current_row.n_p,current_row.abs_pos);
ELSIF current_row.type=2 THEN
	molstatistics_combined_hash_rhs:=molstatistics_combined_hash_rhs || '+' || current_row.hash;
	product_fragment_pos:=product_fragment_pos+1;
--INSERT INTO r_fragments (rxn_id, fragment_type,rel_fragment_pos,molecule,ordered_fragment_pos,n_e,n_p,abs_fragment_pos) VALUES (NEW.rxn_id, current_row.type, current_row.rel_pos, current_row.fragment, product_fragment_pos,current_row.n_e,current_row.n_p,current_row.abs_pos);
--INSERT INTO r_fragments (molecule) VALUES (current_row.fragment);
INSERT INTO r_rf_crossref (rxn_id, fragment_type,rel_fragment_pos,fragment_id,ordered_fragment_pos,n_e,n_p,abs_fragment_pos) VALUES (NEW.rxn_id, current_row.type, current_row.rel_pos, current_iiid, product_fragment_pos,current_row.n_e,current_row.n_p,current_row.abs_pos);
END IF;

END LOOP;

molstatistics_combined_hash_lhs:=btrim(molstatistics_combined_hash_lhs,'+');
molstatistics_combined_hash_rhs:=btrim(molstatistics_combined_hash_rhs,'+');

molstatistics_combined_hash_lhs:=molstatistics_combined_hash_lhs || '->' || molstatistics_combined_hash_rhs;

INSERT INTO r_hash (rxn_id, rxn_hash) VALUES (NEW.rxn_id, molstatistics_combined_hash_lhs);

ELSIF TG_OP='UPDATE' THEN

	DELETE FROM r_rf_crossref WHERE rxn_id=NEW.rxn_id;

FOR current_row IN select type ,abs_pos,rel_pos, fragment, n_e,n_p,n_total,md5(ms_fingerprint_long) as hash from decode_rxnfile(NEW.rxn) ORDER BY type,hash
LOOP

current_iiid=NULL;

molstatistics_values:=ms_fingerprint_short_a(current_row.fragment);

SELECT iiid FROM r_fragments WHERE molstatistics_values = ms_fingerprint_short_a(molecule) AND match_exact_a(current_row.fragment,molecule,false,false,false)=true INTO current_iiid;

IF (current_iiid IS NULL) THEN
INSERT INTO r_fragments (molecule) VALUES (current_row.fragment);
current_iiid:=currval('r_fragments_iiid_seq');
END IF;

IF current_row.type=1 THEN
	molstatistics_combined_hash_lhs:=molstatistics_combined_hash_lhs || '+' || current_row.hash;
	educt_fragment_pos:=educt_fragment_pos+1;
--INSERT INTO r_fragments (rxn_id, fragment_type,rel_fragment_pos,molecule,ordered_fragment_pos,n_e,n_p,abs_fragment_pos) VALUES (NEW.rxn_id, current_row.type, current_row.rel_pos, current_row.fragment, educt_fragment_pos, current_row.n_e,current_row.n_p,current_row.abs_pos);
INSERT INTO r_rf_crossref (rxn_id, fragment_type,rel_fragment_pos,fragment_id,ordered_fragment_pos,n_e,n_p,abs_fragment_pos) VALUES (NEW.rxn_id, current_row.type, current_row.rel_pos, current_iiid, educt_fragment_pos,current_row.n_e,current_row.n_p,current_row.abs_pos);

ELSIF current_row.type=2 THEN
	molstatistics_combined_hash_rhs:=molstatistics_combined_hash_rhs || '+' || current_row.hash;
	product_fragment_pos:=product_fragment_pos+1;
--INSERT INTO r_fragments (rxn_id, fragment_type,rel_fragment_pos,molecule,ordered_fragment_pos,n_e,n_p,abs_fragment_pos) VALUES (NEW.rxn_id, current_row.type, current_row.rel_pos, current_row.fragment, product_fragment_pos,current_row.n_e,current_row.n_p,current_row.abs_pos);
INSERT INTO r_rf_crossref (rxn_id, fragment_type,rel_fragment_pos,fragment_id,ordered_fragment_pos,n_e,n_p,abs_fragment_pos) VALUES (NEW.rxn_id, current_row.type, current_row.rel_pos, current_iiid, product_fragment_pos,current_row.n_e,current_row.n_p,current_row.abs_pos);
END IF;

END LOOP;

molstatistics_combined_hash_lhs:=btrim(molstatistics_combined_hash_lhs,'+');
molstatistics_combined_hash_rhs:=btrim(molstatistics_combined_hash_rhs,'+');

molstatistics_combined_hash_lhs:=molstatistics_combined_hash_lhs || '->' || molstatistics_combined_hash_rhs;

 DELETE FROM r_hash WHERE rxn_id=NEW.rxn_id;
  INSERT INTO r_hash (rxn_id, rxn_hash) VALUES (NEW.rxn_id, molstatistics_combined_hash_lhs);

ELSE

  RAISE EXCEPTION 'PGCHEM INDEX MAINTENANCE TRIGGER CALLED OUTSIDE INSERT OR UPDATE!';

END IF;

RETURN NULL;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE STRICT;
ALTER FUNCTION t_decompose_reaction() OWNER TO postgres;

-- Function: public.t_maintain_r_indextables()

-- DROP FUNCTION public.t_maintain_r_indextables();

CREATE OR REPLACE FUNCTION public.t_maintain_r_indextables()
  RETURNS "trigger" AS
$BODY$
DECLARE molstatistics_values text;
DECLARE molstatistics_hash text;
BEGIN

  molstatistics_values:=ms_fingerprint_short_a(NEW.molecule);

IF TG_OP='INSERT' THEN
  
  EXECUTE 'INSERT INTO r_molstatistics VALUES (' || NEW.iiid || ',' || molstatistics_values || ');';
 

ELSIF TG_OP='UPDATE' THEN

  DELETE FROM r_molstatistics WHERE iiid=OLD.iiid;
  EXECUTE 'INSERT INTO r_molstatistics VALUES (' || NEW.iiid || ',' || molstatistics_values || ');';
  
ELSE

  RAISE EXCEPTION 'PGCHEM INDEX MAINTENANCE TRIGGER CALLED OUTSIDE INSERT OR UPDATE!';

END IF;

RETURN NULL;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE STRICT;
ALTER FUNCTION public.t_maintain_r_indextables() OWNER TO postgres;

