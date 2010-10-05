/**********************************************************************
 * trigger_functions.sql necessary SQL trigger functions
 *
 * Copyright (c) 2004,2005 by Ernst-G. Schmid
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


CREATE OR REPLACE FUNCTION public.t_precompute_properties()
  RETURNS trigger AS
'
BEGIN

IF TG_OP=\'INSERT\' OR TG_OP=\'UPDATE\' THEN

  NEW.molecule:=precompute_properties(NEW.molecule,TRUE);

ELSE

  RAISE EXCEPTION \'PGCHEM PRECOMPUTE MOLECULE PROPERTIES TRIGGER CALLED OUTSIDE INSERT OR UPDATE!\';

END IF;
RETURN NEW;
END;
'
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION public.t_maintain_pgchem_indextables()
  RETURNS trigger AS
'
DECLARE molkeys_values text;
DECLARE fgroup_codes text;
DECLARE fgroup_code text;
DECLARE i integer;
BEGIN

  molkeys_values:=molkeys_short(NEW.molecule);
  fgroup_codes:=public.fgroup_codes(NEW.molecule);
  fgroup_code:=\'\';
  i:=1;

IF TG_OP=\'INSERT\' THEN
  
  INSERT INTO molfingerprints (iiid,fp2bitmap) VALUES (NEW.iiid,fingerprint2(NEW.molecule));
  EXECUTE \'INSERT INTO molkeys VALUES (\' || NEW.iiid || \',\' || molkeys_values || \');\';

	LOOP
		fgroup_code:=split_part(fgroup_codes,\';\',i);
		EXIT WHEN fgroup_code=\'\';
		i:=i+1;
		INSERT INTO molfgroups VALUES (NEW.iiid, fgroup_code);
	END LOOP;

ELSIF TG_OP=\'UPDATE\' THEN

  DELETE FROM molfingerprints WHERE iiid=OLD.iiid;
  INSERT INTO molfingerprints (iiid,fp2bitmap) VALUES (NEW.iiid,fingerprint2(NEW.molecule));
  DELETE FROM molkeys WHERE iiid=OLD.iiid;
  EXECUTE \'INSERT INTO molkeys VALUES (\' || NEW.iiid || \',\' || molkeys_values || \');\';
  DELETE FROM molfgroups WHERE iiid=OLD.iiid;

	LOOP
		fgroup_code:=split_part(fgroup_codes,\';\',i);
		EXIT WHEN fgroup_code=\'\';
		i:=i+1;
		INSERT INTO molfgroups VALUES (NEW.iiid, fgroup_code);
	END LOOP;

ELSE

  RAISE EXCEPTION \'PGCHEM INDEX MAINTENANCE TRIGGER CALLED OUTSIDE INSERT OR UPDATE!\';

END IF;
RETURN NULL;
END;
'
  LANGUAGE 'plpgsql' VOLATILE;
 
