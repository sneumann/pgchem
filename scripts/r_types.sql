-- Type: public.reactionmap

-- DROP TYPE public.reactionmap;

CREATE TYPE public.reactionmap AS
   ("type" int4,
    abs_pos int4,
    rel_pos int4,
    fragment bytea,
    n_e int4,
    n_p int4,
    n_total int4,
    ms_fingerprint_long text);
ALTER TYPE public.reactionmap OWNER TO chembank;
