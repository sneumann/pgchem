#!/bin/sh
PATH=$PATH:/opt/checkmol/bin
export PATH
/opt/postgres/bin/pg_ctl -o -i -D /var/pgsql/data/ start
