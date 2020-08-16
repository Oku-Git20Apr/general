#!/bin/csh -f

foreach f (conf/*.dat)

uniq ${f} ${f}_uniq


end
