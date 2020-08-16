#!/bin/csh

foreach min (1.08 1.09 1.10)
  foreach max (1.13 1.14 1.16 1.20 1.25)

root -b  << EOF
.L decaytime_SNdep.C
decaytime_SNdep(${min},${max});
EOF

  end
end
