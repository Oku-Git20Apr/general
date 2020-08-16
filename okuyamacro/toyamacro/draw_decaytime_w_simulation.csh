#!/bin/tcsh

foreach NL (110 120 130 140 150 160 170 180 190 200)
  foreach mass_min (1.105 1.110 1.113)
    foreach t_width (180 190)
      switch(${mass_min})
        case 1.105 :
          set mass_max = 1.125
          @ NL_a = ${NL}
          breaksw
        case 1.110 :
          set mass_max = 1.120
          @ NL_a = ${NL} - 40
          breaksw
        case 1.113 :
          set mass_max = 1.117
          @ NL_a = ${NL} - 100
          breaksw
        default :
          set mass_max = 1.125
      endsw

echo "NL:"${NL_a} " Mmin:"${mass_min}
echo "min:"${mass_min} " max:"${mass_max}
root -b  << EOF
.L draw_decaytime_w_simulation.C
draw_decaytime_w_simulation(${NL_a},${mass_min},${mass_max},${t_width});
EOF

    end
  end
end
