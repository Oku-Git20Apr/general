#!/bin/csh

foreach t_width (180 190)
  foreach mass_min (1.105 1.110 1.113)
    switch(${mass_min})
      case 1.105 :
        set mass_max = 1.125
        breaksw
      case 1.110 :
        set mass_max = 1.120
        breaksw
      case 1.113 :
        set mass_max = 1.117
        breaksw
      default :
        set mass_max = 1.125
    endsw

root -b  << EOF
.L draw_decaytime_w_simulation.C
draw_result(${mass_min},${mass_max},${t_width});
EOF

  end
end

