#!/bin/csh -f

#plz execute this script at macros/toyamacro


  set period = 2017_04
  set base_dir  = `pwd`/../../
##group1##
#foreach run (10226 10227 10228 \
#             10229 10230 10232 10233 \
#             10234 10236 10237 )

##group2##st10238 
#foreach run (10238 10239 10240 10241 \
#             10242 10243 10245  )
##group3##st10246 
#foreach run (10246 10247 10248 10249 \
#             10250 10251 10252 10253 \
#             10254 10255 10256 )

##group4##st10257  
#foreach run (10257 10258 10261  \
#             10262 10268 10269 )

##group5##st10270 
#foreach run (10270 10271 10272 10273  \
#             10274 10275 10276 10277  \
#             10278 10280  )

##group6##st10282 
#foreach run (10282 10283 10284 10286  \
#             10287 10288 10289 10290  \
#             10291 10292 10293 )

##group7##st10294 
#foreach run (10294 10296 10297 10298  \
#            )

##group8##st10310 
#foreach run (10304 10305 10307 10308  \
#             10309 10310 \
#             10311 10312 10313  \
#             10314 10315 10316 10317  \
#             10318 10320  )

##group9##st10321 
foreach run (10321 10322 10323 10324  \
             10325 10326 10342 10345  \
             10346 10347 10348 10350  )
echo ${run}

root -l -b << EOF
.L copy_param.C
copy_param("param/${run}.param");
EOF

end
