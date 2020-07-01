## Validation procedure

Model has a global setup. Results at validation gauges will just be
retrieved from model run.

Model has 1 optimal parameter sets (same for objective 1 and objective
2).

I have done an analysis on the shape-files and FYI these are the
matching with the old/original (724-subbasin routing network) that I
used: 

04185000  >>> in-channel (subid: 344) + in-basin (subid: 327 *0.8)

04195500  >>> in-channel (subid: 159) 

04201500  >>> in-channel (subid: 111 + 110) +  in-basin (subid: 613 *0.1)

04208000  >>> in-channel (subid: 108 + 109)
Here there is an issue with the provided watershed for 04208000 that
id 4417 and 4419 drain into id 7033 but ID 7033 is not given in the
shapefile as the most downstream (basically there are two downstream
sub-basins). I assume the gages is downstream of both 4417 and 4419 in
the GL network (the sum of them). Right? 

04168000  >>> in-basin (subid: 211)  # no in channel routing, NS
result for this gauges can be low 

04167000  >>> in-basin (subid: 196)  # no in channel routing, NS
result for this gauges can be low 

02GE003   >>> in-channel (subid: 324) + in-basin (subid: 312 *0.6)
