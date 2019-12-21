
load 2OIQ.pdb

select 2OIQ_kin=(2OIQ & POL & chain A)

hide all

show cartoon, 2OIQ_kin
color gray60, 2OIQ_kin

center 2OIQ_kin

# rotate to relevant view point
rotate x, -90
rotate y, 0
rotate z, -15

bg_color white
set ray_shadows, 0
set orthoscopic, on
set ray_opaque_background, off

unset specular
set ray_trace_gain, 0
set ray_trace_mode, 3
set ray_trace_color, black
unset depth_cue

ray
png kinase, dpi=500
