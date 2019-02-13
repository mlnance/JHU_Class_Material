load 3ay4.pdb

space cmyk
bg_color white
set cartoon_side_chain_helper, on
set cartoon_transparency, 0.4
remove resn HOH

# for FcgRIIIa
hide everything, chain C

# for IgG1
hide everything, chain B
select chain A
util.cbas sele
select chain A and organic
util.cbak sele
show sticks, resi 297 and chain A
hide everything, resi 229-237 and chain A
deselect

# pi interactions
show sticks, resi 241+243 and chain A
select !(name n+c+ca+o) and (resi 241+243 and chain A)
util.cbac sele

# view
set_view (\
     0.591056406,    0.246667117,   -0.767980754,\
     0.079328969,   -0.965257525,   -0.248977602,\
    -0.802718103,    0.086240076,   -0.590088427,\
     0.000023870,    0.000215431,  -78.314102173,\
    12.517524719,   21.217731476,  128.566696167,\
  -31372.626953125, 31530.068359375,  -20.000000000 )
dss

# image
set use_shaders, on
set opaque_background, off
# not using ray
draw 4000, 3000
png MLN_IgG1_pi_stacking.png
