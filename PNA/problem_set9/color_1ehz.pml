load 1ehz.pdb
remove resn hoh

select acceptor_helix, resi 1-7+66-76
color red, acceptor_helix
select variable_loop, resi 8-11
color orange, variable_loop
select D_stem-loop, resi 12-23
color yellow, D_stem-loop
select anticodon_stem-loop, resi 24-45
color green, anticodon_stem-loop
select link, resi 46-48
color blue, link
select T_stem-loop, resi 49-65
color purple, T_stem-loop
deselect