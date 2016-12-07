load 1t7p.pdb
hide everything, resn hoh
dist vdw1, /1t7p//P/DG3`23/C5', /1t7p//P/2DA`22/C2'
dist polar1, /1t7p//P/DG3`23/O6, /1t7p//P/2DA`22/N6
dist polar2, /1t7p//P/DG3`23/N2, /1t7p//T/DC`4/O2
dist polar3, /1t7p//P/DG3`23/N1, /1t7p//T/DC`4/N3
dist polar4, /1t7p//P/DG3`23/N2, /1t7p//T/DC`4/N3
dist polar5, /1t7p//P/DG3`23/O6, /1t7p//T/DC`4/N4
dist polar6, /1t7p//P/DG3`23/N2, /1t7p//T/DT`5/O2
dist polar7, /1t7p//P/DG3`23/N2, /1t7p//T/DT`5/N3
dist polar8, /1t7p//P/DG3`23/O4', /1t7p//A/ARG`429/NH2
dist polar9, /1t7p//P/DG3`23/O1G, /1t7p//A/ASP`475/OD1
dist polar10, /1t7p//P/DG3`23/O2A, /1t7p//A/ASP`475/OD1
dist polar11, /1t7p//P/DG3`23/O1G, /1t7p//A/ALA`476/O
dist polar12, /1t7p//P/DG3`23/O2B, /1t7p//A/ALA`476/O
dist polar13, /1t7p//P/DG3`23/O2G, /1t7p//A/GLY`478/N
dist polar14, /1t7p//P/DG3`23/O2B, /1t7p//A/GLY`478/N
dist polar15, /1t7p//P/DG3`23/O2B, /1t7p//A/LEU`479/N
dist vdw2, /1t7p//P/DG3`23/C2', /1t7p//A/GLU`480/CD
dist polar16, /1t7p//P/DG3`23/O1B, /1t7p//A/HIS`506/NE2
dist polar17, /1t7p//P/DG3`23/O3G, /1t7p//A/ARG`518/NH1
dist polar18, /1t7p//P/DG3`23/O2G, /1t7p//A/ARG`518/NH2
dist polar19, /1t7p//P/DG3`23/O3G, /1t7p//A/ARG`518/NH2
dist polar20, /1t7p//P/DG3`23/O3G, /1t7p//A/LYS`522/NZ
dist polar21, /1t7p//P/DG3`23/O1A, /1t7p//A/LYS`522/NZ
dist vdw3, /1t7p//P/DG3`23/C2', /1t7p//A/TYR`526/CZ
dist polar22, /1t7p//P/DG3`23/O1B, /1t7p//A/TYR`526/OH
dist polar23, /1t7p//P/DG3`23/O2A, /1t7p//A/ASP`654/OD1
dist polar24, /1t7p//P/DG3`23/O2B, /1t7p//A/ASP`654/OD2
dist polar25, /1t7p//P/DG3`23/O2A, /1t7p//A/ASP`654/OD2
dist polar26, /1t7p//P/DG3`23/PG, /1t7p//A/MG`4001/MG
dist polar27, /1t7p//P/DG3`23/O1G, /1t7p//A/MG`4001/MG
dist polar28, /1t7p//P/DG3`23/PB, /1t7p//A/MG`4001/MG
dist polar29, /1t7p//P/DG3`23/O2B, /1t7p//A/MG`4001/MG
dist polar30, /1t7p//P/DG3`23/O2A, /1t7p//A/MG`4001/MG
dist polar31, /1t7p//P/DG3`23/PA, /1t7p//A/MG`4002/MG
dist polar32, /1t7p//P/DG3`23/O2A, /1t7p//A/MG`4002/MG
group vdw_contacts, vdw*
group polar_contacts, polar*
show sticks, resn DG3
color magenta, vdw_contacts
color yellow, polar_contacts
hide labels
show spheres, resn MG
set sphere_scale, 0.4
bg_color white
set_view (     0.059504092,    0.576347291,    0.815033674,    -0.203007400,    0.806406200,   -0.555421352,    -0.977367759,   -0.132408097,    0.164987102,     0.000000000,    0.000000000,  -49.289268494,    45.421001434,   24.983999252,    0.598999977,    41.054149628,   61.524391174,  -20.000000000 )