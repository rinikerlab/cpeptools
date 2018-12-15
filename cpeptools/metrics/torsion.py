
#FIXME
indices = get_1_4_pairs(look_up[j], "[O:1]=[C:2]@;-[NX3:3]-[CX4H3:4]") + get_1_4_pairs(look_up[j], "[O:1]=[C:2]@;-[NX3:3]-[H:4]")
out = md.compute_dihedrals(traj, indices)

raise NotImplementedError()
#TODO per-residue Ramanchandran plots
