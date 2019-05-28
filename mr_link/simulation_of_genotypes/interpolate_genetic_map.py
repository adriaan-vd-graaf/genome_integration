import gzip

"""
This script will take the variants from 1000g and interpolate the centimorgan position, based on the genetic map
of HAPMAP3 overlapping positions.
"""


ordered_1000g_positions = []
ordered_1000g_lines = []
with open("genetic_variants_in_region.txt", "r") as f:
    for line in f:
        ordered_1000g_positions.append(int(line.split()[1]))
        ordered_1000g_lines.append(line[:-1])


ordered_gmap_positions = []
ordered_gpos = []
ordered_gmap_cm = []
ordered_gmap_line = []
with open("genetic_map_GRCh37_chr2.txt", "r") as f:
    f.readline()
    for line in f:
        split = line.split()
        ordered_gmap_positions.append(int(split[1]))
        ordered_gmap_cm.append(float(split[2]))
        ordered_gpos.append(float(split[3]))
        ordered_gmap_line.append(line)



interpolated_genetic_map = open("interpolated_genetic_map.100Mb.to.105Mb.map", "w")
interpolated_genetic_map.write("position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n")

thousand_g_indice = 0
gmap_indice = 0

map_position = 0.0
prev_pos = None
map_zero = None

while thousand_g_indice < len(ordered_1000g_positions):
    ref_pos = ordered_1000g_positions[thousand_g_indice]
    gmap_pos = ordered_gmap_positions[gmap_indice]

    if gmap_pos < (len(ordered_gmap_positions)-1) and \
        ref_pos < gmap_pos and \
        ref_pos < ordered_gmap_positions[gmap_indice + 1]:
        gmap_pos +=1

    elif ref_pos == gmap_pos:
        if map_zero is None:
            map_zero = ordered_gpos[gmap_indice]

        interpolated_genetic_map.write("{} {:.6f} {:.6f}\n".format(ref_pos, ordered_gmap_cm[gmap_indice], ordered_gpos[gmap_indice] - map_zero))
        gmap_indice += 1
        thousand_g_indice += 1

        printed_cm = ordered_gmap_cm[gmap_indice]

    elif ref_pos < gmap_pos:
        if gmap_indice == 0:
            interpolated_genetic_map.write("{} {} {}\n".format(ref_pos, ordered_gmap_cm[gmap_indice], map_position))
            thousand_g_indice += 1
        else:

            if map_zero is None:
                map_zero = ordered_gpos[gmap_indice]

            prev_cm = ordered_gmap_cm[gmap_indice - 1]
            curr_cm = ordered_gmap_cm[gmap_indice]
            gmap_pos_prev = ordered_gmap_positions[gmap_indice-1]
            gmap_pos_curr = ordered_gmap_positions[gmap_indice]
            frac = (ref_pos - gmap_pos)  / (gmap_pos_curr - gmap_pos_prev)
            interpolated_cm = curr_cm + frac * (curr_cm - prev_cm )
            printed_cm = interpolated_cm
            interpolated_genetic_map.write("{} {:.6f} {:.6f}\n".format(ref_pos, interpolated_cm, ordered_gpos[gmap_indice] - map_zero))
            thousand_g_indice += 1
    elif ref_pos > gmap_pos:
        gmap_indice += 1


