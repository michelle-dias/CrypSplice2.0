import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks_with_hic = """
[hic matrix]
file = Li_et_al_2015.h5
title = depth = 200000; transform = log1p; min_value = 5
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = false

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 250000; orientation = inverted; customized colormap; min_value = 5; max_value = 70
min_value = 5
max_value = 70
depth = 250000
colormap = ['white', (1, 0.88, .66), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
file_type = hic_matrix
show_masked_bins = false
orientation = inverted

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 300000; transform = log1p; colormap Blues (TADs: overlay_previous = share-y; line_width = 1.5)
colormap = Blues
min_value = 10
max_value = 150
depth = 300000
transform = log1p
file_type = hic_matrix

[tads]
file = tad_classification.bed
#title = TADs color = none; border_color = black
file_type = domains
border_color = black
color = none
height = 5
line_width = 1.5
overlay_previous = share-y
show_data_range = false

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 250000; transform = log1p; colormap = bone_r (links: overlay_previous = share-y; links_type = triangles; color = darkred; line_style = dashed, bigwig: color = red)
colormap = bone_r
min_value = 15
max_value = 200
depth = 250000
transform = log1p
file_type = hic_matrix
show_masked_bins = false

[test arcs]
file = links2.links
title =
links_type = triangles
line_style = dashed
overlay_previous = share-y
line_width = 0.8
color = darkred
show_data_range = false


[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
height = 4
title =
overlay_previous = yes
min_value = 0
max_value = 50
show_data_range = false

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 200000; show_masked_bins = true; colormap = ['blue', 'yellow', 'red']; max_value = 150
depth = 200000
colormap = ['blue', 'yellow', 'red']
max_value = 150
file_type = hic_matrix
show_masked_bins = true

[spacer]
height = 0.1

[x-axis]

"""

with open(os.path.join(ROOT, "browser_tracks_hic.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)

with open(os.path.join(ROOT, "browser_tracks_cool.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic.replace(".h5", ".cool"))

browser_tracks_with_hic = """
[hic matrix]
file = Li_et_al_2015.cool
title = depth = 200000; transform = log1p; min_value = 5; height = 5
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = false
height = 5

[hic matrix 2]
file = Li_et_al_2015.h5
title = same but orientation=inverted; no height
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = false
orientation = inverted

[spacer]
height = 0.5

[hic matrix 3]
file = Li_et_al_2015.h5
title = same rasterize = false
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
rasterize = false
show_masked_bins = false

[x-axis]

"""

with open(os.path.join(ROOT, "browser_tracks_hic_rasterize_height.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)

browser_tracks_with_hic = """
[hic matrix]
file = small_test2.cool
title = cool with few interactions show_masked_bins = false (default)
depth = 200000
file_type = hic_matrix
height = 5

[spacer]

[hic matrix]
file = small_test2.cool
title = cool with few interactions show_masked_bins = true
depth = 200000
file_type = hic_matrix
height = 5
show_masked_bins = true

[x-axis]
"""
with open(os.path.join(ROOT, "browser_tracks_hic_small_test.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)
with open(os.path.join(ROOT, "browser_tracks_hic_small_test_invalid_colormap1.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic.replace('[x-axis]', 'colormap = [\n[x-axis]'))
with open(os.path.join(ROOT, "browser_tracks_hic_small_test_invalid_colormap2.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic.replace('[x-axis]', "colormap = ['a']\n[x-axis]"))
with open(os.path.join(ROOT, "browser_tracks_hic_small_test_invalid_colormap3.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic.replace('[x-axis]', "colormap = fakeColormap\n[x-axis]"))
with open(os.path.join(ROOT, "browser_tracks_hic_small_test_invalid_colormap4.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic.replace('[x-axis]', "colormap = []\n[x-axis]"))
with open(os.path.join(ROOT, "browser_tracks_hic_small_test_invalid_colormap5.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic.replace('[x-axis]', "colormap = ['white', (1, 0.88, 2./3)]\n[x-axis]"))
with open(os.path.join(ROOT, "browser_tracks_hic_small_test_invalid_colormap6.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic.replace('[x-axis]', "colormap = 'white'\n[x-axis]"))


browser_tracks_with_hic = """
[hic matrix log]
file = Li_et_al_2015.h5
title = depth = 200000; transform = log; show_masked_bins = true
depth = 200000
transform = log
min_value = -4
max_value = 4
file_type = hic_matrix
show_masked_bins = true

[spacer]

[hic matrix -log]
file = Li_et_al_2015.h5
title = depth = 200000; transform = -log; show_masked_bins = true
depth = 200000
transform = -log
min_value = -4
max_value = 4
file_type = hic_matrix
show_masked_bins = true

[x-axis]
"""

with open(os.path.join(ROOT, "browser_tracks_hic_log-log.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)

browser_tracks_with_mcool = """
[x-axis]

[mcool1]
file = matrix.mcool::/0
depth = 1000000
file_type = hic_matrix

[mcool2]
file = matrix.mcool::/1
depth = 1000000
file_type = hic_matrix

[mcool3]
file = matrix.mcool::/2
depth = 1000000
file_type = hic_matrix

[mcool4]
file = matrix.mcool::/4
depth = 1000000
file_type = hic_matrix
"""

with open(os.path.join(ROOT, "mcool.ini"), 'w') as fh:
    fh.write(browser_tracks_with_mcool)

browser_tracks_with_hic_small_2 = """
[hic matrix]
file = small_test2.cool
title = cool with few interactions depth=10kb
file_type = hic_matrix
depth = 10000


[x-axis]
"""

with open(os.path.join(ROOT, "browser_tracks_hic_small_test_2.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic_small_2)

browser_tracks_with_one_interaction = """
[hic matrix]
file = one_interaction_4chr.cool
title = cool with one interaction show_masked_bins = false
depth = 200000
file_type = hic_matrix
height = 5
show_masked_bins = false

[spacer]

[hic matrix]
file = one_interaction_4chr.cool
title = cool with one interaction show_masked_bins = true
depth = 200000
file_type = hic_matrix
height = 5
show_masked_bins = true

[spacer]

[hic matrix]
file = one_interaction_4chr.cool
title = cool with one interaction show_masked_bins = true transform = log
depth = 200000
file_type = hic_matrix
height = 5
show_masked_bins = true
transform = log

[newtrack]
file_type = x_axis
"""

with open(os.path.join(ROOT, "browser_tracks_hic_one_interaction_cool.ini"), 'w') as fh:
    fh.write(browser_tracks_with_one_interaction)

# I keep cool in the title but I use the h5 matrix
with open(os.path.join(ROOT, "browser_tracks_hic_one_interaction_h5.ini"), 'w') as fh:
    fh.write(browser_tracks_with_one_interaction.replace('one_interaction_4chr.cool', 'one_interaction_4chr.h5'))
# I keep cool in the title but I use the h5 matrix
with open(os.path.join(ROOT, "browser_tracks_hic_one_interaction_diag_h5.ini"), 'w') as fh:
    fh.write(browser_tracks_with_one_interaction.replace('one_interaction_4chr.cool', 'one_interaction_diag_4chr.h5'))
# >>> from HiCMatrix import hicmatrix
# >>> hic_ma = HiCMatrix.hiCMatrix(os.path.join(ROOT, 'small_test2.cool'))
# >>> instances, features = hic_ma.matrix.nonzero()
# >>> instances
# array([   0,    0,    1,    1,    2,    3,    3,    3, 9166], dtype=int32)
# >>> features
# array([   0,    3,    2,    3,    1,    0,    1,    3, 9166], dtype=int32)
# >>> hic_ma.matrix.data
# array([1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=int32)
# >>> import numpy as np
# >>> hic_ma.matrix.data = np.array([10, 5, 2, 1, 2, 5, 1, 100, -20])
# >>> hic_ma.save('pygenometracks/tests/test_data/small_test3.cool')

browser_tracks_with_hic_small_3 = """
[hic matrix]
file = small_test3.cool
title = cool with few interactions transform = no
depth = 200000
file_type = hic_matrix
min_value = 0
max_value = 50

[hic matrix]
file = small_test3.cool
title = cool with few interactions transform = log
transform = log
depth = 200000
file_type = hic_matrix
min_value = 0
max_value = 5

[hic matrix]
file = small_test3.cool
title = cool with few interactions transform = -log
transform = -log
depth = 200000
file_type = hic_matrix
min_value = 0
max_value = -5

[hic matrix]
file = small_test3.cool
title = cool with few interactions transform = log1p
transform = log1p
depth = 200000
file_type = hic_matrix
min_value = 0.1
max_value = 50

[x-axis]
"""

with open(os.path.join(ROOT, "browser_tracks_hic_small_test_3.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic_small_3)

browser_tracks_with_hic_small_3b = """
[hic matrix]
file = small_test3.cool
title = cool with few interactions
depth = 200000
file_type = hic_matrix
transform = log1p

[x-axis]
"""

with open(os.path.join(ROOT, "browser_tracks_hic_small_test_3_invalid.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic_small_3b)

with open(os.path.join(ROOT, "browser_tracks_hic_small_test_3_invalid2.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic_small_3b.replace('log1p', 'log'))

with open(os.path.join(ROOT, "browser_tracks_hic_small_test_3_invalid3.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic_small_3b.replace('log1p', '-log'))

browser_tracks_with_hic_force_scale = """
[hic matrix]
file = Li_et_al_2015.h5
title = depth = 500000; colormap = PuRd; min_value = 5; max_value = 70
min_value = 5
max_value = 70
depth = 500000
colormap = PuRd
file_type = hic_matrix
show_masked_bins = false

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 1000000; colormap = PuRd; min_value = 0; max_value = 80
min_value = 0
max_value = 80
depth = 1000000
colormap = PuRd
file_type = hic_matrix
show_masked_bins = false

[x-axis]
"""

with open(os.path.join(ROOT, "browser_tracks_hic_force_scale.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic_force_scale)


tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_tracks_with_hic():
    if mpl.__version__ == "3.1.1":
        my_tolerance = 14
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "browser_tracks_hic.ini")
    region = "X:2500000-3500000"
    expected_file = os.path.join(ROOT, 'master_plot_hic.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_dec():
    if mpl.__version__ == "3.1.1":
        my_tolerance = 14
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "browser_tracks_hic.ini")
    region = "X:2500000-3500000"
    expected_file = os.path.join(ROOT, 'master_plot_hic_dec.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           "--decreasingXAxis "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_cool_region():
    if mpl.__version__ == "3.1.1":
        my_tolerance = 14
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "browser_tracks_cool.ini")
    region = "X:2500000-3500000"
    expected_file = os.path.join(ROOT, 'master_plot_hic.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_hic_logmlog():
    if mpl.__version__ == "3.1.1":
        my_tolerance = 21
    else:
        my_tolerance = tolerance
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "browser_tracks_hic_log-log.ini")
    region = "X:2500000-3500000"
    expected_file = os.path.join(ROOT, 'master_plot_hic_log-log.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_rasterize_height_2chr():
    if mpl.__version__ == "3.1.1" or sys.version_info.minor == 6:
        my_tolerance = 18
    else:
        my_tolerance = tolerance
    extension = '.pdf'
    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "browser_tracks_hic_rasterize_height.ini")
    bed_file = os.path.join(ROOT, 'regions_XY.bed')
    args = f"--tracks {ini_file} --BED {bed_file} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 10 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    for region in ['X:2500000-2600000', 'Y:0-1000000']:
        region_str = region.replace(':', '-')
        output_file = outfile.name[:-4] + '_' + region_str + extension
        expected_file = os.path.join(ROOT, 'master_plot_hic_rasterize_height_'
                                     + region_str + extension)
        # matplotlib compare on pdf will create a png next to it.
        # To avoid issues related to write in test_data folder
        # We copy the expected file into a temporary place
        new_expected_file = NamedTemporaryFile(suffix='.pdf',
                                               prefix='pyGenomeTracks_test_',
                                               delete=False)
        os.system(f'cp {expected_file} {new_expected_file.name}')
        expected_file = new_expected_file.name
        res = compare_images(expected_file,
                             output_file, my_tolerance)
        assert res is None, res

        os.remove(output_file)


def test_plot_tracks_with_hic_rasterize_height_2chr_individual():
    extension = '.pdf'
    ini_file = os.path.join(ROOT, "browser_tracks_hic_rasterize_height.ini")
    for region in ['X:2500000-2600000', 'Y:0-1000000']:
        outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                     delete=False)
        expected_file = os.path.join(ROOT, 'master_plot_hic_rasterize_height_'
                                     + region.replace(':', '-') + extension)
        # matplotlib compare on pdf will create a png next to it.
        # To avoid issues related to write in test_data folder
        # We copy the expected file into a temporary place
        new_expected_file = NamedTemporaryFile(suffix='.pdf',
                                               prefix='pyGenomeTracks_test_',
                                               delete=False)
        os.system(f'cp {expected_file} {new_expected_file.name}')
        expected_file = new_expected_file.name
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.23 --width 38 --dpi 10 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
    assert res is None, res


def test_plot_tracks_with_hic_small_test():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "browser_tracks_hic_small_test.ini")
    bed_file = os.path.join(ROOT, 'regions_chr1XY.bed')
    args = f"--tracks {ini_file} --BED {bed_file} "\
           "--trackLabelFraction 0.23 --width 38 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    for region in ['chr1:0-500000', 'chrX:2500000-2600000', 'chrY:0-1000000']:
        region_str = region.replace(':', '-')
        output_file = outfile.name[:-4] + '_' + region_str + extension
        expected_file = os.path.join(ROOT, 'master_plot_hic_small_test_'
                                     + region_str + extension)
        res = compare_images(expected_file,
                             output_file, tolerance)
        assert res is None, res

        os.remove(output_file)


# The tests with individual chromosome does not give the same result:
# For the empty the problem is the colorbar which is different
# when you do not load all data and the transformation of nan to 0
# when the matrix is not empty

def test_plot_tracks_with_hic_small_test_individual():
    extension = '.png'
    ini_file = os.path.join(ROOT, "browser_tracks_hic_small_test.ini")
    for region in ['chr1:0-500000']:  # , 'chrX:2500000-2600000', 'chrY:-0-1000000']:

        outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                     delete=False)
        expected_file = os.path.join(ROOT, 'master_plot_hic_small_test_'
                                     + region.replace(':', '-') + extension)

        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.23 --width 38 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)


def test_plot_tracks_with_hic_small_other_chr_name():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    region = '1:0-200000'
    expected_file = os.path.join(ROOT, 'master_plot_hic_small_test.png')

    for suf in [''] + ['_invalid_colormap' + s for s in ['1', '2', '3', '4', '5', '6']]:
        ini_file = os.path.join(ROOT, f'browser_tracks_hic_small_test{suf}.ini')

        args = f"--tracks {ini_file} --region {region} "\
            "--trackLabelFraction 0.23 --width 38 "\
            f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)
        if 'invalid' in ini_file:
            os.remove(ini_file)


def test_plot_tracks_with_hic_small_above_chr_length():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_small_test.ini')
    region = 'chrM:0-20000'
    expected_file = os.path.join(ROOT, 'master_plot_hic_small_test_chrM.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_small_above_chr_length_other_chrName():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_small_test.ini')
    region = 'Y:90000000-100000000'
    expected_file = os.path.join(ROOT, 'master_hic_small_test_above_chrY.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_hic_depth_smaller_binsize():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_small_test_2.ini')
    region = 'chr1:0-100000'
    expected_file = os.path.join(ROOT, 'master_hic_small_test_2.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_hic_plotting_region_smaller_binsize():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_small_test.ini')
    region = 'chr1:0-5000'
    expected_file = os.path.join(ROOT, 'master_hic_small_test_small_region.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_mcool():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:2500000-3500000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT,
                                  'mcool.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT,
                                      'master_mcool.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_one_interaction_cool():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "browser_tracks_hic_one_interaction_cool.ini")
    bed_file = os.path.join(ROOT, 'regions_chr1XY.bed')
    args = f"--tracks {ini_file} --BED {bed_file} "\
           "--trackLabelFraction 0.23 --width 38 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    for region in ['chr1:0-500000', 'chrX:2500000-2600000', 'chrY:0-1000000']:
        region_str = region.replace(':', '-')
        output_file = outfile.name[:-4] + '_' + region_str + extension
        expected_file = os.path.join(ROOT, 'master_plot_hic_one_interaction_withBED_'
                                     + region_str + extension)
        res = compare_images(expected_file,
                             output_file, tolerance)
        assert res is None, res

        os.remove(output_file)


def test_plot_tracks_with_hic_one_interaction_h5():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "browser_tracks_hic_one_interaction_h5.ini")
    bed_file = os.path.join(ROOT, 'regions_chr1XY.bed')
    args = f"--tracks {ini_file} --BED {bed_file} "\
           "--trackLabelFraction 0.23 --width 38 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    for region in ['chr1:0-500000', 'chrX:2500000-2600000', 'chrY:0-1000000']:
        region_str = region.replace(':', '-')
        output_file = outfile.name[:-4] + '_' + region_str + extension
        expected_file = os.path.join(ROOT, 'master_plot_hic_one_interaction_withBED_'
                                     + region_str + extension)
        res = compare_images(expected_file,
                             output_file, tolerance)
        assert res is None, res

        os.remove(output_file)


def test_plot_tracks_with_hic_one_interaction_h5_individual():
    extension = '.png'
    ini_file = os.path.join(ROOT, "browser_tracks_hic_one_interaction_h5.ini")
    for region in ['chr1:0-500000', 'chrX:2500000-2600000', 'chrY:0-1000000']:

        outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                     delete=False)
        expected_file = os.path.join(ROOT, 'master_plot_hic_one_interaction_withBED_'
                                     + region.replace(':', '-') + extension)

        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.23 --width 38 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)


# The tests with individual chromosome with cool does not give the same result:
# For the empty the problem is the colorbar which is different
# when you do not load all data and the transformation of nan to 0
# when the matrix is not empty
# In addition here because there was no data on diagonal it has been modified.
# When you load only one empty chromosome it is not empty anymore
def test_plot_tracks_with_hic_one_interaction_individual():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_one_interaction_cool.ini')
    region = 'chr1:0-500000'
    expected_file = os.path.join(ROOT, f"master_plot_hic_one_interaction_withBED_{region.replace(':', '-')}.png")
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_one_interaction_cool.ini')
    region = 'chrY:0-1000000'
    expected_file = os.path.join(ROOT, 'master_plot_hic_one_interaction_withRegion_chrY-0-1000000.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_one_interaction_diag():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_one_interaction_diag_h5.ini')
    regions = ['chrY:0-1000000', 'chrY:130000000-131000000']
    for region in regions:
        expected_file = os.path.join(ROOT, 'master_plot_hic_one_interaction_diag_'
                                     + region.replace(':', '-') + '.png')

        args = f"--tracks {ini_file} --region {region} "\
            "--trackLabelFraction 0.23 --width 38 "\
            f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)


def test_plot_tracks_with_hic_small3():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_small_test_3.ini')
    region = 'chr1:0-200000'
    expected_file = os.path.join(ROOT, 'master_hic_small_test_3.png')

    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_small3_invalid():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    region = "chrM:0-20000"
    for suf in ['', '2', '3']:
        ini_file = os.path.join(ROOT, f"browser_tracks_hic_small_test_3_invalid{suf}.ini")
        args = f"--tracks {ini_file} --region {region} "\
               f"--outFileName {outfile.name}".split()
        try:
            pygenometracks.plotTracks.main(args)
        except Exception as e:
            assert 'transformation can not be applied to' in str(e)
        else:
            raise Exception("The plot_tracks_with_hic_small3_invalid should fail.")
        os.remove(ini_file)


def test_plot_tracks_with_hic_force_scale():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'browser_tracks_hic_force_scale.ini')
    region = 'X:2500000-3500000'
    expected_file = os.path.join(ROOT, 'master_plot_hic_force_scale.png')

    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           "--fontSize 16 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
