# Maaike made a better version. This is WIDLY outdated.


from hcipy import *
import hcipy
import numpy as np
import matplotlib.pyplot as plt


def make_keck_aperture(normalized=False, with_spiders=False, with_segment_gaps=False, gap_padding=0, segment_transmissions=1, return_header=False, return_segments=False):
    pupil_diameter = 10.9 #m actual circumscribed diameter
    actual_segment_flat_diameter = np.sqrt(3)/2 * 1.8 #m actual segment flat-to-flat diameter
    # iris_ao_segment = np.sqrt(3)/2 * .7 mm (~.606 mm)
    actual_segment_gap = 0.003 #m actual gap size between segments
    # (3.5 - (3 D + 4 S)/6 = iris_ao segment gap (~7.4e-17)
    spider_width = 1*0.02450 #m actual strut size
    gap_padding = 10.
    segment_gap = actual_segment_gap * gap_padding #padding out the segmentation gaps so they are visible and not sub-pixel
    segment_transmissions = 1.

    segment_flat_diameter = actual_segment_flat_diameter - (segment_gap - actual_segment_gap)
    segment_circum_diameter = 2 / np.sqrt(3) * segment_flat_diameter #segment circumscribed diameter

    num_rings = 3 #number of full rings of hexagons around central segment

    segment_positions = make_hexagonal_grid(actual_segment_flat_diameter + actual_segment_gap, num_rings)
    segment_positions = segment_positions.subset(lambda grid: ~(circular_aperture(segment_circum_diameter)(grid) > 0))

    segment = hexagonal_aperture(segment_circum_diameter, np.pi / 2)

    spider1 = make_spider_infinite([0, 0], 0, spider_width)
    spider2 = make_spider_infinite([0, 0], 60, spider_width)
    spider3 = make_spider_infinite([0, 0], 120, spider_width)
    spider4 = make_spider_infinite([0, 0], 180, spider_width)
    spider5 = make_spider_infinite([0, 0], 240, spider_width)
    spider6 = make_spider_infinite([0, 0], 300, spider_width)

    segmented_aperture = make_segmented_aperture(segment, segment_positions, segment_transmissions, return_segments=return_segments)

    def func(grid):
        res = segmented_aperture(grid) * spider1(grid) * spider2(grid) * spider3(grid)* spider4(grid) * spider3(grid)* spider5(grid) * spider6(grid) # * coro(grid)
        return Field(res, grid)

    return func


wavelength = 638e-9 # who knows if this is right?
keck_aperture = evaluate_supersampled(make_keck_aperture(), pupil_grid, 6)
wf = hcipy.Wavefront(keck_aperture, wavelength)
pupil_grid = make_pupil_grid(256, 1.5) # this might be totally wrong whoops
focal_grid = make_focal_grid(8, 12)  # same tho
prop = FraunhoferPropagator(pupil_grid, focal_grid)
img_ref = prop(wf).intensity


telescope_diameter = 10.9
num_pupil_pixels = 240 #based on PWFS sim
# we should chose this based on image sampling
# make sure our sampling is fine enough to relate well to science plane
pupil_grid_diameter = telescope_diameter
pupil_grid = make_pupil_grid(num_pupil_pixels, pupil_grid_diameter)

keck_aperture  =evaluate_supersampled(make_keck_aperture(telescope_diameter), pupil_grid, 6)

imshow_field(keck_aperture)
plt.xlabel('x position(m)')
plt.ylabel('y position(m)')
plt.colorbar()
plt.show()
# check the PDF -- for no segment gaps + zero intensity
