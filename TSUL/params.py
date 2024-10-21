resolution = 10

lake_alt = 1767
dam_alt = lake_alt + 10
flood_seed = (2.67026e6, 1.171490e6)

bounds = dict(
    xmin = 2669850.,
    xmax = 2670900.,
    ymin = 1170300.,
    ymax = 1171960.,
)


out_format = "binary"

nx = 20
ny = int(1.62*nx)

amr_ratios = dict(
    x=(2, 2),
    y=(2, 2),
    t=(2, 2),
)

