resolution = 5

lake_alt = 1767
dam_alt = lake_alt + 10
flood_seed = (2.67026e6, 1.171490e6)

bounds = dict(
    xmin = 2669700.,
    xmax = 2671000.,
    ymin = 1170100.,
    ymax = 1172100.,
)


out_format = "binary"

nx = 50
ny = int(1.62*nx)

amr_ratios = dict(
    x=(2,),
    y=(2,),
    t=(2,),
)

