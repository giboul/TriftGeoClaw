Resolution = 5

lake_level = 1767

xmin = 2669700.
xmax = 2671000.
ymin = 1170100.
ymax = 1172100.

dam_z = lake_level + 10
dam_thk = 30
limit = 2670561

out_format = "binary"

nx = 50
ny = int(1.62*nx)

amr_ratios = dict(
    x=(2,),
    y=(2,),
    t=(2,),
)

