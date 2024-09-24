xRes = yRes = 5

lake_level = 1767

xmin = 2669700.
xmax = 2671000.
ymin = 1170100.
ymax = 1172100.

dam_z = lake_level + 10
dam_thk = 30
limit = 2670561

out_format = "binary"

nx = 20
ny = int(1.62*nx)

amr_ratios = dict(
    x=(1,),
    y=(1,),
    t=(1,),
)

def dam_upstream(x, y, l=limit):
  yd = y.max() - 0.3*(x-x.min()) - 50000/(x-l) - 300
  yd[x > limit] = float("nan")
  return yd

def dam_downstream(x, y, thk=dam_thk):
  return dam_upstream(x+thk, y) + thk

